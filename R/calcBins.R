#' Sample driven calculation of bins for counting
#'
#' This function will randomly select a sample of bam files to calculate a range of bins to calculate read-depth on
#' @param pD A table in the format of the output of pData()
#' @param n The number of subsamples to use
#' @param rl The read length of the experiment
#' @param med The median number of reads across sub-samples to reach before creating a new bin
#' @param min The miminum number of coverage across all sub-samples required to create the proto-region
#' @param genome The BSGenome object that assists in calculations of the GC content of the bins
#' @param map_file A path to the bigwig file of 100mer mappability of the corresponding genome
#' @param seed Sets the seed so results are reproducible. Defaults to 1337
#' @keywords calcBins
#' @export
calcBins <- function(pD, n, rl, med, min, genome, map_file, seed=1337){
	set.seed(seed)
	pD_sub = pD[sample(1:dim(pD)[1], n, replace=F),]
	hasChr = FALSE

	message("Reading coverage information of subsamples")
	covs = list()
	for(i in 1:length(pD_sub$bam_path)){
		GR=BAM2GRanges(pD_sub$bam_path[i])
		reads = coverage(GR)
		if(sum(str_detect(names(reads), "chr"))>0){
			hasChr=TRUE
		}
		if(hasChr==TRUE){
			names(reads) = str_replace(names(reads), "chr", "")
		}
		reads = reads[names(reads) %in% 1:22]
		covs[[i]] = reads
	}

	message("Calculating Proto-regions")
	rle_threshold = sapply(covs, function(x) x>=min)
	rle_track = rle_threshold[[1]]
		for(i in 2:length(covs)){rle_track = rle_track + rle_threshold[[i]]}
	red = which(rle_track>0)
	      red = red[which(sapply(red, length)>0)]

	bins = NULL
	for(chromosome in names(red)){
		bins = c(bins, .processChr(as.character(chromosome), red, covs, rl, med))
	}
	bins_out = suppressWarnings(do.call('c', bins))
	message("Bin segmentation complete")

	seqlevels(bins_out) = paste0("chr", seqlevels(bins_out))
	seqs = getSeq(genome, bins_out)
	
	message("Calculating GC content")
	GC = sapply(seqs, .calcGC)
	bins_out$GC = GC

	message("Calculating mappability")
	map_track = import(map_file)
	ol = findOverlaps(bins_out, map_track)
	segment_sums = by(rep(1, length(ol)), queryHits(ol), sum)
		segment_sums_1 = as.numeric(names(segment_sums)[which(segment_sums==1)])
	map = rep(0, length(segment_sums))
		map[segment_sums_1] = map_track$score[subjectHits(ol)[match(segment_sums_1, queryHits(ol))]]
	ol_remain = ol[!queryHits(ol) %in% segment_sums_1]
	map[unique(queryHits(ol_remain))] = sapply(unique(queryHits(ol_remain)), .mappabilityHelper, ol_remain, bins_out, map_track)
	bins_out$mappability = map

	message("Filtering bins based on GC and Mappability")
	rm_ind = unique(c(which(bins_out$mappability<0.75), which(bins_out$GC>0.85 | bins_out$GC<0.15)))
	if(length(rm_ind)>0){bins_out = bins_out[-rm_ind]}

	if(hasChr==FALSE)
		seqlevels(bins_out) = str_replace(seqlevels(bins_out), "chr", "")
	
	return(bins_out)
}

## Helper functions
.processChr = function(chr, proto_info, covs, rl, med){
      message(paste0("Selecting Proto-regions in Chr ", chr))
      proto_region = proto_info[[chr]]
      proto_gr = reduce(GRanges(seqnames=chr, IRanges(start=proto_region, end=proto_region)))
      
      proto_gr_covs_rle = lapply(covs, function(x) x[proto_gr])
      proto_gr_covs = lapply(proto_gr_covs_rle, function(x) lapply(x, sum))
      proto_gr_covs_mat = apply(do.call(rbind, proto_gr_covs), 2, as.numeric)
      proto_gr_covs_mat_normed = t(t(proto_gr_covs_mat)/rl)
      proto_gr_covs_mat_med = apply(proto_gr_covs_mat_normed, 2, median)
      proto_gr$reads = proto_gr_covs_mat_med
      proto_gr_select = proto_gr[proto_gr$reads>=med]
      
      if(length(proto_gr_select)>0){
            message(paste0("Segmenting Chr ", chr, " Proto-regions"))
            pb = txtProgressBar(min = 0, max = length(proto_gr_select), style = 3)
            chr_out = NULL
            for(i in 1:length(proto_gr_select)){
                  setTxtProgressBar(pb, i)
                  chr_out = c(chr_out, .divideSegs(proto_gr_select[i], covs, rl, med))
            }
            chr_out = do.call('c', chr_out)
      }else{
            return(NULL)
      }
      return(chr_out)
}
.extractCounts = function(cov_list){
      return(suppressWarnings((as.numeric(unlist(cov_list)))))
}
.divideSegs = function(seg, covs, rl, med){
      output = NULL
      j = 1
      count = 1
      reads_remain = seg$read
      num_segs = floor(seg$reads/med)
      cov_sub = lapply(covs, function(x) x[seg])
      cov_sub_mat = do.call(rbind, lapply(cov_sub, .extractCounts))
      
      cov_sub_cs_normed = t(apply(cov_sub_mat, 1, cumsum))/rl
      medz = apply(cov_sub_cs_normed, 2, median)
      while(max(medz)>=med){
            cut = which(medz>=med)[1]
            output = rbind(output, c(cut, medz[cut]))
            count = count+1; j = sum(output[,1])+1
            cov_sub_cs = t(apply(cov_sub_mat[,j:dim(cov_sub_mat)[2]], 1, cumsum))
            cov_sub_cs_normed = cov_sub_cs/rl
            medz = apply(cov_sub_cs_normed, 2, median)
      }
      output = rbind(output, c(width(seg)-sum(output[,1]), max(medz)))

      cs = cumsum(output[,1])
      coords = cbind(c(0, cs[-length(cs)]), cs-1)
      gr_out = GRanges(seqnames=seqnames(seg), IRanges(start(seg)+coords[,1], start(seg)+coords[,2]))
      gr_out$median_count = output[,2]
      return(gr_out)
}
.calcGC = function(biostring){
      string = as.character(biostring)
      return(str_count(string, "G|C")/nchar(string))
}
.mappabilityHelper = function(i, ol, bins, map){
      overlaps <- pintersect(map[subjectHits(ol)[queryHits(ol)==i]], bins[i])
      percentOverlap <- width(overlaps) / width(bins[i])
      return(sum(percentOverlap*overlaps$score))
}
