#' Sample driven calculation of bins for counting
#'
#' This function will randomly select a sample of bam files to calculate dynamic
#' MDTS bins for subsequent read-depth analysis.
#' @param metaData A table in the format of the output of getMetaData().
#' @param n The number of subsamples to use.
#' @param readLength The read length of the experiment.
#' @param medianCoverage The median number of reads across sub-samples to reach 
#' before creating a new bin.
#' @param minimumCoverage The miminum number of coverage across all sub-samples 
#' required to create the proto-region.
#' @param genome The BSGenome object that assists in calculations of the GC 
#' content of the bins.
#' @param mappabilityFile A path to the bigwig file of 100mer mappability of the
#' corresponding genome.
#' @param seed Sets the seed so results are reproducible. Defaults to 1337.
#' @importFrom Rsamtools scanBamFlag
#' @importFrom  Biostrings getSeq
#' @importFrom GenomicAlignments readGAlignments
#' @importFrom GenomicAlignments coverage
#' @keywords calcBins
#' @examples 
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	bins
#' @export
#' @return Returns a \code{GRanges} object depicting the dynamic bins that MDTS 
#' calculates.
calcBins <- function(metaData, n, readLength, medianCoverage, minimumCoverage,
                     genome, mappabilityFile, seed=1337){
	set.seed(seed)
      metaSub <- metaData[sample(seq_len(n), replace=FALSE),]
	hasChr <- FALSE

	message("Reading coverage information of subsamples")
	covs <- vector('list', n)
	names(covs) <- metaSub[['bam_path']]
	for(i in metaSub[['bam_path']]){
	      reads = .extractCov(i)
		if(sum(stringr::str_detect(names(reads), "chr"))>0)
			hasChr=TRUE
		if(hasChr)
			names(reads) <- stringr::str_replace(names(reads), "chr", "")
		reads <- reads[names(reads) %in% 1:22]
		covs[[i]] <- reads
	}

	message("Calculating Proto-regions")
	rle_threshold <- lapply(covs, function(x) x>=minimumCoverage)
	rle_track <- rle_threshold[[1]]
		for(i in seq_len(n)[-1]){rle_track = rle_track + rle_threshold[[i]]}
	rle_track0 <- rle_track>0

	bins <- lapply(seq_len(22), .processChr, 
	              rle_track0, covs, readLength, medianCoverage)
	bins_out <- suppressWarnings(
	      do.call(c, bins[vapply(bins, length, integer(1))>0]))
	message("Bin segmentation complete")

      seqlevels(bins_out) <- paste0("chr", seqlevels(bins_out))
	seqs <- Biostrings::getSeq(genome, bins_out)
	
	message("Calculating GC content")
	GC <- vapply(seqs, .calcGC, double(1))
	bins_out$GC <- GC

	message("Calculating mappability")
	map_track <- import(mappabilityFile)
	ol <- findOverlaps(bins_out, map_track)
	segment_sums <- by(rep(1, length(ol)), queryHits(ol), sum)
	segment_sums_1 <- as.numeric(names(segment_sums)[which(segment_sums==1)])
	map <- rep(0, length(segment_sums))
	### For mappability=1 segments, no need to process
	map[segment_sums_1] <- map_track$score[subjectHits(ol)
            [match(segment_sums_1, queryHits(ol))]]
	### For mappability<1, process
	ol_remain = ol[!queryHits(ol) %in% segment_sums_1]
	map[unique(queryHits(ol_remain))] = vapply(unique(queryHits(ol_remain)), 
            .mappabilityHelper, double(1), ol_remain, bins_out, map_track)
	### Put the composite results together
	bins_out$mappability = map

	message("Filtering bins based on GC and Mappability")
	rm_ind <- unique(c(which(bins_out$mappability<0.75), 
            which(bins_out$GC>0.85 | bins_out$GC<0.15)))
	if(length(rm_ind)>0){bins_out = bins_out[-rm_ind]}

	if(!hasChr)
            seqlevels(bins_out) <- stringr::str_replace(seqlevels(bins_out), 
            "chr", "")
	return(bins_out)
}

## Helper functions
.extractCov <- function(path){
      message(path)
      flag <- Rsamtools::scanBamFlag(isUnmappedQuery = FALSE, 
                                     isDuplicate = FALSE)
      param <- Rsamtools::ScanBamParam(what = character(), flag=flag)      
      ga <- GenomicAlignments::readGAlignments(path, param=param)
      return(GenomicAlignments::coverage(ga))
}
.processChr = function(chr, proto_info, covs, rl, med){
      message(paste0("Selecting Proto-regions in Chr ", chr))
      proto_region <- proto_info[[chr]]
      if(sum(proto_region)>0){
            proto_gr <- GRanges(seqnames=chr, IRanges::IRanges(proto_region))
            proto_gr_covs_rle <- lapply(covs, function(x) x[proto_gr])
            proto_gr_covs <- lapply(proto_gr_covs_rle, function(x) 
                  lapply(x, sum))
            proto_gr_covs_mat <- apply(do.call(rbind, proto_gr_covs), 2, 
                  as.numeric)
            proto_gr_covs_mat_normed <- t(t(proto_gr_covs_mat)/rl)
            proto_gr_covs_mat_med <- apply(proto_gr_covs_mat_normed, 2, 
                  stats::median)
            proto_gr$reads <- proto_gr_covs_mat_med
            proto_gr_select <- proto_gr[proto_gr$reads>=med]
            
            if(length(proto_gr_select)>0){
                  message(paste0("Segmenting Chr ", chr, " Proto-regions"))
                  chr_out <- lapply(proto_gr_select, .divideSegs, covs, rl, med)
                  chr_out <- do.call(c, chr_out)
                  return(chr_out)
            }else{
                  return(NULL)
            }
      }else{
            return(NULL)
      }
}
.extractCounts <- function(cov_list){
      return(suppressWarnings((as.numeric(unlist(cov_list)))))
}
.divideSegs <- function(seg, covs, rl, med){
      cov_sub <- lapply(covs, function(x) x[seg])
      cov_sub_mat <- do.call(rbind, lapply(cov_sub, .extractCounts))
      cov_sub_cs_normed <- t(apply(cov_sub_mat, 1, cumsum))/rl
      medz <- apply(cov_sub_cs_normed, 2, stats::median)
      medzCut <- floor(medz/med)
      changePoint <- vapply(seq_len(max(medzCut)), function(x){
            min(which(medzCut==x))}, integer(1))
      changePoint <- c(0, changePoint)
            changePoint[length(changePoint)] <- width(seg)
      gr_out <- GRanges(seqnames=seqnames(seg),
            IRanges::IRanges(start(seg)+changePoint[-length(changePoint)],
            start(seg)+changePoint[-1]-1))
      gr_medz <- medz[changePoint]
      gr_out$median_count <- gr_medz-med*(seq_len(length(gr_medz))-1)
      return(gr_out)
}
.calcGC <- function(biostring){
      string <- as.character(biostring)
      return(stringr::str_count(string, "G|C")/nchar(string))
}
.mappabilityHelper <- function(i, ol, bins, map){
      overlaps <- pintersect(map[subjectHits(ol)[queryHits(ol)==i]], bins[i])
      percentOverlap <- width(overlaps) / width(bins[i])
      return(sum(percentOverlap*overlaps$score))
}


