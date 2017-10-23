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

	print("Reading coverage information of subsamples"); flush.console()
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

	print("Calculating Proto-regions"); flush.console()
	rle_threshold = sapply(covs, function(x) x>=min)
	rle_track = rle_threshold[[1]]
		for(i in 2:length(covs)){rle_track = rle_track + rle_threshold[[i]]}
	red = which(rle_track>0)

	bins = NULL
	for(chromosome in names(red)){
		bins = c(bins, processChr(as.character(chromosome), red, covs, rl, med))
	}
	bins_out = suppressWarnings(do.call('c', bins))
	print("Bin segmentation complete"); flush.console()

	seqlevels(bins_out) = paste0("chr", seqlevels(bins_out))
	seqs = getSeq(genome, bins_out)
	print("Calculating GC content"); flush.console()
	GC = sapply(seqs, calcGC)
	bins_out$GC = GC

	print("Calculating mappability"); flush.console()
	map_track = import(map_file)
	ol = findOverlaps(bins_out, map_track)
	segment_sums = by(rep(1, length(ol)), queryHits(ol), sum)
		segment_sums_1 = as.numeric(names(segment_sums)[which(segment_sums==1)])
	map = rep(0, length(segment_sums))
		map[segment_sums_1] = map_track$score[subjectHits(ol)[match(segment_sums_1, queryHits(ol))]]
	ol_remain = ol[!queryHits(ol) %in% segment_sums_1]
	map[unique(queryHits(ol_remain))] = sapply(unique(queryHits(ol_remain)), mappabilityHelper, ol_remain, bins_out, map_track)
	bins_out$mappability = map

	print("Filtering bins based on GC and Mappability"); flush.console()
	rm_ind = unique(c(which(bins_out$mappability<0.75), which(bins_out$GC>0.85 | bins_out$GC<0.15)))
	if(length(rm_ind)>0){bins_out = bins_out[-rm_ind]}

	if(hasChr==FALSE){
		seqlevels(bins_out) = str_replace(seqlevels(bins_out), "chr", "")
	}
	return(bins_out)
}
