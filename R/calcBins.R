#' Sample driven calculation of bins for counting
#'
#' This function will randomly select a sample of bam files to calculate a range of bins to calculate read-depth on
#' @pData A table in the format of the output of pData()
#' @n The number of subsamples to use
#' @med The median number of reads across sub-samples to reach before creating a new bin
#' @min
#' @genome
#' @map_file
#' @seed Sets the seed so results are reproducible. Defaults to 1337
#' @keywords calcBins
#' @export
#' @examples
#' calcBins() 
calcBins <- function(pD, n, rl, med, min, genome, map_file, seed=1337){
	set.seed(seed)
	pD_sub = pD[sample(1:dim(pD)[1], n, replace=F),]
	
	print("Reading coverage information of subsamples"); flush.console()
	covs = list()
	for(i in 1:length(pD_sub$bam_path)){
		GR=BAM2GRanges(pD_sub$bam_path[i])
		reads = coverage(GR)
		reads = reads[names(reads) %in% 1:22]
		covs[[i]] = reads
	}

	print("Calculating Proto-regions"); flush.console()
	track = lapply(covs, function(x) which(x>=min))
	red = base::Reduce(intersect, track)

	bins = NULL
	for(chromosome in 1:22){
		bins = c(bins, processChr(as.character(chromosome), red, covs, rl, med))
	}
	bins_out = suppressWarnings(do.call('c', bins))
	print("Bin segmentation complete"); flush.console()

	if(sum(str_detect(seqlevels(bins_out), "chr"))==0){
		seqlevels(bins_out) = paste0("chr", seqlevels(bins_out))
	}
	seqs = getSeq(genome, bins_out)
	print("Calculating GC content"); flush.console()
	GC = sapply(seqs, calcGC)
	bins_out$GC = GC

	print("Calculating mappability"); flush.console()
	map_track = import(map_file)
	map = sapply(bins_out, calcMap, map_track)
	bins_out$mappability = map

	return(bins_out)
}
