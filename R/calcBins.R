#' Sample driven calculation of bins for counting
#'
#' This function will randomly select a sample of bam files to calculate a range of bins to calculate read-depth on
#' @pData A table in the format of the output of pData()
#' @n The number of subsamples to use
#' @med The median number of reads across sub-samples to reach before creating a new bin
#' @seed Sets the seed so results are reproducible. Defaults to 1337
#' @keywords calcBins
#' @export
#' @examples
#' calcBins() 
calcBins <- function(pData, n, rl, med, min, seed=1337){
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
	red = Reduce(intersect, track)

	bins = NULL
	for(chromosome in 1:22){
		bins = c(bins, processChr(chromosome, red, covs, rl, med))
	}
	bins_out = suppressWarnings(do.call('c', bins))
}

