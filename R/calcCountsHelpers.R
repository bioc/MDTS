countHelper = function(path, bins){
	GR=BAM2GRanges(path)
	reads = coverage(GR)
	reads = reads[match(seqlevels(bins), names(reads))]
	binned = binnedAverage(bins, reads, "count")
	return(binned$count)
}
