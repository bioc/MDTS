countHelper = function(path, bins){
	GR=BAM2GRanges(path)
	reads = coverage(GR)
	reads = reads[match(seqlevels(bins), names(reads))]
	binned = binnedAverage(bins, reads, "count")
	return(binned$count)
}

# BAM2GRanges2 = function(path, bins){
# 	what = character()
# 	flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE)
# 	filters <- ScanBamParam(what = what, flag = flag, which=bins)	
# 	aligns <- readGAlignments(pD$bam_path[1], param = filters)
# 	aligns.info <- values(aligns)
# 	  aligns.gr <- as(aligns, "GRanges")
# 	values(aligns.gr) <- aligns.info
# 	return(aligns.gr)
# }
