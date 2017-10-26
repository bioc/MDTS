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
# 	aligns <- readGAlignments(path, paste0(path, ".bai"), param = filters)
# 	aligns.info <- values(aligns)
# 	  aligns.gr <- as(aligns, "GRanges")
# 	values(aligns.gr) <- aligns.info
# 	return(aligns.gr)
# }

# bv = BamViews(pD$bam_path[1:5], bamRanges=bins[1])
# test = countBam(bv)

# 		  what = character()
# 		  flag = scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE)
#           filters <- ScanBamParam(what = what, flag = flag)
# 	      aligns <- readGAlignments(path, param = filters)
# 	      aligns.info <- values(aligns)
#           aligns.gr <- as(aligns, "GRanges")
# 	      values(aligns.gr) <- aligns.info
# 	      aligns.gr
