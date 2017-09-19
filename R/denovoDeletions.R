#' This function will report the denovo deletions based on the copy number segment GRanges object with appropriate filters
#'
#' This function will return a single GRanges object containing all denovo deletions that passed filtering
#' @param cbs The GRanges output from segmentMD() that contains the copy number segments of the families
#' @param mCounts The normalized counts matrix output by normalizeCounts()
#' @param bins The set of bins determined by calcBins
#' @keywords denovoDeletions
#' @export
denovoDeletions = function(cbs, mCounts, bins){
	print("Selecting Ccndidate deNnvo deletions")
	candidate = cbs[abs(cbs$m+1)<0.3]
	candidate_fam_count = by(rep(1, length(candidate)), candidate$famid, sum)

	print("Filtering out families with probable sequencing failure")
	bad_families = names(candidate_fam_count)[candidate_fam_count>4]
	if(length(bad_families)>0){
		mCounts = mCounts[-str_replace(colnames(mCounts), "_[0-9]*", "") %in% bad_families]
	}
	candidate = candidate[!candidate$famid %in% bad_families]

	print("Calculating problematic bins")
	win = 0.5
	raw = abs(mCounts)<win
	bins_filter = bins[which(apply(raw, 1, mean)<0.95)]

	print("Filtering candidates by problematic bins")
	ol_filter = findOverlaps(candidate, bins_filter)
	ol_bins = findOverlaps(candidate, bins)
	count_filter = by(rep(1, length(ol_filter)), queryHits(ol_filter), sum)
	count_bins = by(rep(1, length(ol_bins)), queryHits(ol_bins), sum)
	filtering_info = data.frame(id = names(count_bins), count_base = as.numeric(count_bins))
	filtering_info2 = data.frame(id = names(count_filter), count_filter = as.numeric(count_filter))
	filtering = merge(filtering_info, filtering_info2, by="id", sort=F)
	ratios = filtering[,3]/filtering[,2]
	cut = which(ratios>=0.5)
	drop_ids = as.numeric(as.character(filtering$id[cut]))
	
	dels = candidate[-drop_ids]
	return(dels)
}