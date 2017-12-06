#' Denovo Deletion Calling
#'
#' This function will return a single GRanges object containing all denovo deletions 
#' 		that passed filtering from a Circular Binary Segmentation object with supplementary information.
#' @param cbs The GRanges output from segmentMD() that contains the copy number segments of the families.
#' @param mCounts The normalized counts matrix output by normalizeCounts().
#' @param bins The set of bins determined by calcBins().
#' @keywords denovoDeletions
#' @examples 
#'	setwd(system.file('extdata', package='MDTS'))
#'	load('bins.RData')
#'	load('counts.RData')
#'	load('pD.RData')
#'	mCounts = normalizeCounts(counts, bins)
#'	md = calcMD(mCounts, bins, pD)
#'	cbs = segmentMD(md, bins)
#'	denovo = denovoDeletions(cbs, mCounts, bins)
#' @export
denovoDeletions = function(cbs, mCounts, bins){
	dels = GRanges()

	print("Selecting Candidate de Novo deletions")
	candidate = cbs[abs(cbs$m+1)<0.3]
	candidate_fam_count = by(rep(1, length(candidate)), candidate$famid, sum)
	mCounts_local = mCounts

	print("Calculating problematic bins")
	win = 0.5
	raw = abs(mCounts_local)<win
	raw_perc = apply(raw, 1, mean)
	raw_perc = sapply(raw_perc+1/dim(mCounts)[2], min, 1)

	# ### New version
	# ol = findOverlaps(candidate, bins)
	# percs = by(subjectHits(ol), queryHits(ol), function(x) mean(raw_perc[x]))
	# keep = which(percs>=0.95)
	# dels = candidate[keep]

	# ### Old version
	bins_filter = bins[which(raw_perc<0.95)]

	print("Filtering candidates by problematic bins")
	ol_filter = findOverlaps(candidate, bins_filter)
	ol_bins = findOverlaps(candidate, bins)
	if(length(ol_filter)>0){
		count_filter = by(rep(1, length(ol_filter)), queryHits(ol_filter), sum)
		count_bins = by(rep(1, length(ol_bins)), queryHits(ol_bins), sum)
		filtering_info = data.frame(id = names(count_bins), count_base = as.numeric(count_bins))
		filtering_info2 = data.frame(id = names(count_filter), count_filter = as.numeric(count_filter))
		filtering = merge(filtering_info, filtering_info2, by="id", sort=F)
		ratios = filtering[,3]/filtering[,2]
		cut = which(ratios>=0.5)
		drop_ids = as.numeric(as.character(filtering$id[cut]))
		if(length(drop_ids)>0){
		      dels = candidate[-drop_ids]
		}else{
		      dels = candidate
		}
	}else{
		dels = candidate
	}
	return(dels)
}