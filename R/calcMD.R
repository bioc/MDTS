calcMD = function(counts, bins, pData){
	normedCounts = normalizeCounts(counts, bins)

	proband_ind = which(pData$father_id %in% pData$subj_id)
	proband = pData$subj_id[proband_ind]

	father_ind = match(pData$father_id[proband_ind], pData$subj_id)
	mother_ind = match(pData$mother_id[proband_ind], pData$subj_id)

	md1 = normedCounts[,proband_ind] - normedCounts[,mother_ind]
	md2 = normedCounts[,proband_ind] - normedCounts[,father_ind]
	md = md1; md[abs(md1)>abs(md2)] = md2[abs(md1)>abs(md2)]

	colnames(md) = proband
	return(md)
}

