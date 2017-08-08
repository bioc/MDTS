processChr = function(chr, proto_info, covs, rl, med){
	print(paste0("Selecting Proto-regions in Chr ", chr)); flush.console()
	proto_region = proto_info[[chr]]
	proto_gr = reduce(GRanges(seqnames=chr, IRanges(start=proto_region, end=proto_region)))
	proto_gr_covs_rle = lapply(covs, function(x) x[proto_gr])
	proto_gr_covs = lapply(proto_gr_covs_rle, function(x) lapply(x, sum))
	proto_gr_covs_mat = apply(do.call(rbind, proto_gr_covs), 2, as.numeric)
	proto_gr_covs_mat_normed = t(t(proto_gr_covs_mat)/rl)
	proto_gr_covs_mat_med = apply(proto_gr_covs_mat_normed, 2, median)
		proto_gr$reads = proto_gr_covs_mat_med
	proto_gr_select = proto_gr[proto_gr$reads>=med]

	if(length(proto_gr_select)>0){
		print(paste0("Segmenting Chr ", chr, " Proto-regions")); flush.console()
		pb = txtProgressBar(min = 0, max = length(proto_gr_select), style = 3)
		chr_out = NULL
		for(i in 1:length(proto_gr_select)){
			setTxtProgressBar(pb, i)
			chr_out = c(chr_out, divideSegs(proto_gr_select[i], covs, rl, med))
		}
		print(""); flush.console()
		chr_out = do.call('c', chr_out)
	}else{
		return(NULL)
	}
	return(chr_out)
}

extractCounts = function(cov_list){
	return(suppressWarnings((as.numeric(unlist(cov_list)))))
}

divideSegs = function(seg, covs, rl, med){
	output = NULL
	j = 1
	count = 1
	reads_remain = seg$read
	num_segs = floor(seg$reads/med)
	cov_sub = lapply(covs, function(x) x[seg])
	cov_sub_mat = do.call(rbind, lapply(cov_sub, extractCounts))
	while(count<num_segs){
		print(paste0(count, "-", num_segs)); flush.console()
		cov_sub_cs = t(apply(cov_sub_mat[,j:dim(cov_sub_mat)[2]], 1, cumsum))
		cov_sub_cs_normed = cov_sub_cs/rl
		medz = apply(cov_sub_cs_normed, 2, median)
		cut = which(medz>=med)[1]
		if(!is.na(cut)){
			output = rbind(output, c(cut, medz[cut]))
		}
		count = count+1; j = sum(output[,1])+1
	}
	cov_sub_cs = t(apply(cov_sub_mat[,j:dim(cov_sub_mat)[2]], 1, cumsum))
	cov_sub_cs_normed = cov_sub_cs/rl
	medz = apply(cov_sub_cs_normed, 2, median)
	output = rbind(output, c(width(seg)-sum(output[,1]), max(medz)))
	cs = cumsum(output[,1])
	coords = cbind(c(0, cs[-length(cs)]), cs-1)
	gr_out = GRanges(seqnames=seqnames(seg), IRanges(start(seg)+coords[,1], start(seg)+coords[,2]))
		gr_out$median_count = output[,2]
	return(gr_out)
}


calcGC = function(biostring){
	string = as.character(biostring)
	return(str_count(string, "G|C")/nchar(string))
}

calcMap = function(bins, map_track){
	ol = findOverlaps(bins, map_track)
	segment_counts = by(rep(1, length(ol)), queryHits(ol), sum)
	map = rep(0, length(bins))
	map[as.numeric(names(segment_counts)[which(segment_counts==1)])] =
		 map_track$score[subjectHits(ol)[which(segment_counts==1)]]
	
	ind_compute = names(segment_counts)[which(segment_counts>1)]
	map[as.numeric(names(segment_counts)[which(segment_counts>1)])] =
		 sapply(bins[as.numeric(ind_compute)], calcMapHelper, map_track)
	return(map)
}

calcMapHelper = function(bin, map_track){
	ol = findOverlaps(bin, map_track)
	map_track_sub = map_track[subjectHits(ol)]
	map_track_red = do.call(c, lapply(map_track_sub[c(1, length(ol))], intersect, bin))
	map_track_red$score = map_track_sub$score[c(1, length(ol))]
	map_track_sub[c(1, length(ol))] = map_track_red
	return(sum(width(map_track_sub)*map_track_sub$score)/width(bin))
}