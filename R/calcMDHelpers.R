normalizeCounts = function(countMat, bins){
	print("Log Transforming Counts"); flush.console()
	log_counts = log(countMat+1, 2)
	intermediate = t(t(log_counts) - apply(log_counts, 2, median))
	res = intermediate - apply(intermediate, 1, median)

	print("GC Adjust"); flush.console()
	res_gc = do.call(cbind,mclapply(1:dim(countMat)[2], fitLoess, bins, res, "GC", mc.cores=detectCores()))
	print("Mappability Adjust"); flush.console()
	res_gc_map = do.call(cbind, mclapply(1:dim(countMat)[2], fitLoess, bins, res_gc, "Mappability", mc.cores=detectCores()))
	colnames(res_gc_map) = colnames(countMat)
	return(res_gc_map)
}

fitLoess = function(i, bins, full_data, adjust){
	if(adjust=="GC"){
		res = residuals(loess(full_data[,i]~bins$GC))
	}
	if(adjust=="Mappability"){
		loess_ind = which(bins$mappability<1)
		res = full_data[,i]
		res[loess_ind] = residuals(loess(full_data[loess_ind,i]~bins$mappability[loess_ind]))
	}
	return(res)
}


