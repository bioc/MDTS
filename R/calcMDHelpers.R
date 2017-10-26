fitLoess = function(i, bins, full_data, adjust){
	if(adjust=="GC"){
		control = loess.control(trace.hat="approximate")
		res = residuals(loess(full_data[,i]~bins$GC, control = control))
	}
	if(adjust=="Mappability"){
		loess_ind = which(bins$mappability<1)
		res = full_data[,i]
		res[loess_ind] = residuals(loess(full_data[loess_ind,i]~bins$mappability[loess_ind]))
	}
	return(res)
}


