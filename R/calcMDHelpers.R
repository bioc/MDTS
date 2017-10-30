fitLoess = function(i, bins, full_data, adjust){
	print(i); flush.console()
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

segmentMDCore = function(i, md, bins, alpha, undo.splits, undo.SD){
	print(paste0("Processing family number: ", i)); flush.console()
	cna = CNA(genomdat=md[,i,drop=FALSE], chrom=as.vector(seqnames(bins)),
		maploc=start(bins), data.type="logratio", sampleid=colnames(md), presorted=T)
	cbs = segment(cna, alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD, verbose=0)

	segRows = cbs$segRows
	gr = GRanges(seqnames=cbs$output$chrom, 
		IRanges(start(bins)[segRows[,1]], end(bins)[segRows[,2]]), 
		m=cbs$output$seg.mean, 
		famid=colnames(md)[i])
	return(gr)
}