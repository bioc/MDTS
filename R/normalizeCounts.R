#' Calculating the normalized M scores
#'
#' This function will return a matrix of normalized M scores where ecah column 
#' is a sample, and each row is a bin.
#' @param counts A matrix of raw coverage output by calcCounts().
#' @param bins The set of bins determined by calcBins().
#' @param mc.cores The number of cores to use for multi-threaded analysis. 
#' Defaults to 1.
#' @keywords normalizeCounts
#' @examples 
#'	setwd(system.file('extdata', package='MDTS'))
#'	load('bins.RData')
#'	load('counts.RData')
#'	load('pD.RData')
#'	mCounts = normalizeCounts(counts, bins)
#' @export
#' @return A \code{data.frame} of normalized counts. Each column is a sample,
#' and each row is a entry of \code{bins}.
normalizeCounts = function(counts, bins, mc.cores=1){
	message("Log Transforming Counts")
	log_counts = log(counts+1, 2)
	intermediate = t(t(log_counts) - apply(log_counts, 2, median))
	res = intermediate - apply(intermediate, 1, median)

	message("GC Adjust")
	res_gc = do.call(cbind, mclapply(1:dim(counts)[2], .fitLoess, bins, res, "GC", mc.cores=mc.cores))
	
	message("Mappability Adjust")
	res_gc_map = do.call(cbind, mclapply(1:dim(counts)[2], .fitLoess, bins, res_gc, "Mappability", mc.cores=mc.cores))
	
	colnames(res_gc_map) = colnames(counts)
	return(res_gc_map)
}

## Helper functions
.fitLoess = function(i, bins, full_data, adjust){
      if(adjust=="GC"){
            control = stats::loess.control(trace.hat="approximate")
            res = stats::residuals(stats::loess(full_data[,i]~bins$GC, control = control))
      }
      if(adjust=="Mappability"){
            loess_ind = which(bins$mappability<1)
            res = full_data[,i]
            res[loess_ind] = stats::residuals(stats::loess(full_data[loess_ind,i]~bins$mappability[loess_ind]))
      }
      return(res)
}

