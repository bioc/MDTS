#' Calculating the normalized M scores
#'
#' This function will return a matrix of normalized M scores where ecah column 
#' is a sample, and each row is a bin.
#' @param counts A matrix of raw coverage output by calcCounts().
#' @param bins The set of bins determined by calcBins().
#' @param mc.cores The number of cores to use for multi-threaded analysis. 
#' @param GC Whether to loess adjust for GC. Defaults to TRUE.
#' @param map Whether to loess adjust for mappability. Defaults to TRUE.
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
normalizeCounts = function(counts, bins, GC=TRUE, map=TRUE, mc.cores=1){
	message("Log Transforming Counts")
	log_counts = log(counts+1, 2)
	intermediate = t(t(log_counts) - apply(log_counts, 2, stats::median))
	res = intermediate - apply(intermediate, 1, stats::median)
	if(GC==TRUE){
	      message("GC Adjust")
	      res = do.call(cbind, mclapply(1:dim(counts)[2], .fitLoess, bins, res, "GC", mc.cores=mc.cores))      
	}
	if(map==TRUE){
	      message("Mappability Adjust")
	      res = do.call(cbind, mclapply(1:dim(counts)[2], .fitLoess, bins, res, "Mappability", mc.cores=mc.cores))
	}
	colnames(res) = colnames(counts)
	return(res)
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

