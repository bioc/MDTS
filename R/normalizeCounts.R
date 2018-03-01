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
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	load(system.file("extdata", 'pD.RData', package = "MDTS"))
#'	mCounts <- normalizeCounts(counts, bins)
#' @export
#' @return A \code{data.frame} of normalized counts. Each column is a sample,
#' and each row is a entry of \code{bins}.
normalizeCounts <- function(counts, bins, GC=TRUE, map=TRUE, mc.cores=1){
	message("Log Transforming Counts")
	log_counts <- log(counts+1, 2)
	intermediate <- t(t(log_counts) - apply(log_counts, 2, stats::median))
	res <- intermediate - apply(intermediate, 1, stats::median)
	
	if(GC){
	      message("GC Adjust")
	      res <- do.call(cbind, mclapply(seq_len(dim(res)[2]), function(i){
	            control <- stats::loess.control(trace.hat="approximate")
                  res <- stats::residuals(stats::loess(res[,i]~bins$GC, 
                        control = control))
	      }, mc.cores=mc.cores))
	}
	if(map){
	      message("Mappability Adjust")
            loess_ind <- which(bins$mappability<1)
	      res <- do.call(cbind, mclapply(seq_len(dim(res)[2]), function(i){
                  tmp <- res[,i]
                  tmp[loess_ind] <- stats::residuals(stats::loess(
                        tmp[loess_ind]~bins$mappability[loess_ind]))
                  return(tmp)
	      }, mc.cores=mc.cores))
	}
	colnames(res) <- colnames(counts)
	return(res)
}
