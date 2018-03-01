#' Circular Binary Segmentation on Minimum Distances
#'
#' This function will return a GRanges object containing the copy number 
#'	segments of all families in the input minimum distance matrix. It calls
#'	segment() from DNAcopy (alpha=0.001, undo.splits="sdundo", undo.SD=4).
#' @param md The minimum distance matrix produced by calcMD.
#' @param bins The set of bins determined by calcBins.
#' @param alpha Controls the alpha option in calling DNAcopy::segment()
#' @param undo.splits Controls the undo.splits option in DNAcopy::segment()
#' @param undo.SD Controls the undo.SD option in calling DNAcopy::segment()
#' @param mc.cores The number of cores to use for multi-threaded analysis. 
#' Defaults to 1.
#' @keywords segmentMD
#' @examples 
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	load(system.file("extdata", 'pD.RData', package = "MDTS"))
#'	mCounts <- normalizeCounts(counts, bins)
#'	md <- calcMD(mCounts, metaData)
#'	cbs <- segmentMD(md, bins)
#' @export
#' @return A \code{data.frame} containing the segmented regions based to be
#' parsed by denovoDeletions()
#' minimum distance.
segmentMD <- function(md, bins, alpha=0.001, undo.splits='sdundo', undo.SD=4, 
                     mc.cores=1){
	segs <- mclapply(1:(dim(md)[2]), .segmentMDHelper, 
            md=md, bins=bins, alpha=alpha, 
            undo.splits = undo.splits, undo.SD=undo.SD, 
            mc.cores=mc.cores)
	segs.out <- do.call(rbind, segs)
	return(segs.out)
}

## Helper function
.segmentMDHelper <- function(i, md, bins, alpha, undo.splits, undo.SD){
      set.seed(137)
      message(paste0("Processing family number: ", i))
      md_sub <- md[,i,drop=FALSE]
      cna <- CNA(genomdat=md_sub, chrom=as.vector(seqnames(bins)),
                maploc=start(bins), data.type="logratio", 
                sampleid=colnames(md_sub), presorted=TRUE)
      cbs <- segment(cna, alpha=alpha, undo.splits=undo.splits, undo.SD=undo.SD)
      
      segRows <- cbs$segRows
      out <- data.frame(start=segRows[,1], end=segRows[,2], 
                       seg.mean=cbs$output$seg.mean, family=colnames(md_sub))
      return(out)
}
