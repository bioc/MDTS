#' Visualization for deletions
#'
#' This function plots the raw read information from the location of interest 
#' for a trio.
#' @param deletion A GRanges object in the format of the output of 
#' denovoDeletions().
#' @param bins The set of bins determined by calcBins().
#' @param pD A table in the format of the output of pData().
#' @param mCounts A matrix of normalized coverage output by normalizedCounts().
#' @param md The minimum distance matrix output by calcMD()
#' @param save If TRUE will save plot to current working directory instead of 
#' rendering.
#' @keywords visualizeDeletion
#' @examples 
#' \dontrun{
#'	load(system.file("extdata", 'bins.RData', package = "MDTS"))
#'	load(system.file("extdata", 'counts.RData', package = "MDTS"))
#'	load(system.file("extdata", 'pD.RData', package = "MDTS"))
#'	mCounts <- normalizeCounts(counts, bins)
#'	md <- calcMD(mCounts, pD)
#'	cbs <- segmentMD(md, bins)
#'	denovo <- denovoDeletions(cbs, mCounts, bins)
#'	visualizeDeletion(denovo[1], bins, pD, mCounts, md)
#'	}
#' @return The file name if the plot was saved.
#' @export
visualizeDeletion = function(deletion, bins, pD, mCounts, md, save=FALSE){
	window <- 1000
	famid <- pD$family_id[pD$subj_id==deletion$famid]
	pD_sub <- pD[stringr::str_detect(pD$family_id, famid),]
	row_inds <- subjectHits(findOverlaps(deletion, bins))
	col_inds <- match(pD_sub$subj_id, colnames(mCounts))
	
	id <- paste0(famid, "-", seqnames(deletion), ":", start(deletion), 
	            "-", end(deletion), ".pdf")
	if(save==TRUE){
		pdf(width=30, height=8, file=id)	
	}else{
		dev.new(width=30, height=8)	
	}
	.visualizeFamily(famid, deletion, window, row_inds, col_inds, pD, mCounts)
	dev.off()
	return(id)
}

### Helper functions
.pullReads = function(bam_path, bait){
      flag1 <- scanBamFlag(isPaired=TRUE, isFirstMateRead=TRUE, 
                          hasUnmappedMate=FALSE, isDuplicate=FALSE, 
                          isSecondaryAlignment=FALSE)
      flag2 <- scanBamFlag(isPaired=TRUE, isSecondMateRead=TRUE, 
                          hasUnmappedMate=FALSE, isDuplicate=FALSE, 
                          isSecondaryAlignment=FALSE)
      sbp <- ScanBamParam(flag=flag1, what=c("qname", "rname", "pos", "isize"), 
                       which=bait)
      bam1 <- scanBam(bam_path, index=paste0(bam_path, ".bai"), param=sbp)[[1]]
      bam2 <- scanBam(bam_path, index=paste0(bam_path, ".bai"), param=sbp)[[1]]
      
      counts <- by(rep(1, length(bam1$rname)+length(bam2$rname)), 
                  c(bam1$qname, bam2$qname), sum)
      purge <- names(counts)[which(as.vector(counts)==2)]
      ind <- which(bam1$qname %in% purge)
      GR1 <- GRanges(seqnames=as.character(seqnames(bait)[1]), 
                    IRanges(bam1$pos[ind], bam1$pos[ind]+99), 
                    isize=bam1$isize[ind],rn=bam1$qname[ind])
      x1 <- data.frame(ord=GR1$rn)
      x2 <- data.frame(ord=bam2$qname, 1:length(bam2$qname))
      x <- merge(x1, x2, sort=FALSE, all=FALSE, by="ord")
      ord <- x[,2]
      GR2 <- GRanges(seqnames=as.character(seqnames(bait)[1]), 
                    IRanges(bam2$pos[ord], bam2$pos[ord]+99), 
                    isize=bam2$isize[ord], rn=bam2$qname[ord])
      return(GRangesList(GR1, GR2))
}

.visualize = function(GRL, track, window, title, denote, coll, scale){
      GR1 <- GRL[[1]]
      GR2 <- GRL[[2]]
      del <- track
      start(del) <- start(del) + window
      end(del) <- end(del)- window
      frags <- countOverlaps(del, GR1)+ countOverlaps(del, GR2)
      coords <- c(start(GR1), start(GR2), end(GR1), end(GR2))
      plot(c(start(track), end(track)), c(0, length(GR1)+1), ty="n", 
           main=paste0(title, " M= ", denote), 
           ylab=yl, xaxt="n", cex.lab=scale, cex.axis=scale*.8, cex.main=scale, 
           cex.sub=scale, xlab="", col=coll, col.axis=coll, col.sub=coll, 
           col.main=coll, col.lab=coll)
      j <- 0.1
      tmp <- (start(GR1)<start(GR2))*1
      cls <- rep("#6495ed", length(tmp))
      cls[which(tmp==1)] <- "#ffb90f"
      rect(xleft=start(GR1), xright=start(GR1)+99, ytop = 1:length(GR1)-j, 
           ybot = 0:(length(GR1)-1)-j, col=cls, border=cls)
      rect(xleft=start(GR2), xright=start(GR2)+99, ytop = 1:length(GR2)-j, 
           ybot = 0:(length(GR2)-1)-j, col=cls, border=cls)
      segments(x0 = start(GR1)+99, x1=start(GR2), y0 = 1:length(GR1)-j-0.5, 
               y1 = 1:length(GR2)-j-0.5, lwd=0.1, col= cls)
      
      z <- width(track)-2*window
      text(x=(start(track)+end(track))/2, y=length(GR1)-length(GR1)/20, 
           labels=paste0(z, "bp"), cex=scale, col=coll)
      text(x=(start(track)+end(track))/2, length(GR1)/20, 
           labels=paste0(frags, " reads"), cex=scale, col=coll)
      ylims <- par("usr")[3:4]
      rect(xleft=start(track)+window, xright=end(track)-window, ytop = ylims[2],
           ybot=ylims[1], col=rgb(0, 0, 0, 0.2), border=rgb(0,0,0,0)) 
}
.visualizeFamily = function(famid, bait, window, row_inds, col_inds, pD, 
            mCounts, scale=4, col1 = rgb(0, 146, 146, max=255), 
            col2 = rgb(255, 109, 182, max=255)){
      sub <- pD[pD$family_id==famid,]
      label <- paste0(seqnames(bait), ":", start(bait), "-", end(bait))
      bamids <- as.character(sub$bamid)
      bamids <- sub$bam_path
      start(bait) <- start(bait)-window; end(bait) = end(bait)+window
      
      cands1 <- .pullReads(bamids[1], bait)
      cands2 <- .pullReads(bamids[2], bait)
      cands3 <- .pullReads(bamids[3], bait)
      
      res <- mCounts[row_inds, ,drop=FALSE]
      i1 <- which(colnames(res)==sub$subj_id[1])
      i2 <- which(colnames(res)==sub$subj_id[2])
      i3 <- which(colnames(res)==sub$subj_id[3])
      md1 <- res[,i1,drop=FALSE]-res[,i2,drop=FALSE]
      md2 <- res[,i1,drop=FALSE]-res[,i3,drop=FALSE]; md = md1
      md[abs(md1)>abs(md2)] <- md2[abs(md1)>abs(md2)]
      md <- mean(md); md = round(md, 3)
      res <- apply(res, 2, mean)
      res <- round(res, 3)
      
      layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, rep(10, 3)), 2, 6, byrow=TRUE),
             widths=c(0.5, 0.5, rep(5, 4)), heights = c(4, 0.1))
      par(mar=rep(0,4))
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", 
           type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=0.5, srt=90, paste0(famid, " - MD = ", md), cex=4, pos=3)
      plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", 
           type="n", main="", xaxs="i", yaxs="i", bty="n")
      text(x=0.6, y=.5, srt=90, label, cex=4, pos=3)
      par(xpd=NA)
      text(x=1.2, y=.5, srt=90, paste0(length(row_inds), " bins"), pos=3, cex=4)
      par(xpd=FALSE)
      par(mar=c(scale, scale*3, scale, scale))
      minM <- stats::quantile(res, 0.3)-1
      maxM <- stats::quantile(res, 0.97)+1
      plot_wid = (maxM-minM)
      hist(res, main="", xlim=c(minM, maxM), breaks=seq(-10, 10, by=0.3), 
           ylab="Hist M Scores", col="grey", cex.lab=scale*.8, 
           cex.axis=scale*.7, cex.main=scale, cex.sub=scale, xlab="")
      plot_lims <- par("usr")
      ylim <- plot_lims[3:4]
      lens <- (ylim[2]-ylim[1])/30
      ynotches <- seq(ylim[1], ylim[2], by = lens)
      rect(xleft=res[i1]-plot_wid/100, xright=res[i1]+plot_wid/100, 
           ytop=ynotches[seq(1, length(ynotches), by=3)]+lens, 
           ybot=ynotches[seq(1, length(ynotches), by=3)], col=1)
      rect(xleft=res[i2]-plot_wid/100, xright=res[i2]+plot_wid/100, 
           ytop=ynotches[seq(2, length(ynotches), by=3)]+lens, 
           ybot=ynotches[seq(2, length(ynotches), by=3)], col=col1, bord=col1)
      rect(xleft=res[i3]-plot_wid/100, xright=res[i3]+plot_wid/100, 
           ytop=ynotches[seq(3, length(ynotches), by=3)]+lens, 
           ybot=ynotches[seq(3, length(ynotches), by=3)], col=col2, bord=col2)
      
      par(xpd=FALSE)
      par(mar=rep(5,4))
      .visualize(cands1, bait, window, "Proband", 
                 denote=res[i1], coll=1, scale=scale)
      .visualize(cands2, bait, window, "Parent 1", 
                 denote=res[i2], coll=col1, scale=scale)
      .visualize(cands3, bait, window, "Parent 2", 
                 denote=res[i3], coll=col2, scale=scale)	
      
      par(mar=rep(0,4))
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      plot(c(0,1), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      par(xpd=NA)
      plot(c(0,3), c(0,1), type="n", axes=FALSE, xlab="", ylab="")
      rect(xleft=0.6, xright=0.7, ybot=1.2, ytop=1.4, 
           bord="#ffb90f", col="#ffb90f")
      rect(xleft=1.7, xright=1.8, ybot=1.2, ytop=1.4, 
           bord="#6495ed", col="#6495ed")
      text(x=1, y=1.3, labels="+ Strand Reads", cex=scale*.8)
      text(x=2.1, y=1.3, labels="- Strand Reads", cex=scale*.8)
      par(xpd=FALSE)
}
