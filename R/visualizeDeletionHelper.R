pullReads = function(bam_path, bait){
	target = bait
	flag1 = scanBamFlag(isPaired=T, isFirstMateRead=T, hasUnmappedMate=F, isDuplicate=F, isSecondaryAlignment=F)
	flag2 = scanBamFlag(isPaired=T, isSecondMateRead=T, hasUnmappedMate=F, isDuplicate=F, isSecondaryAlignment=F)
	bam1 = scanBam(bam_path, index=paste0(bam_path, ".bai"), param=ScanBamParam(flag=flag1, what=c("qname", "rname", "pos", "isize"), which=target))[[1]]
	bam2 = scanBam(bam_path, index=paste0(bam_path, ".bai"), param=ScanBamParam(flag=flag2, what=c("qname", "rname", "pos", "isize"), which=target))[[1]]
	
	counts = by(rep(1, length(bam1$rname)+length(bam2$rname)), c(bam1$qname, bam2$qname), sum)
	purge = names(counts)[which(as.vector(counts)==2)]
	ind = which(bam1$qname %in% purge)
	GR1 = GRanges(seqnames=as.character(seqnames(bait)[1]), IRanges(bam1$pos[ind], bam1$pos[ind]+99), isize=bam1$isize[ind],rn=bam1$qname[ind])
	x1 = data.frame(ord=GR1$rn)
	x2 = data.frame(ord=bam2$qname, 1:length(bam2$qname))
	x = merge(x1, x2, sort=F, all=F, by="ord")
	ord=x[,2]
	GR2 = GRanges(seqnames=as.character(seqnames(bait)[1]), IRanges(bam2$pos[ord], bam2$pos[ord]+99), isize=bam2$isize[ord], rn=bam2$qname[ord])
	return(GRangesList(GR1, GR2))
}
visualize = function(GRL, track, window, title, yl, denote, coll, bait, scale){
	GR1 = GRL[[1]]
	GR2 = GRL[[2]]
	coords = c(start(GR1), start(GR2), end(GR1), end(GR2))
	plot(c(start(bait), end(bait)), c(0, length(GR1)+1), ty="n", main=title, ylab=yl, xaxt="n", cex.lab=scale, cex.axis=scale*.8, cex.main=scale, cex.sub=scale, xlab="", col=coll, col.axis=coll, col.sub=coll, col.main=coll, col.lab=coll)
	j=0.1
	tmp = (start(GR1)<start(GR2))*1
	cls = rep("#6495ed", length(tmp))
	cls[which(tmp==1)] = "#ffb90f"
	rect(xleft=start(GR1), xright=start(GR1)+99, ytop = 1:length(GR1)-j, ybot = 0:(length(GR1)-1)-j, col=cls, border=cls)
	rect(xleft=start(GR2), xright=start(GR2)+99, ytop = 1:length(GR2)-j, ybot = 0:(length(GR2)-1)-j, col=cls, border=cls)
	segments(x0 = start(GR1)+99, x1=start(GR2), y0 = 1:length(GR1)-j-0.5, y1 = 1:length(GR2)-j-0.5, lwd=0.1, col= cls)
	
	segments(x0=start(track), x1=end(track), y0 = length(GR1)-length(GR1)/20, y1 = length(GR1)-length(GR1)/20, lty=2, col=coll)
	text(x=(start(track)+end(track))/2, y=length(GR1)-length(GR1)/20, labels=width(track), cex=scale, col=coll)
	abline(v=c(start(track), end(track)), lty=2, col=coll)
	text(x=min(start(bait)), y=length(GR1), labels = denote, cex=scale, col=coll)
}

visualizeFamily = function(famid, bait, track, window, row_inds, col_inds, pD, scale=4, col1 = rgb(0, 146, 146, max=255), col2 = rgb(255, 109, 182, max=255)){
	sub = pD[pD$family_id==famid,]
	bamids = as.character(sub$bamid)
	bamids = sub$bam_path
		
	cands1 = pullReads(bamids[1], bait)
	cands2 = pullReads(bamids[2], bait)
	cands3 = pullReads(bamids[3], bait)

	res = mCounts[row_inds, ,drop=F]
	res = apply(res, 2, mean)
	res = round(res, 3)

	i1 = which(names(res)==sub$subj_id[1])
	i2 = which(names(res)==sub$subj_id[2])
	i3 = which(names(res)==sub$subj_id[3])

	layout(matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, rep(10, 3)), 2, 6, byrow=T), widths=c(0.5, 0.5, rep(5, 4)), heights = c(4, 0.1))
		par(mar=rep(0,4))
		plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
			text(x=0.6, y=0.5, srt=90, famid, cex=2, pos=3)
		plot(c(0, 1), c(0,1), xlab="", ylab="", xaxt="n", yaxt="n", type="n", main="", xaxs="i", yaxs="i", bty="n")
			text(x=0.6, y=0.5, srt=90, paste0(seqnames(bait), ":", start(bait), "-", end(bait)), cex=2, pos=3)
		par(mar=c(scale, scale*3, scale, scale))
		hist(res, main="", xlim=c(-6.5, 2.5), breaks=seq(-10, 5, by=0.35), ylab="", col="grey", cex.lab=scale*.8, cex.axis=scale*.7, cex.main=scale, cex.sub=scale, xlab="")
		rect(xleft=res[i1]-0.07, xright=res[i1]+0.07, ytop=seq(100, 3000, by=200), ybot= seq(200, 3000, by=200), col=1)
		rect(xleft=res[i2]-0.07, xright=res[i2]+0.07, ytop=seq(0, 3000, by=200), ybot= seq(100, 3000, by=200), col=col1, bord=col1)
		rect(xleft=res[i3]-0.07, xright=res[i3]+0.07, ytop=seq(100, 3000, by=200), ybot= seq(200, 3000, by=200), col=col2, bord=col2)

	par(xpd=F)
	par(mar=rep(5,4))
	visualize(cands1, track, window=window, "Proband", yl="", denote="B", coll=1, bait=bait, scale=scale)
	visualize(cands2, track, window=window, "Parent 1", yl="", denote="C", coll=col1, bait=bait, scale=scale)
	visualize(cands3, track, window=window, "Parent 2", yl="", denote="D", coll=col2, bait=bait, scale=scale)	

	par(mar=rep(0,4))
	plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
	plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
	plot(c(0,1), c(0,1), type="n", axes=F, xlab="", ylab="")
	par(xpd=NA)
	plot(c(0,3), c(0,1), type="n", axes=F, xlab="", ylab="")
	rect(xleft=0.6, xright=0.7, ybot=1.2, ytop=1.4, bord="#ffb90f", col="#ffb90f")
	rect(xleft=1.7, xright=1.8, ybot=1.2, ytop=1.4, bord="#6495ed", col="#6495ed")
	text(x=1, y=1.3, labels="+ Strand Reads", cex=scale*.8)
	text(x=2.1, y=1.3, labels="- Strand Reads", cex=scale*.8)
	par(xpd=F)
}
