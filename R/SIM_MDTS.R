#####################################################################################################################
### Calling - MDTS + MDTS
rm(list=ls())
library(GenomicRanges); library(stringr); library(DNAcopy); library(rtracklayer); library(Repitools)
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
load('/dcl01/beaty/data/ingo/10_160_bins.rda')
load('/dcl01/beaty/data/ingo/10_160_counts.rda')
load(file='/dcl01/beaty/data/ingo/10_160_mCounts.rda')
load(file='/dcl01/beaty/data/ingo/10_160_md.rda')

counts_clean = counts
log_counts = log(counts_clean+1, 2)
col_med = apply(log_counts, 2, median)
test = t(t(log_counts)-col_med)
row_med = apply(test, 1, median)
	
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/c/coverage_US")
samples = list.files(); samples = samples[str_detect(samples, ".count$")]
paths = str_replace(samples, ".count", "")
bams = str_replace(paths, ".count", "")
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/")

readTab = function(path){
	print(which(bams==path)); flush.console()
	return(read.table(paste0(path, ".count"))[,1])
}
getGR = function(i){
	new = do.call(cbind, lapply(paste0(c("c/coverage_US/", "p1/coverage_US/", "p2/coverage_US/"), bams[i]), readTab))
	logNew = log(new+1, 2)
	fam = logNew
	famid = str_replace(bams[i], "-[0-9]*", "")
	col_inds = match(paste0(famid, c("_1", "_2", "_3")), colnames(counts))
	fam_normed = t(t(fam) - col_med[col_inds])
	fam_normed = fam_normed - row_med

	res_gc = fam_normed
	for(j in 1:3){
		res_gc[,j] = residuals(loess(fam_normed[,j]~bins$GC))
	}
	res_gc_map = res_gc
	loess_ind = which(bins$mappability<1)
	for(j in 1:3){
		res_gc_map[loess_ind,j] = residuals(loess(fam_normed[loess_ind,j]~bins$mappability[loess_ind]))
	}
	md1 = res_gc_map[,1] - res_gc_map[,2]
	md2 = res_gc_map[,1] - res_gc_map[,3]
	md = md1; md[abs(md1)>abs(md2)] = md2[abs(md1)>abs(md2)]
	cna = CNA(genomdat=md, chrom=as.vector(seqnames(bins)), maploc=(start(bins)+end(bins))/2, data.type="logratio", sampleid=bams[i], presorted=T)
	cbs = segment(cna, alpha=0.001, undo.splits="sdundo", undo.SD=4)
	gr = GRanges(seqnames=cbs$output$chrom, IRanges(start(bins[cbs$segRows[,1]]), end(bins[cbs$segRows[,2]])) , m=cbs$output$seg.mean, famid=bams[i])
	return(gr)
}
getDels = function(i){
	print(i); flush.console()
	out = paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/MDTS/MDTS/output.", i, ".rda")
	if(file.exists(out)==FALSE){
		file.create(out)
		output = getGR(i)
		save(output, file=out)
		rm(output)
	}
}

info = mclapply(1:1000, getDels, mc.cores=10)


#####################################################################################################################
### Calling - MDTS + probes
rm(list=ls())
library(GenomicRanges); library(stringr); library(DNAcopy); library(MDTS)
source("http://www.columbia.edu/~ys2411/canoes/CANOES.R")
	source("/dcl01/beaty/data/ingo/alternative_methods/canoes.R")
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
load('/dcl01/beaty/data/ingo/alternative_methods/CANOES:p/counts.rda')
load("/dcl01/beaty/data/ingo/bins_probes.rda")
	seqlevels(bins) = stringr::str_replace(seqlevels(bins), "chr", "")

setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/c/coverage_CANOES")
samples = list.files(); samples = samples[str_detect(samples, ".count$")]
paths = str_replace(samples, ".count", "")
bams = str_replace(paths, ".count", "")
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/")

counts_clean = counts
log_counts = log(counts_clean+1, 2)
col_med = apply(log_counts, 2, median)
test = t(t(log_counts)-col_med)
row_med = apply(test, 1, median)
	
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/c/coverage_CANOES")
samples = list.files(); samples = samples[str_detect(samples, ".count$")]
paths = str_replace(samples, ".count", "")
bams = str_replace(paths, ".count", "")
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/")

readTab = function(path){
	print(which(bams==path)); flush.console()
	return(read.table(paste0(path, ".count"))[,1])
}
getGR = function(i){
	new = do.call(cbind, lapply(paste0(c("c/coverage_CANOES/", "p1/coverage_CANOES/", "p2/coverage_CANOES/"), bams[i]), readTab))
	logNew = log(new+1, 2)
	fam = logNew
	famid = str_replace(bams[i], "-[0-9]*", "")
	col_inds = match(paste0(famid, c("_1", "_2", "_3")), colnames(counts))
	fam_normed = t(t(fam) - col_med[col_inds])
	fam_normed = fam_normed - row_med

	res_gc = fam_normed
	for(j in 1:3){
		res_gc[,j] = residuals(loess(fam_normed[,j]~bins$GC))
	}
	res_gc_map = res_gc
	loess_ind = which(bins$mappability<1)
	for(j in 1:3){
		res_gc_map[loess_ind,j] = residuals(loess(fam_normed[loess_ind,j]~bins$mappability[loess_ind]))
	}
	md1 = res_gc_map[,1] - res_gc_map[,2]
	md2 = res_gc_map[,1] - res_gc_map[,3]
	md = md1; md[abs(md1)>abs(md2)] = md2[abs(md1)>abs(md2)]
	cna = CNA(genomdat=md, chrom=as.vector(seqnames(bins)), maploc=(start(bins)+end(bins))/2, data.type="logratio", sampleid=bams[i], presorted=T)
	cbs = segment(cna, alpha=0.001, undo.splits="sdundo", undo.SD=4)
	gr = GRanges(seqnames=cbs$output$chrom, IRanges(start(bins[cbs$segRows[,1]]), end(bins[cbs$segRows[,2]])) , m=cbs$output$seg.mean, famid=bams[i])
	return(gr)
}
getDels = function(i){
	print(i); flush.console()
	out = paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/MDTS/probes/output.", i, ".rda")
	if(file.exists(out)==FALSE){
		file.create(out)
		output = getGR(i)
		save(output, file=out)
		rm(output)
	}
}

info = mclapply(1:1000, getDels, mc.cores=10)