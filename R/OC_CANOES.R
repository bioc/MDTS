########################################################################################################################
## MDTS Bins
rm(list=ls())
library(GenomicRanges); library(stringr); library(DNAcopy); library(MDTS)
source("http://www.columbia.edu/~ys2411/canoes/CANOES.R")
	source("/dcl01/beaty/data/ingo/alternative_methods/canoes.R")
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
load('/dcl01/beaty/data/ingo/10_160_bins.rda')
load('/dcl01/beaty/data/ingo/10_160_counts.rda')

gc = bins$GC
canoes.reads = counts
gr = cbind(str_replace(as.character(seqnames(bins)), "chr", ""), start(bins), end(bins))
target = seq(1, nrow(canoes.reads))
	target = data.frame(target); names(target) = "target"
	gr = data.frame(gr); gr = apply(gr, 2, as.numeric); names(gr) = c("chromosome", "start", "end")
	gc = data.frame(gc); names(gc) = "gc"
	canoes.reads = data.frame(canoes.reads)
	sample.names = pD$subj_id
	names(canoes.reads) = sample.names
canoes.reads = cbind(target, gr, gc, canoes.reads)

	mean.counts = rep(0, length(sample.names))
	for(i in 1:length(sample.names)){
		mean.counts[i] = mean(canoes.reads[,sample.names[i]])
	}
	mean.counts = mean(mean.counts)
	for(i in 1:length(sample.names)){
		tmp = canoes.reads[,sample.names[i]]
		canoes.reads[,sample.names[i]] = round(tmp * mean.counts / mean(tmp))	
	}
	colnames(canoes.reads)[2:4] = c("chromosome", "start", "end")
	# cov_out = cor(canoes.reads[,sample.names], canoes.reads[,sample.names])
	# save(cov_out, file='/dcl01/beaty/data/ingo/alternative_methods/CANOES/covariance_probes.rda')
load(file='/dcl01/beaty/data/ingo/alternative_methods/CANOES/covariance_probes.rda')

infer = function(sample.name){
	print(which(sample.names==sample.name)); flush.console()
	if(file.exists(sample.name)==F){
		file.create(sample.name)
		xcnv = CallCNVs(sample.name, canoes.reads, cov=cov_out)
		save(xcnv, file=sample.name)
	}
}
setwd('/dcl01/beaty/data/ingo/alternative_methods/CANOES')
diag = mclapply(sample.names, infer, mc.cores = 15)

# info = lapply(xcnv.list, getInfo)
# grz = lapply(info, getGR, bins_clean)
# info = lapply(info, convertInfo, bins_clean)
# lens = unlist(lapply(info, length))
# output = do.call(c, info)
# output$ID = rep(colnames(canoes.reads)[inds_list], times=lens)
# save(output, file=outpath)
# canoes.reads[,inds_list] = old 

########################################################################################################################
## Probes
rm(list=ls())
library(GenomicRanges); library(stringr); library(DNAcopy); library(MDTS)
source("http://www.columbia.edu/~ys2411/canoes/CANOES.R")
	source("/dcl01/beaty/data/ingo/alternative_methods/canoes.R")
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
load('/dcl01/beaty/data/ingo/alternative_methods/CANOES:p/counts.rda')
load("/dcl01/beaty/data/ingo/bins_probes.rda")
	seqlevels(bins) = stringr::str_replace(seqlevels(bins), "chr", "")

gc = bins$GC
canoes.reads = counts
gr = cbind(str_replace(as.character(seqnames(bins)), "chr", ""), start(bins), end(bins))
target = seq(1, nrow(canoes.reads))
	target = data.frame(target); names(target) = "target"
	gr = data.frame(gr); gr = apply(gr, 2, as.numeric); names(gr) = c("chromosome", "start", "end")
	gc = data.frame(gc); names(gc) = "gc"
	canoes.reads = data.frame(canoes.reads)
	sample.names = pD$subj_id
	names(canoes.reads) = sample.names
canoes.reads = cbind(target, gr, gc, canoes.reads)

	mean.counts = rep(0, length(sample.names))
	for(i in 1:length(sample.names)){
		mean.counts[i] = mean(canoes.reads[,sample.names[i]])
	}
	mean.counts = mean(mean.counts)
	for(i in 1:length(sample.names)){
		tmp = canoes.reads[,sample.names[i]]
		canoes.reads[,sample.names[i]] = round(tmp * mean.counts / mean(tmp))	
	}
	colnames(canoes.reads)[2:4] = c("chromosome", "start", "end")
	# cov_out = cor(canoes.reads[,sample.names], canoes.reads[,sample.names])
# save(cov_out, file='/dcl01/beaty/data/ingo/alternative_methods/CANOES:p/covariance_probes.rda')
load('/dcl01/beaty/data/ingo/alternative_methods/CANOES:p/covariance_probes.rda')

infer = function(sample.name){
	print(which(sample.names==sample.name)); flush.console()
	if(file.exists(sample.name)==F){
		file.create(sample.name)
		xcnv = CallCNVs(sample.name, canoes.reads, cov=cov_out)
		save(xcnv, file=sample.name)
	}
}
setwd('/dcl01/beaty/data/ingo/alternative_methods/CANOES:p')
diag = mclapply(sample.names, infer, mc.cores = 10)
