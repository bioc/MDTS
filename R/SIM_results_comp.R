####################################################################################
### Probes:MDTS 
rm(list=ls())
library(GenomicRanges); library(stringr); library(rtracklayer)
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/")
load("/users/jmfu/trios/good_lanes/pD_good.rda")
denovo_files = list.files("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo")
canoes_MDTS_tp = matrix(0, ncol=5, nrow=1000)
canoes_MDTS_in = matrix(0, ncol=5, nrow=1000)
canoes_MDTS_fp = GRanges()
rm_list = list()
fp_list = list()
cand_list = list()
fams = str_replace(denovo_files, "-.*", "")
rdaz = str_replace(denovo_files, ".bed", "")

percentOverlap = function(target, hits){
	ints = intersect(hits, target, ignore.strand=T)
	perc = 0
	if(length(ints)>0){perc = sum(width(ints))/width(target)}
	return(perc)
}

reduceSegs = function(seg, gr){
	gr_sub = gr[gr$segs==seg]
	return(c(as.numeric(as.character(seqnames(gr_sub)[1])), min(start(gr_sub)), max(end(gr_sub))))
}

for(i in c(1:601, 603:1000)){
	print(i); flush.console()
	load(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/CANOES/MDTS/", rdaz[i], ".rda"))	
	famid = str_extract(rdaz[i], "DS[0-9]*")
	proband = info[[1]]
	parents = do.call(c, info[2:3])
	parents = reduce(parents)

	# Removing inherited deletions
	ol_parents = ((sapply(proband, percentOverlap, parents))>0.25)*1
	# ol_parents = ((sapply(proband, percentOverlap, parents))>0)*1
	if(sum(ol_parents)>0){
		gr = proband[-which(ol_parents==1)]
	}else{
		gr = proband
	}

	denovo = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo/", rdaz[i], ".bed")); seqlevels(denovo) = str_replace(seqlevels(denovo), "chr", "")
	inherited = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/inherited/", rdaz[i], ".bed")); seqlevels(inherited) = str_replace(seqlevels(inherited), "chr", "")

	ol_denovo_perc = sapply(denovo, percentOverlap, gr)
	ol_inherited_perc = sapply(inherited, percentOverlap, gr)

	canoes_MDTS_tp[i,] = ol_denovo_perc
	canoes_MDTS_in[i,] = ol_inherited_perc

	ol_rm = queryHits(findOverlaps(gr, c(denovo, inherited)))
	fp_ind = which(! 1:length(gr) %in% ol_rm)

	cand_list[[i]] = gr
	rm_list[[i]] = ol_rm
	fp_list[[i]] = fp_ind

	tmp = gr[fp_ind]
	if(length(tmp)>0){
		tmp$id = rdaz[i]	
	}	
	canoes_MDTS_fp = c(canoes_MDTS_fp, tmp)
}
apply(canoes_MDTS_tp, 2, mean)
apply(canoes_MDTS_in, 2, mean)
canoes_MDTS_fp
save(canoes_MDTS_tp, canoes_MDTS_in, canoes_MDTS_fp, file="/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/CANOES/MDTS/performance.rda")

####################################################################################
### CANOES:PROBES
rm(list=ls())
library(GenomicRanges); library(stringr); library(rtracklayer)
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/")
load("/users/jmfu/trios/good_lanes/pD_good.rda")
denovo_files = list.files("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo")
canoes_probes_tp = matrix(0, ncol=5, nrow=1000)
canoes_probes_in = matrix(0, ncol=5, nrow=1000)
canoes_probes_fp = GRanges()
fp_list = list()
rm_list = list()
cand_list = list()
fams = str_replace(denovo_files, "-.*", "")
rdaz = str_replace(denovo_files, ".bed", "")

percentOverlap = function(target, hits){
	ints = intersect(hits, target, ignore.strand=T)
	perc = 0
	if(length(ints)>0){perc = sum(width(ints))/width(target)}
	return(perc)
}

reduceSegs = function(seg, gr){
	gr_sub = gr[gr$segs==seg]
	return(c(as.numeric(as.character(seqnames(gr_sub)[1])), min(start(gr_sub)), max(end(gr_sub))))
}

for(i in c(1:601, 603:1000)){
	print(i); flush.console()
	load(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/CANOES/PROBES/", rdaz[i], ".rda"))	
	famid = str_extract(rdaz[i], "DS[0-9]*")
	proband = info[[1]]
	parents = do.call(c, info[2:3])
	parents = reduce(parents)

	# Removing inherited deletions
	ol_parents = ((sapply(proband, percentOverlap, parents))>0.25)*1
	# ol_parents = ((sapply(proband, percentOverlap, parents))>0)*1
	if(sum(ol_parents)>0){
		gr = proband[-which(ol_parents==1)]
	}else{
		gr = proband
	}

	denovo = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo/", rdaz[i], ".bed")); seqlevels(denovo) = str_replace(seqlevels(denovo), "chr", "")
	inherited = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/inherited/", rdaz[i], ".bed")); seqlevels(inherited) = str_replace(seqlevels(inherited), "chr", "")

	ol_denovo_perc = sapply(denovo, percentOverlap, gr)
	ol_inherited_perc = sapply(inherited, percentOverlap, gr)

	canoes_probes_tp[i,] = as.numeric(ol_denovo_perc)
	canoes_probes_in[i,] = as.numeric(ol_inherited_perc)

	ol_rm = queryHits(findOverlaps(gr, c(denovo, inherited)))

	fp_ind = which(! 1:length(gr) %in% ol_rm)
	rm_list[[i]] = ol_rm
	fp_list[[i]] = fp_ind
	cand_list[[i]] = gr
	tmp = gr[fp_ind]
	if(length(tmp)>0){
		tmp$id = rdaz[i]	
	}	
	canoes_probes_fp = c(canoes_probes_fp, tmp)
}
apply(canoes_probes_tp, 2, mean)
apply(canoes_probes_in, 2, mean)
canoes_probes_fp
save(canoes_probes_tp, canoes_probes_in, canoes_probes_fp, file="~/canoes_probes.rda")
