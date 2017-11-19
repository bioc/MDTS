####################################################################################
### MDTS
####################################################################################
### MDTS:MDTS
rm(list=ls())
library(GenomicRanges); library(stringr)
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/calls/MDTS/MDTS")
files = list.files(pattern = "output")
output_all = GRanges()
for(f in files){
	load(f)
	output_all = c(output_all, output)
}
fams = unique(output_all$famid)
dels = output_all[output_all$m<= -0.7 & output_all$m >= -1.3]
seqlevels(dels) = str_replace(seqlevels(dels), "chr", "")

### Filtering by 50%-overlap-bad-bins definition
load("/dcl01/beaty/data/ingo/bins_empirical/bin_info.rda")
load("/dcl01/beaty/data/ingo/bins_empirical/pD_good_res_gc_map.rda")
win = 0.5
raw = abs(res_gc_map)<win
bins = bins_clean
bins_filter = bins[which(apply(raw, 1, mean)<0.95)]
	seqlevels(bins_filter) = str_replace(seqlevels(bins_filter), "chr", "")
	seqlevels(bins) = str_replace(seqlevels(bins), "chr", "")
ol_filter = findOverlaps(dels, bins_filter)
ol_bins = findOverlaps(dels, bins)
count_filter = by(rep(1, length(ol_filter)), queryHits(ol_filter), sum)
count_bins = by(rep(1, length(ol_bins)), queryHits(ol_bins), sum)
filtering_info = data.frame(id = names(count_bins), count_base = as.numeric(count_bins))
filtering_info2 = data.frame(id = names(count_filter), count_filter = as.numeric(count_filter))
filtering = merge(filtering_info, filtering_info2, by="id", sort=F)
ratios = filtering[,3]/filtering[,2]
cut = which(ratios>=0.5)
drop_ids = as.numeric(as.character(filtering$id[cut]))
dels = dels[-drop_ids]

percentOverlap = function(target, hits){
	ints = intersect(hits, target, ignore.strand=T)
	perc = 0
	if(length(ints)>0){perc = sum(width(ints))/width(target)}
	return(perc)
}

mdts_mdts_tp = matrix(0, ncol=5, nrow=1000)
mdts_mdts_in = matrix(0, ncol=5, nrow=1000)
mdts_mdts_fp = GRanges()
rm_list = list()
fp_list = list()
paths = fams
for(i in 1:length(paths)){
	print(i); flush.console()
	gr = dels[dels$famid==paths[i]]
	if(length(gr)>0){
		seqlevels(gr) = str_replace(seqlevels(gr), "chr", "")
		denovo = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo/", paths[i], ".bed")); seqlevels(denovo) = str_replace(seqlevels(denovo), "chr", "")
		inherited = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/inherited/", paths[i], ".bed")); seqlevels(inherited) = str_replace(seqlevels(inherited), "chr", "")

		ol_denovo_perc = sapply(denovo, percentOverlap, gr)
		ol_inherited_perc = sapply(inherited, percentOverlap, gr)

		mdts_mdts_tp[i,] = as.numeric(ol_denovo_perc)
		mdts_mdts_in[i,] = as.numeric(ol_inherited_perc)

		ol_rm = queryHits(findOverlaps(gr, c(denovo, inherited)))
		fp_ind = which(! 1:length(gr) %in% ol_rm)
		rm_list[[i]] = ol_rm
		fp_list[[i]] = fp_ind

		mdts_mdts_fp = c(mdts_mdts_fp, gr[fp_ind])
	}
}
apply(mdts_mdts_tp, 2, mean)
apply(mdts_mdts_in, 2, mean)
mdts_mdts_fp

save(mdts_mdts_tp, mdts_mdts_in, mdts_mdts_fp, file="~/mdts_mdts.rda")

#####################################################################################################################
### MDTS:PROBES
rm(list=ls())
library(GenomicRanges); library(stringr)
setwd("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/calls/MDTS/probes")
files = list.files(pattern = "output")
output_all = GRanges()
for(f in files){
	load(f)
	output_all = c(output_all, output)
}
dels = output_all[output_all$m<= -0.7 & output_all$m >= -1.3]
seqlevels(dels) = str_replace(seqlevels(dels), "chr", "")
fams = unique(dels$famid)

### Filtering by 50%-overlap-bad-bins definition
load("/dcl01/beaty/data/ingo/bins_probes.rda")
load(file="/dcl01/beaty/data/ingo/alternative_methods/CANOES.mCounts.rda")
win = 0.5
raw = abs(mCounts)<win
bins_filter = bins[which(apply(raw, 1, mean)<0.95)]
	seqlevels(bins_filter) = str_replace(seqlevels(bins_filter), "chr", "")
	seqlevels(bins) = str_replace(seqlevels(bins), "chr", "")
ol_filter = findOverlaps(dels, bins_filter)
ol_bins = findOverlaps(dels, bins)
count_filter = by(rep(1, length(ol_filter)), queryHits(ol_filter), sum)
count_bins = by(rep(1, length(ol_bins)), queryHits(ol_bins), sum)
filtering_info = data.frame(id = names(count_bins), count_base = as.numeric(count_bins))
filtering_info2 = data.frame(id = names(count_filter), count_filter = as.numeric(count_filter))
filtering = merge(filtering_info, filtering_info2, by="id", sort=F)
ratios = filtering[,3]/filtering[,2]
cut = which(ratios>=0.5)
drop_ids = as.numeric(as.character(filtering$id[cut]))
dels = dels[-drop_ids]

					dels = dels[width(dels)>=300]

percentOverlap = function(target, hits){
	ints = intersect(hits, target, ignore.strand=T)
	perc = 0
	if(length(ints)>0){perc = sum(width(ints))/width(target)}
	return(perc)
}

mdts_probes_tp = matrix(0, ncol=5, nrow=1000)
mdts_probes_in = matrix(0, ncol=5, nrow=1000)
mdts_probes_fp = GRanges()
rm_list = list()
fp_list = list()
paths = fams
for(i in 1:length(paths)){
	print(i); flush.console()
	gr = dels[dels$famid==paths[i]]
	if(length(gr)>0){
		seqlevels(gr) = str_replace(seqlevels(gr), "chr", "")
		denovo = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/denovo/", paths[i], ".bed")); seqlevels(denovo) = str_replace(seqlevels(denovo), "chr", "")
		inherited = rtracklayer::import.bed(paste0("/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/inherited/", paths[i], ".bed")); seqlevels(inherited) = str_replace(seqlevels(inherited), "chr", "")

		ol_denovo_perc = sapply(denovo, percentOverlap, gr)
		ol_inherited_perc = sapply(inherited, percentOverlap, gr)

		mdts_probes_tp[i,] = as.numeric(ol_denovo_perc)
		mdts_probes_in[i,] = as.numeric(ol_inherited_perc)

		ol_rm = queryHits(findOverlaps(gr, c(denovo, inherited)))
		fp_ind = which(! 1:length(gr) %in% ol_rm)
		rm_list[[i]] = ol_rm
		fp_list[[i]] = fp_ind

		mdts_probes_fp = c(mdts_probes_fp, gr[fp_ind])
	}
}
apply(mdts_probes_tp, 2, mean)
apply(mdts_probes_in, 2, mean)
mdts_probes_fp

save(mdts_probes_tp, mdts_probes_in, mdts_probes_fp, file="~/mdts_probes_300.rda")

####################################################################################
### CANOES
####################################################################################
### CANOES:MDTS 
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
save(canoes_probes_tp, canoes_probes_in, canoes_probes_fp, file="/dcl01/beaty/data/ingo/alternative_methods/simulations/1000/CANOES/PROBES/performance.rda")
