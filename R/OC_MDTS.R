########################################################################################################################
### MDTS:MDTS
rm(list=ls())
library(MDTS)
# library(BSgenome.Hsapiens.UCSC.hg19)
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
# genome = BSgenome.Hsapiens.UCSC.hg19
# map_file = "~/trios/wgEncodeCrgMapabilityAlign100mer.bigWig"
# bins = calcBins(pD, 25, 100, 160, 10, genome, map_file) # 2.7 hours
# 	save(bins, file='/dcl01/beaty/data/ingo/10_160_bins.rda')
load('/dcl01/beaty/data/ingo/10_160_bins.rda')
# counts = calcCounts(pD, bins, 100, mc.cores=15) # 1.5 hours
# 	save(counts, file='/dcl01/beaty/data/ingo/10_160_counts.rda')
load('/dcl01/beaty/data/ingo/10_160_counts.rda')
# mCounts = normalizeCounts(counts, bins, mc.cores=15) 
# 	save(mCounts, file='/dcl01/beaty/data/ingo/10_160_mCounts.rda')
# md = calcMD(mCounts, bins, pD)
# 	save(md, file='/dcl01/beaty/data/ingo/10_160_md.rda')
load(file='/dcl01/beaty/data/ingo/10_160_mCounts.rda')
load(file='/dcl01/beaty/data/ingo/10_160_md.rda')
	vars = apply(md, 2, var)
	# lag10 = apply(md, 2, function(x) acf(x, 10)$acf[10])
	bad_families = names(vars)[which(vars>0.05)]
set.seed(3350)
cbs = segmentMD(md, bins, mc.cores = 15)
denovo = denovoDeletions(cbs, mCounts, bins)
	denovo = denovo[-which(denovo$famid %in% bad_families)]
save(denovo, file='/dcl01/beaty/data/ingo/10_160_denovo.rda')

# pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
# load('/dcl01/beaty/data/ingo/10_160_mCounts.rda')
# load('/dcl01/beaty/data/ingo/10_160_bins.rda')
# load('/dcl01/beaty/data/ingo/10_160_denovo.rda')
# cnps = findCNP(mCounts, bins, pD)
# visualizeDeletion(denovo[1], bins, pD, mCounts, save=TRUE)

########################################################################################################################
### MDTS:p
rm(list=ls())
library(MDTS)
pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
load("/dcl01/beaty/data/ingo/bins_probes.rda")
	seqlevels(bins) = stringr::str_replace(seqlevels(bins), "chr", "")
counts = calcCounts(pD, bins, rl=100, mc.cores=1) # 1.5 hours
	save(counts, file='/dcl01/beaty/data/ingo/alternative_methods/CANOES:p/counts.rda')
mCounts = normalizeCounts(counts, bins, mc.cores=20) 
md = calcMD(mCounts, bins, pD)
	vars = apply(md, 2, var)
	bad_families = names(vars)[which(vars>0.05)]
set.seed(3350)
cbs = segmentMD(md, bins, mc.cores = 20)
	save(cbs, file='/dcl01/beaty/data/ingo/MDTS_probes_cbs.rda')
denovo = denovoDeletions(cbs, mCounts, bins)
	denovo = denovo[-which(denovo$famid %in% bad_families)]
save(denovo, file='/dcl01/beaty/data/ingo/MDTS_probes_denovo.rda')

# pD = pData('/dcl01/beaty/data/ingo/pD_MDTS.tab')
# load('/dcl01/beaty/data/ingo/10_160_mCounts.rda')
# load('/dcl01/beaty/data/ingo/10_160_bins.rda')
# load('/dcl01/beaty/data/ingo/10_160_denovo.rda')
# cnps = findCNP(mCounts, bins, pD)
# visualizeDeletion(denovo[1], bins, pD, mCounts, save=TRUE)