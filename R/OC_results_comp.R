### TrioCNV
library(GenomicRanges)
setwd('/dcl01/beaty/data/ingo/alternative_methods/TrioCNV/OurBins')

output = GRanges()
files=list.files()
for(f in files){
	tab = read.table(f, header=T)
	gr = GRanges(seqnames=tab[,1], IRanges(tab[,2], tab[,3]), state=paste0(tab[,4], '-', tab[,5], '-', tab[,6]))
	gr$fam = stringr::str_replace(f, ".vcf", "")
	output = c(output, gr)
}
load(file='/dcl01/beaty/data/ingo/10_160_md.rda')
	vars = apply(md, 2, var)
	lag10 = apply(md, 2, function(x) acf(x, 10)$acf[10])
bad_families = names(vars)[which(vars>0.05)]


DS10826 = output[output$fam=="DS10826"]
DS11025 = output[output$fam=="DS11025"]
targets = GRanges(c(1,8,8), IRanges(c(209945655, 129614522, 130113612), c(209947210, 129616078, 130132753)))
findOverlaps(targets, DS10826)

### TrioCNV
java -jar ~/downloads/TrioCNV-0.1.2/TrioCNV-0.1.2.jar call -I preprocess -P phenotype.ped -O DS10826_g.vcf --min_distance 1000 -M ~/downloads/wgEncodeCrgMapabilityAlign100mer.bigWig

tab = read.table('/dcl01/beaty/data/ingo/alternative_methods/DS10826_g.vcf', header=TRUE)
# 296,772 hits || 20,721 DDD, 275,524 DDDup
gr = GRanges(seqnames=tab[,1], IRanges(tab[,2], tab[,3]), state=paste0(tab[,4], '-', tab[,5], '-', tab[,6]))
gr[gr$state=='REF-REF-DEL']
findOverlaps(targets, gr)
