## ---- echo=F, message=F, warning=F---------------------------------------
library(MDTS); library(BSgenome.Hsapiens.UCSC.hg19)
setwd(system.file("extdata", package="MDTS"))
load('pD.RData')
pD

## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github('jmf47/MDTSdata')

## ---- echo=FALSE, message=FALSE------------------------------------------
library(MDTS)

## ---- eval=F-------------------------------------------------------------
#  library(MDTS); library(BSgenome.Hsapiens.UCSC.hg19)
#  # Using the raw data from MDTSData
#  setwd(system.file("data", package="MDTSData"))
#  
#  # Importing the pedigree file that includes information on where to locate the raw bam files
#  pD = pData("pData.ped")
#  
#  # Information on the GC content and mappability to estimate GC and mappability for the MDTS bins
#  genome = BSgenome.Hsapiens.UCSC.hg19; map_file = "chr22_100mer_map.bw"
#  
#  # This command now subsets 5 samples to determine MDTS bins
#  # pD is the phenotype matrix
#  # n is the number of samples to examine to calculate the bins
#  # rl is the sequencing read length
#  # min is the minimum read depth before a location is to be included in a proto region
#  # med is the minimum number of the median number of reads across the n samples in a bin.
#  bins = calcBins(pD, n=5, rl=100, med=100, min=5, genome, map_file)

## ---- eval=F-------------------------------------------------------------
#  # pD is the phenotype matrix
#  # bins is the previously calculated MDTS bins
#  # rl is the sequencing read length
#  counts = calcCounts(pD, bins, rl=100)

## ---- message=FALSE------------------------------------------------------
setwd(system.file("extdata", package="MDTS"))
load('bins.RData')
load('counts.RData')
load('pD.RData')

## ------------------------------------------------------------------------
bins

## ------------------------------------------------------------------------
head(counts)

## ------------------------------------------------------------------------
# counts is the raw read depth of [MDTS bins x samples]
# bins is the previously calculated MDTS bins
mCounts = normalizeCounts(counts, bins)

## ------------------------------------------------------------------------
# mCounts is the normalized read depth of [MDTS bins x samples]
# bins is the previously calculated MDTS bins
# pD is the phenotype matrix
md = calcMD(mCounts, bins, pD)

## ---- warning=FALSE, message=FALSE---------------------------------------
# md is the Minimum Distance of [MDTS bins x trio]
# bins is the previously calculated MDTS bins
# mCounts is the normalized read depth of [MDTS bins x samples]
cbs = segmentMD(md, bins)
denovo = denovoDeletions(cbs, mCounts, bins)

## ------------------------------------------------------------------------
denovo

