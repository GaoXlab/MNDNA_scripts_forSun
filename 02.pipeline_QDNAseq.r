library(QDNAseq)
library(QDNAseq.hg38)
library(QDNAseq.mm10)
library(matrixStats)
library(Biobase)
library(stringr)
library(ggplot2)


# General
args = commandArgs(trailingOnly=TRUE)

data_path = args[1]
binsize = args[3]
species = args[4]
bam_in = args[2]

#definded binsize
bin.size <- as.numeric(binsize)
print(paste("bam_in:", bam_in))
## 1. generate annotation file either by preloading calculated files or generating new one.
if(species=='hg38' || species=='mm10'){
    binFilename=paste0(species, '.bins.log')
    bins <- getBinAnnotations(binSize=bin.size, genome = species)
    if ( !file.exists(binFilename) ) {
        write.table(pData(bins), binFilename, sep='\t')
    }
}
## 2. caculate readcounts.
readCounts <- binReadCounts(bins, bamfile=paste0(data_path,bam_in))
exportBins(readCounts, file=paste0(bam_in,'.',bin.size,"kb_readCounts.txt"),logTransform = FALSE) ##export raw readCounts
## 3. Apply filters and Estimate the correction for GC content and mappability, and for the relationship between the observed standard deviation in the data and its read depth.
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)
exportBins(readCountsFiltered, file=paste0(bam_in,'.',bin.size,"kb_copyNumbersCorrect.txt"),logTransform = FALSE) ##export raw readCounts
## 4. Apply the correction for GC content and mappability.
copyNumbers <- correctBins(readCountsFiltered)
exportBins(copyNumbers, file=paste0(bam_in,'.',bin.size,"kb_copyNumbers.txt"),logTransform = FALSE) ##export raw readCounts
## 5. Apply the normalization (median value normalization).
copyNumbersNormalized <- normalizeBins(copyNumbers)
exportBins(copyNumbersNormalized, file=paste0(bam_in,'.',bin.size,"kb_copyNumbersNormalized.txt"),logTransform = FALSE) ##export raw readCounts
## 6. Apply the smooth outliers.
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)
exportBins(copyNumbersSmooth, file=paste0(bam_in,'.',bin.size,"kb_copyNumbersSmooth.txt"),logTransform = FALSE) ##export raw readCounts

