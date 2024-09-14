#!/usr/bin/env Rscript

## Passing BAM files for ATAC QC analysis. 
## An index file per BAM file must be included in the same directory

## Loading all required packages

suppressPackageStartupMessages({
  library(ATACseqQC)
  library(futile.logger)
  library(ChIPpeakAnno)
})

invisible(flog.threshold(INFO))


## Getting the BAM file path, sample ID and GFF file path
args <- commandArgs(trailingOnly=TRUE)
BAMfile <- args[1]
BAMfile.sample.ID <- gsub(".bam", "", basename(BAMfile))
GFFfile <- args[2]

flog.info("ATACseqQC analysis for sample: %s", BAMfile.sample.ID)

outPath <- BAMfile.sample.ID
if (!dir.exists(outPath)){
  dir.create(outPath)
}


## Plotting size distribution of fragments
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".fragment.size.distribution.pdf")),
    width=10, height=8)
fragSize <- ATACseqQC::fragSizeDist(bamFiles=BAMfile, 
                                    bamFiles.labels=BAMfile.sample.ID)
invisible(dev.off())

flog.info("Fragment size distribution analysis completed successfully")


## Reading in paired-end read alignment
gal <- ATACseqQC::readBamFile(BAMfile, 
                              asMates=TRUE, 
                              bigFile=TRUE)

flog.info("Paired-end read alignment (BAM) import completed successfully")


## Shifting the coordinates of 5' ends of the aligned reads in the bam file,
## +4 for reads mapping to the positve strand, -5 for reads mapping to the negative strand
gal1 <- ATACseqQC::shiftGAlignmentsList(gal)

flog.info("5' end coordinate shift completed successfully")


## Getting information of transcripts
txdb <- GenomicFeatures::makeTxDbFromGFF(GFFfile, format="gff3")
txs <- GenomicFeatures::transcripts(txdb)

flog.info("Transcript information (GFF) import completed successfully")


## Calculating Promoter/Transcript body (PT) score
pt <- ATACseqQC::PTscore(gal1, txs)
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".PTscore.pdf")),
    width=10, height=8)
plot(pt$log2meanCoverage, pt$PT_score,
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript")
invisible(dev.off())

flog.info("Promoter/Transcript body (PT) score analysis completed successfully")


## Calculating Nucleosome Free Regions (NFR) score
nfr <- ATACseqQC::NFRscore(gal1, txs)
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".NFRscore.pdf")),
    width=10, height=8)
plot(nfr$log2meanCoverage, nfr$NFR_score,
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10,0), ylim=c(-5,5))
invisible(dev.off())

flog.info("Nucleosome Free Regions (NFR) score analysis completed successfully")


## Calculating Transcription Start Site (TSS) Enrichment Score
tsse <- ATACseqQC::TSSEscore(gal1, txs)
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".TSSEscore.pdf")),
    width=10, height=8)
plot(100*(-9:10-.5), tsse$values, type="b",
     xlab="distance to TSS",
     ylab="aggregate TSS socre")
invisible(dev.off())

flog.info("Transcription Start Site (TSS) Enrichment Score analysis completed successfully")


## Splitting reads into nucleosome-free, mono-, di- and tri-nucleosome
objs <- ATACseqQC::splitGAlignmentsByCut(gal1, 
                                         txs=txs, 
                                         outPath=outPath)

flog.info("Nucleosome classifiction analysis completed successfully")


## Heatmap and coverage curve for nucleosome-free and mono-, di- and tri-nucleosome regions
bamfiles <- file.path(outPath,
                      c("NucleosomeFree.bam",
                        "mononucleosome.bam",
                        "dinucleosome.bam",
                        "trinucleosome.bam"))

# Extracting TSSs coordinates
TSS <- GenomicFeatures::promoters(txs, 
                                  upstream=0, 
                                  downstream=1)
TSS <- unique(TSS)

# Estimating the library size for normalization
librarySize <- ChIPpeakAnno::estLibSize(bamfiles)

# Calculating the signals around TSSs
NTILE <- 101
dws <- ups <- 1010

sigs <- ATACseqQC::enrichedFragments(bamfiles,
                                     TSS=TSS,
                                     librarySize=librarySize,
                                     TSS.filter=0.5,
                                     n.tile=NTILE,
                                     upstream=ups,
                                     downstream=dws)

# log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

# Plotting heatmap showing signals for nucleosome-free and mononucleosome regions around TSSs
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".TSSsignal.heatmap.pdf")),
    width=10, height=16)
ChIPpeakAnno::featureAlignedHeatmap(sigs.log2, 
                                    ChIPpeakAnno::reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt=.5, 
                                    n.tile=NTILE)
invisible(dev.off())

flog.info("TSS signals heatmap visualization completed successfully")

# Plotting curve showing signals for nucleosome-free and mononucleosome regions around TSSs
out <- ChIPpeakAnno::featureAlignedDistribution(sigs.log2,
                                                ChIPpeakAnno::reCenterPeaks(TSS, width=ups+dws),
                                                zeroAt=.5, 
                                                n.tile=NTILE,
                                                type="l",
                                                ylab="Averaged coverage")
range01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}
out <- apply(out, 2, range01)
pdf(file.path(outPath, paste0(BAMfile.sample.ID, ".TSSsignal.curve.pdf")),
    width=10, height=8)
matplot(out, 
        type="l", 
        xaxt="n", 
        xlab="Position (bp)", 
        ylab="Fraction of signal",
        main=BAMfile.sample.ID)
axis(1, at=seq(0, 100, by=10)+1, 
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")

invisible(dev.off())

flog.info("TSS signals curve visualization completed successfully")
flog.info("ATACseqQC analysis for sample (%s) completed successfully", BAMfile.sample.ID)