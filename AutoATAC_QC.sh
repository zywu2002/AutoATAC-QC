#!/bin/bash

# ATACseqQC pipeline requires:
# 1. ATACseqQC
# 2. futile.logger
# 3. ChIPpeakAnno
# Please ensure that the packages mentioned above have been installed in your R library

GFFfile=/path/to/your/gff/file
BAMdir=/path/to/your/BAM/directory
rJobs=10	# Number of parallel jobs

export R_LIBS=/path/to/your/R/library

process_file() {
    local BAMfile=$1
    local GFFfile=$2
    local sample=$(basename "${BAMfile%.bam}")

    /path/to/your/Rscript /path/to/AutoATAC_PIPE.R "$BAMfile" "$GFFfile" \
     > "${sample}.out" 2> "${sample}.err"
}

export -f process_file

find $BAMdir -name *.bam | parallel --jobs $rJobs process_file {} $GFFfile