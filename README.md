# AutoATAC-QC
**AutoATAC-QC** is a streamlined, one-step pipeline designed to efficiently process multiple datasets in parallel, leveraging the power of the `ATACseqQC` R package.

## Requirements
Before running the pipeline, ensure the following R packages are installed:
* `ATACseqQC`
* `futile.logger`
* `ChIPpeakAnno`

## Usage
First, configure the parameters in the `AutoATAC_QC.sh` script:
* `GFFfile`: Path to the GFF annotation file.
* `BAMdir` Directory containing the input BAM file(s).
* `rJobs`: Number of parallel jobs to run with the `AutoATAC_PIPE.R` script.

Then, run the pipeline to perform ATAC-seq quality control across multiple jobs in a single command:

```{shell}
chmod +x AutoATAC_QC.sh
bash AutoATAC_QC.sh
```
