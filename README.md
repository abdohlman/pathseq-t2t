# pathseq-t2t

A host subtraction and filtering pipeline for identifying low-biomass microbial signals in sequencing data.

**Required software**
- `samtools >=1.16`
- `bwa`
- `gatk >=4.0`
- `picard`
- `kraken2`
- `metaphlan4`
- `bowtie2`

**Required databases**
* PathSeq host kmer database (`pathseq_host.bfi` and `pathseq_host.fa.img`)
* T2T-CHM13 reference 
* Kraken2 database
* MetaPhlAn database

## Description

PathSeq-T2T is broken into four steps, given by four commands: (1) `prefilter`, which selects of hg38-unmapped reads; (2) `qcfilter`, which performs quality/complexity filtering and additional host k-mer filtering; (3) `t2tfilter` which performs subtractive alignment against a complete human genome reference; (4) `classify`, which performs microbial classification on filtered sequences using Kraken2 and/or MetaPhlAn4.

## Usage examples

* `pathseq-t2t prefilter --input-bam sample.bam --decoy decoys.bed`

* `pathseq-t2t qcfilter --input-unaligned pst2t_output/bams/sample.prefilter.unaligned.bam \
                     --input-decoy pst2t_output/bams/sample.prefilter.decoy.bam \
                     --hostdir /refs/pathseq_host`
                     
* `pathseq-t2t t2tfilter --input-paired pst2t_output/bams/sample.qcfilt_paired.bam \
                      --input-single pst2t_output/bams/sample.qcfilt_single.bam \
                      --reference /refs/t2t.fa`
                      
* `pathseq-t2t classify --input-paired pst2t_output/bams/sample.t2t.paired.bam \
                     --input-single pst2t_output/bams/sample.t2t.single.bam \
                     --classifier both \
                     --kraken-db /db/kraken \
                     --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202403 \
                     --bowtie2db /db/bowtie2`

## Step 1. Prefilter

Performs pre-filtering on an hg38-aligned bam file.

`pathseq-t2t qcfilter \
  --input-unaligned pst2t_output/bams/sample.prefilter.unaligned.bam \
  --input-decoy pst2t_output/bams/sample.prefilter.decoy.bam \
  --hostdir /refs/pathseq_host`

**Options:**
* `--aligner bwa|dragen (default: dragen)` Specify the aligner used to generate the BAM
* `--threads <int> (auto-detected if omitted)` Number of threads to use
* `--sample-id <string> (default: basename of input)` Sample ID used to track filtering

**Outputs:**
* `<sample>.prefilter.unaligned.bam`
* `<sample>.prefilter.decoy.bam`
* `<sample>.input.flagstat.txt`

## Step 2. QC filter

Performs quality and complexity filtering, additional host k-mer filtering.

`pathseq-t2t qcfilter \
  --input-unaligned pst2t_output/bams/sample.prefilter.unaligned.bam \
  --input-decoy pst2t_output/bams/sample.prefilter.decoy.bam \
  --hostdir /refs/pathseq_host`

**Options:**
* `--sample-id <string> (default: basename of input)` Sample ID used to track filtering
* `--ram-gb <int> (default: 16)` Amount of memory to use for filtering. For large BAMs, we recommend increasing this.
* `--threads <int> (auto-detected)` Number of threads to use
* `--min-clipped-read-length <int> (default: 60)` Masked length at which a read is excluded
* `--dont-overwrite` Skip step if outputs already exist
* `--keep-intermediate` Retain intermediate files
* `--psfilterspark-args "<extra args to pass to GATK PathSeqFilterSpark>"` Other arguments to pass to GATK

**Outputs:**
* `<sample>.qcfilt_paired.bam`
* `<sample>.qcfilt_single.bam`
* `<sample>.decoy.pathseq_filter_metrics.txt` 
* `<sample>.unaligned.pathseq_filter_metrics.txt` 

## Step 3. T2T filter

Performs filtering using a compute human genome reference, T2T-CHM13

`pathseq-t2t t2tfilter \
  --input-paired pst2t_output/bams/sample.qcfilt_paired.bam \
  --input-single pst2t_output/bams/sample.qcfilt_single.bam \
  --reference /refs/t2t.fa`

**Options:**
* `--sample-id <string> (default: basename of input)` Sample ID used to track filtering
* `--threads <int> (default: auto)` Number of threads to use
* `--dont-overwrite` Skip step if outputs already exist
* `--keep-intermediate` Retain intermediate files

**Outputs:**
* `<sample>.t2t.paired.bam`
* `<sample>.t2t.single.bam`
* `<sample>.t2t.unaligned.paired.flagstat.txt`
* `<sample>.t2t.unaligned.single.flagstat.txt`


## Step 4. Classification

Performs microbial classification on filtered sequences

`pathseq-t2t classify \
  --input-paired pst2t_output/bams/sample.t2t.paired.bam \
  --input-single pst2t_output/bams/sample.t2t.single.bam \
  --classifier both \
  --kraken-db /db/kraken \
  --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202403 \
  --bowtie2db /db/bowtie2`

**Options:**
* `--classifier kraken|metaphlan|both (default: kraken)` Classifier to use
* `--sample-id <string> (default: basename of input)` Sample ID used to track filtering
* `--threads <int> (default: auto)` Number of threads to use
* `--dont-overwrite` Skip step if outputs already exist
* `--keep-intermediate` Retain intermediate files
* `--kraken-args` Other arguments to pass to Kraken2
* `--metaphlan-args` Other arguments to pass to MetaPhlAn4


**Kraken2 outputs:**
* `<sample>.kraken.paired.txt`
* `<sample>.kraken.paired.report`
* `<sample>.kraken.single.txt`
* `<sample>.kraken.single.report`

**MetaPhlAn4 outputs:**
* `<sample>.metaphlan.txt`
* `<sample>.metaphlan.bowtie2.bz2`

