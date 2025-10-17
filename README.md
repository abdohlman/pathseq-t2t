# pathseq-t2t

A host subtraction, quality control, and classification pipeline for identifying low-biomass microbial signals in sequencing data.

---

## Required software
- `samtools >=1.16`
- `gatk >=4`
- `java 17`
- `picard`
- `bwa`
- `kraken2`
- `metaphlan4`
- `bowtie2`

---

## Required databases

* **PathSeq host k-mer database** \
  Must contain both `pathseq_host.bfi` (8.6G) and `pathseq_host.fa.img` (6.6G). Available from GATK Best Practices resource bundles:  
  <https://console.cloud.google.com/storage/browser/gatk-best-practices/pathseq/resources>

* **T2T-CHM13 reference genome** \
  Available via NCBI FTP (1G):  
  <https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz>

* **Kraken2 database** \
  There are several prebuilt versions, but we recommend PlusPF (~100G), available here:  
  <https://benlangmead.github.io/aws-indexes/k2>

* **MetaPhlAn4 databases** \
  Requires the CHOCOPhlAn SGB database (~3G) and bowtie2 index (~20G), available here: 
  <http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/>


---

## Installation

### Option 1: Clone and run directly
Clone this repository and add the CLI to your `PATH`:

```
git clone https://github.com/<your-org>/pathseq-t2t.git
cd pathseq-t2t
chmod +x src/pathseq-t2t
export PATH=$PWD/src:$PATH
pathseq-t2t --help
```

### Option 2: Conda (coming soon)


---

## Description

PathSeq-T2T is broken into five steps, given by five commands:

1. **`prefilter`**: selects hg38-unmapped reads (with optional decoy regions excluded).  
2. **`qcfilter`**: performs quality/complexity filtering and host k-mer subtraction using GATK PathSeqFilterSpark.  
3. **`t2tfilter`**: subtractive alignment against the complete T2T-CHM13 human genome reference.  
4. **`classify`**: microbial classification using Kraken2 and/or MetaPhlAn4 (writes raw reports to `classification_stats/`).  
5. **`summarize`**: generates filtering and classification metrics and writes normalized classification tables to `results/`.

Please note that each of these steps requires different amounts of memory, so if you are running this on an HPC environment, its recommended you run these steps separately. Memory recommendations are given below:

prefilter: 50-100Mb
qcfilter: 32-64Gb
t2tfilter: 4-8Gb
classify: 128Gb for kraken2, 32-64Gb for MetaPhlAn4

---

## Usage examples

### Step 1. Prefilter
`pathseq-t2t prefilter --input-bam sample.bam --regions-to-exclude decoys.bed --aligner bwa`

### Step 2. QC filter

```
pathseq-t2t qcfilter \ 
  --input-unaligned pst2t_output/bams/sample.prefilter.unaligned.bam \ 
  --input-excluded pst2t_output/bams/sample.prefilter.excluded.bam \ 
  --hostdir /refs/pathseq_host
```

### Step 3. T2T filter
```
pathseq-t2t t2tfilter \ 
  --input-paired pst2t_output/bams/sample.qcfilt_paired.bam \ 
  --input-unpaired pst2t_output/bams/sample.qcfilt_unpaired.bam \ 
  --reference /refs/t2t.fa
```

### Step 4. Classification
```
pathseq-t2t classify \ 
  --input-paired pst2t_output/bams/sample.t2tfilt_paired.bam \ 
  --input-unpaired pst2t_output/bams/sample.t2tfilt_unpaired.bam \ 
  --classifier both \ 
  --kraken-index /db/kraken \ 
  --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202403 \ 
  --bowtie2-index /db/bowtie2
```

### Step 5. Summarize
```
pathseq-t2t summarize \ 
  --sample-id SAMPLE123 \ 
  --filter-stats-dir $OUTDIR/filter_stats \ 
  --classification-stats-dir $OUTDIR/classification_stats \ 
  --results-dir $OUTDIR/results \ 
  --verbose
```

---

## Step 1. Prefilter

This step selects reads not aligned in the input BAM file. Typically your BAM will be aligned to hg38 prior to delivery. Please check the exact aligner was used (e.g. BWA versus DRAGEN) and which alignment parameters were used, because this will have significant affects on pre-filtering.

```
pathseq-t2t prefilter \ 
  --input-bam sample.bam \ 
  --regions-to-exclude decoys.bed \ 
  --aligner bwa
```

**Options**
* `--aligner bwa|dragen    The aligner used to generate the unput bam (required)`
* `--regions-to-exclude    <bed> (required; use None to disable)`
* `--threads               <int> (auto-detected if omitted)`
* `--sample-id             <string> (default: basename of input)`  

**Outputs**
* `<sample>.prefilter.unaligned.bam` 
* `<sample>.prefilter.excluded.bam`
* `<sample>.flagstat.tsv`

---

## Step 2. QC filter

This step performs quality/complexity filtering, masking, and screens reads for host-matching k-mers.

```
pathseq-t2t qcfilter \ 
  --input-unaligned pst2t_output/bams/sample.prefilter.unaligned.bam \ 
  --input-excluded pst2t_output/bams/sample.prefilter.excluded.bam \ 
  --hostdir /refs/pathseq_host
```

**Options**
* `--sample-id              <string> (default: basename of input)`
* `--ram-gb                 <int> (default: 16)`
* `--threads                <int> (auto-detected)`
* `--min-clipped-read-length<int> (default: 60)`
* `--dont-overwrite         Skip step if outputs already exist`
* `--keep-intermediate      Retain intermediate files`
* `--psfilterspark-args     "<extra args>" to pass to GATK PathSeqFilterSpark`
* `--picard-jar             </path/picard.jar>`

**Outputs**
* `<sample>.qcfilt_paired.bam`
* `<sample>.qcfilt_unpaired.bam`
* `<sample>.excluded.filter_metrics.txt`
* `<sample>.unaligned.filter_metrics.txt`

---

## Step 3. T2T filter

This step performs subtractive alignment to T2T-CHM13.

```
pathseq-t2t t2tfilter \ 
  --input-paired pst2t_output/bams/sample.qcfilt_paired.bam \ 
  --input-unpaired pst2t_output/bams/sample.qcfilt_unpaired.bam \ 
  --reference /refs/t2t.fa
```

**Options**
* `--sample-id          <string> (default: basename of input)`
* `--threads            <int> (auto-detected)`
* `--dont-overwrite     Skip step if outputs already exist`
* `--keep-intermediate  Retain intermediate files`
* `--picard-jar         </path/picard.jar>`

**Outputs**
* `<sample>.t2tfilt_paired.bam`
* `<sample>.t2tfilt_unpaired.bam`
* `<sample>.qcfilt_paired.t2t_aln.flagstat.tsv`
* `<sample>.qcfilt_paired.t2t_unaln.flagstat.tsv`
* `<sample>.qcfilt_unpaired.t2t_aln.flagstat.tsv`
* `<sample>.qcfilt_unpaired.t2t_unaln.flagstat.tsv`

---

## Step 4. Classification

This step classifies filtered, putatively non-human sequencing reads with Kraken2 and MetaPhlAn4. Raw classifier outputs are written to `$OUTDIR/classification_stats/`.

```
pathseq-t2t classify \ 
  --input-paired pst2t_output/bams/sample.t2tfilt_paired.bam \ 
  --input-unpaired pst2t_output/bams/sample.t2tfilt_unpaired.bam \ 
  --classifier both \ 
  --kraken-index /db/kraken \ 
  --metaphlan-index mpa_vJun23_CHOCOPhlAnSGB_202403 \ 
  --bowtie2-index /db/bowtie2
```

**Options**
* `--classifier         kraken|metaphlan|both (default: kraken)`
* `--sample-id          <string> (default: basename of input)`
* `--threads            <int> (auto-detected)`
* `--dont-overwrite     Skip step if outputs already exist`
* `--keep-intermediate  Retain intermediate files`
* `--kraken-index       <dir> or $KRAKEN_INDEX`
* `--metaphlan-index    <name> or $METAPHLAN_INDEX`
* `--bowtie2-index      <dir> or $BOWTIE2_INDEX`
* `--kraken-args        "<extra args>" for Kraken2`
* `--metaphlan-args     "<extra args>" for MetaPhlAn4`

**Kraken2 outputs (in `$OUTDIR/classification_stats/`)**
* `<sample>.paired.kraken.output.txt`
* `<sample>.paired.kraken.report.txt`
* `<sample>.unpaired.kraken.output.txt`
* `<sample>.unpaired.kraken.report.txt`

**MetaPhlAn4 outputs (in `$OUTDIR/classification_stats/`)**
* `<sample>.metaphlan.report.txt`
* `<sample>.metaphlan.bowtie2.bz2`

---

## Step 5. Summarize

This step generates a combined filtering/classification summary and normalized classification tables. Final, outputs are written to `$OUTDIR/results/`.

```
pathseq-t2t summarize   \
  --sample-id <sample>   \
  --filter-stats-dir $OUTDIR/filter_stats  \
  --classification-stats-dir $OUTDIR/classification_stats \
  --results-dir $OUTDIR/results \
  --verbose
```

**Options**
* `--sample-id                 <string> (required)`
* `--filter-stats-dir          <dir> (default: $OUTDIR/filter_stats)`
* `--classification-stats-dir  <dir> (default: $OUTDIR/classification_stats)`
* `--results-dir               <dir> (default: $OUTDIR/results)`
* `-v|--verbose`

**Outputs (to `$OUTDIR/results/`)**
* `<sample>.summary.tsv`  (combined filtering + classification totals)
* `<sample>.kraken.txt`   (normalized Kraken2 table with RPM calculations)
* `<sample>.metaphlan.txt` (normalized MetaPhlAn table with RPM calcuations)


