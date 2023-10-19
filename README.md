# Brotman Baty Institute Viral MIP Assembly Pipeline

## Introduction
Assembly pipeline used by the Brotman Baty Institute to build consensus genomes for sequenced samples.

This pipeline was originally written for assembly of [Seattle Flu Study](https://seattleflu.org/) flu genomes, but has been modified to support assembly of SARS-CoV-2, RSV, and flu genomes sequenced with MIPs technology (CRIseq).

Adapted from Louise Moncla's [Illumina pipline](https://github.com/lmoncla/illumina_pipeline) for influenza snp calling.

## Installation
Install dependencies via Conda by running:
```
conda env create -f envs/seattle-flu-environment.yaml
```
Once installed, the environment needs to be activated by running:
```
conda activate seattle-flu
```

## Set Up
### FASTQ Files
The pipeline expects 2 or more FASTQ files for each individual sample, with Read 1 (R1) and Read 2 (R2) from 1 or more lanes.
Merging across lanes will be performed automatically.

FASTQs should already be trimmed of their 5’ MIP targeting arms, which can be performed during demultiplexing with bcl2fastq.

Expected filename format: `318375_S12_L001_R1_001.fastq.gz` where `318375` is the NWGC sample ID.


## Configuration
Currently, `config/config.criseq.json` has the most updated version of the configuration needed to run assembly.

Recommended file structure:  
Create a work directory, which contains a subdirectory with input fastq files.
The output directory can be specified in the config file and will be created as another subdirectory in the work directory.

Edit the following parameters in the config file:

* `work_dir`: path to the work directory, which will contain both input and output subdirectories
* `input_subdir`: relative path from the work directory to the subdirectory containing input fastqs,
* `assembly_subdir`: relative path from the work directory to the output subdirectory (will be created by Snakemake if it does not already exist)
* `ignored_samples`: an object containing keys that specify samples to be ignored (optional)
* `sample_reference_pairs`: an object containing specific sample/reference pairs, with the keys refering to the samples and an array of references as the values. (optional; if left empty, the pipeline will attempt to create an assembly for every reference in `reference_viruses` for each sample)
* `reference_virusus`: an object containing keys that specify references to be used
* `params`: parameters for the various tools used in the pipeline

## Usage
The Snakemake pipeline generates one consensus FASTA for each assembly.
To run the Snakemake pipeline:
```
snakemake --configfile config/config.criseq.json -k --cores <number_of_cores> --rerun-incomplete --snakefile Snakefile-cri
```

To combine all assemblies into one FASTA per pathogen and gather whole genome sequence metrics with Picard, run these 3 scripts manually after Snakemake.
The first 2 steps should be run in the seattle-flu conda environment, and the last step should be run in base conda, with Picard installed locally.
```
# conda activate seattle-flu (if not already activated)
bash scripts/fasta_to_single_line.sh
bash scripts/sort_index_bams.sh
conda deactivate
bash scripts/picard.sh
```

## Snakemake pipeline basic steps

### 1. Index reference genomes
Using [bowtie2-build](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) to build a Bowtie index for each reference genome. These will later be used in the mapping step.

### 2. Merge lanes
Concatenates FASTQ files for each sample into 2 files (R1 and R2).

### 3. Trim fastqs
Trim raw FASTQ files with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), which cuts out adapter/illumina-specific sequences, trims the reads based on quality scores, and removes short reads. Trim MIP 3’ overhangs with NGmerge.

### 4. Post trim fastqc
Generates summary statistics using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Allows for some quality control checks on raw sequence data coming from high throughput sequencing pipelines.

### 5. Map
Map trimmed reads to reference genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Outputs a BAM file that represents aligned sequences and a log file that contains the alignment summary.

### 6. Sort
Sort the BAM file into "genome order" using [samtools sort](http://www.htslib.org/doc/samtools.html).

### 7. Bamstats
Generates statistics for coverage using [BAMStats](http://bamstats.sourceforge.net/).

### 8. Mapped reads checkpoint
[Checkpoints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) allow for data-dependent conditional execution, forcing Snakemake to re-evaluate the DAG at this point.

Determines if all segments within the reference genome has coverage summary generated from BAMStats.
Then calculates the minimum number of reads required to meet the minimum coverage depth set for `varscan mpileup2snp` and checks the number of mapped reads from the bowtie2 alignment summary. The minimum number of reads required is calculated as:
```
min_reads = (reference_genome_length * min_cov) / (raw_read_length)
```
If not all segments are represented or if `mapped_reads` is less than `min_reads`, then stop the process for this sample/reference pair.

### 9. Not mapped
This rule is necessary for the checkpoint above to work, because it must send the data down one of two paths. This is the "dead end" where a consensus genome is __not__ generated.

Also a placeholder for how to handle sample/reference pairs that do not meet the `min_align_rate` threshold. Currently it just outputs the overall align rate of the sample/reference pair so that we can check that it did not meet the minimum.

### 10. Pileup
Generates [Pileup](https://en.wikipedia.org/wiki/Pileup_format) for BAM file using [samtools mpileup](http://www.htslib.org/doc/samtools.html).

  _Important Flags_:
  * `-a` to print out all positions, including zero depth (this is necessary for generating the low coverage Bedfile later)
  * `-A` to not discard anomalous read pairs
  * `-Q` to set the minimum base quality to consider a read (set to match the minimum in `varscan mpileup2snp` so that their coverage depths match)

### 11. Create bed file
Creates a BED file for positions that need to be masked in the consensus genome. Positions need to be masked if they are below the minimum coverage depth.

### 12. Call Variants
Calls variants (SNPs/indels) from the Pileup based on parameters set in the config file using [varscan mpileup2cns](http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2cns)

### 13. Zip VCF
Compress VCF using [bgzip](http://www.htslib.org/doc/bgzip.html), which allows indexes to be built against the file and allows the file be used without decompressing it.

### 14. Index BCF
Creates index for compressed VCF using [bcftools index](https://samtools.github.io/bcftools/bcftools.html#index). This index is necessary for creating the consensus genome using the compressed VCF.

### 15. VCF to consensus
Create consensus genome by applying VCF variants to the reference genome using [bcftools consensus](https://samtools.github.io/bcftools/bcftools.html#consensus).
This does not account for coverage, so we pass it the BEDfile of low coverage sites with the `--mask` option to mask these sites with `N`.

### 16. FASTA headers
Edits the FASTA headers to fit the pattern needed for downstream analysis.
Example FASTA header: `>SampleID|SampleID-PB2|H1N1pdm|PB2`
1. Replace reference sequence name with the NWGC sample ID using Perl to perform "lookaround" regex matches

### 17. Aggregate
This is the last rule of the pipeline that prints out the final result of each sample/reference pair.

The input of this rule differs based on the result of the checkpoint, so this rule dictates the final outcome of each sample/reference pair.
