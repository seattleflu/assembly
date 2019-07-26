# Seattle Flu Assembly Pipeline

## Introduction
Assembly pipeline used in the [Seattle Flu Study](https://seattleflu.org/) to build consensus genomes for sequenced samples.

Adapted from Louise Moncla's [Illumina pipline](https://github.com/lmoncla/illumina_pipeline) for influenza snp calling.

## Installation
Install dependencies via Conda by running:
```
conda env create -f envs/seattle-flu-environment.yaml
```
Once installed, the environment needs to be activated to use assembly by running:
```
conda activate seattle-flu
```

## Set Up
### FASTQ Files
The pipeline expects 2 or 8 total FASTQ files for each individual sample, with Read 1 (R1) and Read 2 (R2) from 1 lane or 4 different lanes.

Expected filename format: `318375_S12_L001_R1_001.fastq.gz` where `318375` is the NWGC sample ID.

If the filename does not contain the NWGC sample ID, rename the FASTQ files by using:
```
python scripts/rename_fastq.py
usage: rename_fastq.py [-h]
                       [Sample-ID-Map] [Barcode-column] [ID-column]
                       [FASTQ-directory]

Rename fastq files to contain NWGC sample ID

positional arguments:
  Sample-ID-Map    File path to csv file that matches random barcode to NWGC
                   ID (default: None)
  Barcode-column   Column name of column containing random barcodes (default:
                   None)
  ID-column        Column name of column containing NWGC sample IDs (default:
                   None)
  FASTQ-directory  File path to directory of FASTQ files that need to be
                   renamed (default: None)

optional arguments:
  -h, --help       show this help message and exit
```


### NWGC ID/SFS UUID key-value pair
Consensus genome FASTA headers are generated with the NWGC sample ID which then needs to be converted to the SFS UUID in order to be paired with stored metadata.

This is done with the [seqkit replace](https://bioinf.shenwei.me/seqkit/usage/#replace) method, which requires a tab-delimited key-value file for replacing the key with the value. The key-value file can be generated using the following:
```
python scripts/id_barcode_key_value.py
usage: id_barcode_key_value.py [-h]
                               [NWGC-SFS-Match] [NWGC-column] [SFS-column]
                               [Output]

Create key/value pairs of NWGC ID and SFS ID

positional arguments:
  NWGC-SFS-Match  File path to Excel file that matches NWGC ID to SFS ID
                  (default: None)
  NWGC-column     Column name of column containing NWGC sample IDs (default:
                  None)
  SFS-column      Column name of column containing SFS UUIDs (default: None)
  Output          File path for output key/value file (default: None)

optional arguments:
  -h, --help      show this help message and exit
```
__Note__: Remember to add the path of the output file to the config file after running this script.

## Batch Assembly Jobs
Currently, `config/config-flu-only.json` has the most updated version of the configuration needed to run batch assembly jobs. This means that there are several samples that are run at the same time.

* `fastq_directory`: the path to the directory containing the fastq files
* `ignored_samples`: an object containing keys that specify samples to be ignored
* `sample_reference_pairs`: an object containing specific sample/reference pairs, with the keys refering to the samples and an array of references as the values.
* `barcode_match`: the path to the file containing tab-delimited key-value pairs used to replace NWGC sample IDs with SFS UUIDs
* `min_align_rate`: the minimum overall align rate needed to generate a consensus genome (currently arbitrarily set to `1.00` so that `vcf_to_consensus` does not error out)
* `reference_virusus`: an object containing keys that specify references to be used
* `params`: parameters for the various tools used in the pipeline (currently set by Louise's recommendations)

**Setup Config File**
By default, the pipeline will generate combinations of all samples and references and try to create consensus genomes for all combinations.

Avoid this by specifying specific sample-reference pairs and ignored samples in the config file by using:
```
python scripts/setup_config_file.py
usage: setup_config_file.py [-h] [--sample [Sample column]]
                            [--target [Target column]]
                            [--max_difference [Max percent difference]]
                            [--lane [Lane]]
                            [FASTQ directory] [Config template] [Config file]
                            [Sample/target map]

Creates the Snakemake config file by copying the config template then adding
sample reference pairs and ignored samples. Currently need to download Excel
file from Metabase that contains the NWGC sample ID and the target identifiers
where present is True for the sample.

positional arguments:
  FASTQ directory       File path to directory containing fastq files
                        (default: None)
  Config template       File path to config file template used to generate
                        config files. (default: None)
  Config file           File path to config file to be used for snakemake.
                        (default: None)
  Sample/target map     File path to Excel file containing samples and present
                        targets (default: None)

optional arguments:
  -h, --help            show this help message and exit
  --sample [Sample column]
                        Column name of column containg NWGC sample IDs
                        (default: sample)
  --target [Target column]
                        Column name of column containing target identifiers
                        (default: target)
  --max_difference [Max percent difference]
                        The maximum difference acceptable between R1 and R2
                        reads (default: 0)
  --lane [Lane]         Specify the specific lane if samples only come from a
                        single lane. (default: None)
```
__Note__: This currently only works for flu-positive samples. To add more targets/references, edit `references/target_reference_map.json`.

__Check Read Pairs__
Embedded in this script is a check to ensure that the reads in R1 and R2 of each sample match. This is to avoid the demultiplexing error described in [this paper](https://www.nature.com/articles/s41588-019-0349-3). If the reads do not match 100% then the sample will be added to the ignored samples. Redirect the output to a file to keep track of read differences and ignored samples.

__Running Batch Jobs__:
```
snakemake --configfile config/config-flu-only.json -k
```
OR (if on Rhino)
```
snakemake -w 60 --configfile config/config-flu-only.json --cluster-config config/cluster.json --cluster "sbatch --nodes=1 --tasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} -o all_output.out" -j 20 -k
```

You can use `--use-conda` if there are rule-specific environment files.

## Single Sample/Target Jobs
Single sample/target pair jobs template is `config/one-sample-template.json`

* `fastq_files`: an array of absolute paths for FASTQ files of a single sample
* `sample`: the NWGC ID assigned to the sample, must be contained in the filename of the fastq_files
* `references`: an array of references
* `barcode_match`: the filepath to the key/value matching of NWGC ID to SFS UUID

**Setup Config File**
Single sample/target pair jobs require a different config file than the one for batch jobs:
```
usage: setup_config_for_one_sample_target_pair.py [-h] --config-template
                                                  <config_template.json>
                                                  --config-file <config.json>
                                                  --fastq-file
                                                  <sample.fastq.gz>
                                                  [<sample.fastq.gz> ...]
                                                  --nwgc-id [NWGC_ID]
                                                  --sfs-uuid [SFS_UUID]
                                                  --target [TARGET]
                                                  [--target-ref-map [TARGET_REF_MAP]]

Creates a config file for Snakemake for one sample/reference pair.

optional arguments:
  -h, --help            show this help message and exit
  --config-template <config_template.json>
                        File path to config file template (default: None)
  --config-file <config.json>
                        File path to config file to be used for snakemake
                        (default: None)
  --fastq-file <sample.fastq.gz> [<sample.fastq.gz> ...]
                        File path to fastq file. Expected format of fastq
                        filename: '346757_10DB_288_S281_L001_R1_001.fastq.gz
                        (default: None)
  --nwgc-id [NWGC_ID]   NWGC ID for the sample that matches ID in fastq
                        filename (default: None)
  --sfs-uuid [SFS_UUID]
                        The SFS UUID of the sample. (default: None)
  --target [TARGET]     The target to be used as reference in the assembly
                        pipeline (default: None)
  --target-ref-map [TARGET_REF_MAP]
                        The JSON file that contains the target/reference map
                        (default: references/target_reference_map.json)
```
__Check Read Pairs__
Embedded in the set up config script is a check that reads in R1 and R2 of each lane match 100%. If they do not match, then the script will error out.

__Running Single Sample/Target Jobs__:
```
snakemake --snakefile Snakefile_one_sample --configfile config/config-flu-only.json -k
```
OR (if on Rhino)
```
snakemake --snakefile Snakefile_one_sample -w 60 --configfile config/config-flu-only.json --cluster-config config/cluster.json --cluster "sbatch --nodes=1 --tasks=1 --mem={cluster.memory} --cpus-per-task={cluster.cores} --tmp={cluster.disk} --time={cluster.time} -o all_output.out" -j 20 -k
```

You can use `--use-conda` if there are rule-specific environment files.


## Basic Steps

### 1. Index reference genomes
Using [bowtie2-build](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#the-bowtie2-build-indexer) to build a Bowtie index for each reference genome. These will later be used in the mapping step.

### 2. Merge lanes
Concatenates 8 FASTQ files for each sample into 2 files (R1 and R2).

### 3. Trim fastqs
Trim raw FASTQ files with [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic), which cuts out adapter/illumina-specific sequences, trims the reads based on quality scores, and removes short reads.

### 4. Post trim fastqc
Generates summary statistics using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/). Allows for some quality control checks on raw sequence data coming from high throughput sequencing pipelines.

### 5. Map
Map trimmed reads to reference genome using [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml). Outputs a BAM file that represents aligned sequences and a log file that contains the alignment summary.

### 6. Sort
Sort the BAM file into "genome order" using [samtools sort](http://www.htslib.org/doc/samtools.html).

### 7. Bamstats
Generates statistics for coverage using [BAMStats](http://bamstats.sourceforge.net/).

### 8. Align rate checkpoint
[Checkpoints](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) allow for data-dependent conditional execution, forcing Snakemake to re-evaluate the DAG at this point.

Checks the alignment summary produced by the bowtie2 mapping. If the overall align rate is higher than the `min_align_rate` specified in the config, then continue with process to build consensus genome. If not, then stop the process for this sample/reference pair.

### 9. Not mapped
This rule is necessary for the checkpoint above to work, because it must send the data down one of two paths. This is the "dead end" where a consensus genome is __not__ generated.

Also a placeholder for how to handle sample/reference pairs that do not meet the `min_align_rate` threshold. Currently it just outputs the overall align rate of the sample/reference pair so that we can check that it did not meet the minimum.

### 10. Pileup
Generates [Pileup](https://en.wikipedia.org/wiki/Pileup_format) for BAM file using [samtools mpileup](http://www.htslib.org/doc/samtools.html).

  _Important Flags_:
  * `-a` to print out all positions, including zero depth (this is necessary for generating the low coverage Bedfile later)
  * `-A` to not discard anomalous read pairs
  * `-Q` to set the minimum base quality to consider a read (set to match the minimum in `varscan mpileup2snp` so that their coverage depths match)

### 11. Call SNPs
Calls SNPs from the Pileup based on parameters set in the config file using [varscan mpileup2snp](http://varscan.sourceforge.net/using-varscan.html#v2.3_mpileup2snp)

### 12. Zip VCF
Compress VCF using [bgzip](http://www.htslib.org/doc/bgzip.html), which allows indexes to be built against the file and allows the file be used without decompressing it.

### 13. Index BCF
Creates index for compressed VCF using [bcftools index](https://samtools.github.io/bcftools/bcftools.html#index). This index is necessary for creating the consensus genome using the compressed VCF.

### 14. VCF to consensus
Create consensus genome by applying VCF variants to the reference genome using [bcftools consensus](https://samtools.github.io/bcftools/bcftools.html#consensus). This does not account for coverage, so it will just fill in blanks with the base from the reference genome.

### 15. Create bed file
Creates a BED file for positions that need to be masked in the consensus genome. Positions need to be masked if they are below the minimum coverage depth and if they are disproportionally supported by one strand (>90%).

### 16. Mask consensus
Masks the consensus genome with "N" at bases where coverage is below the `min_cov` parameter using [bedtools maskfasta](https://bedtools.readthedocs.io/en/latest/content/tools/maskfasta.html).

### 17. FASTA headers
Edits the FASTA headers to fit the pattern needed for downstream analysis.
Example FASTA header: `>SFS-UUID|SFS-UUID-PB2|H1N1pdm|PB2`
1. Replace reference sequence name with the NWGC sample ID using Perl to perform "lookaround" regex matches
2. Uses [seqkit replace](https://bioinf.shenwei.me/seqkit/usage/#replace) to replace the NWGC sample ID with the SFS UUID.
3.  Create the UUID-gene combination using AWK

### 18. Combined FASTA
Creates the final combined FASTA file that will have all the consensus genomes generated in a day.

### 19. Aggregate
This is the last rule of the pipeline that prints out the final result of each sample/reference pair. If a consensus genome is generated, then this will also add it to the final combined FASTA.

The input of this rule differs based on the result of the checkpoint, so this rule dictates the final outcome of each sample/reference pair.

### Clean
Deletes all of the intermediate directories and files generated from running the pipeline. The only one that it does not delete is the `consensus_genome` directory.
