"""Snakefile

Processes input fastq files to create consensus genomes and their associated summary statistics for the Seattle Flu Study.

To run:
    $ snakemake

Basic steps:
    1. Trim raw fastq's with Trimmomatic
    2. Map trimmed reads to each reference genomes in the reference panel using bowtie2 # This step may change with time
    3. Remove duplicate reads using Picard
    4. Call SNPs using varscan
    5. Use SNPs to generate full consensus genomes for each sample x reference virus combination
    6. Compute summary statistics for each sample x refernce virus combination

Adapted from Louise Moncla's illumina pipeline for influenza snp calling:
https://github.com/lmoncla/illumina_pipeline
"""
import sys, os
import glob
configfile: "config.json"

#### Helper functions and variable def'ns
def generate_sample_names(cfg):
    sample_names = []
    for f in glob.glob("{}/*".format(cfg['fastq_directory'])):
        if f.endswith('.fastq.gz'):
            f = f.split('.')[0].split('/')[-1]
            sample_names.append(f)
    return sample_names

all_sample_names = generate_sample_names(config)
all_references = [ v for v in  config['reference_viruses'].keys() ]


#### Main pipeline
rule all:
    input:
        consensus_genome = expand("consensus_genomes/{reference}/{sample}.consensus.fasta",
               sample=all_sample_names,
               reference=all_references)
        # summary_statistics = ()

rule trim_fastqs:
    input:
        fastq = "test_data/{sample}.fastq.gz"
    output:
        trimmed_fastq = "process/trimmed/{reference}/{sample}.trimmed.fastq"
    params:
        paired_end = config["params"]["trimmomatic"]["paired_end"],
        adapters = config["params"]["trimmomatic"]["adapters"],
        illumina_clip = config["params"]["trimmomatic"]["illumina_clip"],
        window_size = config["params"]["trimmomatic"]["window_size"],
        trim_qscore = config["params"]["trimmomatic"]["trim_qscore"],
        minimum_length = config["params"]["trimmomatic"]["minimum_length"]
    shell:
        """
        trimmomatic \
            {params.paired_end} \
            -phred33 \
            {input.fastq} \
            {output.trimmed_fastq} \
            ILLUMINACLIP:{params.adapters}:{params.illumina_clip} \
            SLIDINGWINDOW:{params.window_size}:{params.trim_qscore} \
            MINLEN:{params.minimum_length}
        """

rule map:
    input:
        fastq = rules.trim_fastqs.output.trimmed_fastq
    output:
        mapped_sam_file = "process/mapped/{reference}/{sample}.sam"
    shell:
        """
        bowtie2 \
            -x references/{wildcards.reference} \
            -U {input.fastq} \
            -S {output.mapped_sam_file} \
            --local
        """

rule sort:
    input:
        mapped_sam = rules.map.output.mapped_sam_file
    output:
        sorted_sam_file = "process/sorted/{reference}/{sample}.sorted.sam"
    shell:
        """
        samtools view \
            -bS {input.mapped_sam} | \
            samtools sort | \
            samtools view -h > {output.sorted_sam_file}
        """

rule remove_duplicate_reads:
    input:
        sorted_sam = rules.sort.output.sorted_sam_file
    output:
        deduped = "process/deduped/{reference}/{sample}.nodups.sam"
    params:
        picard_params = "file.params.txt"
    shell:
        """
        picard \
            MarkDuplicates \
            I={input.sorted_sam} \
            O={output.deduped} \
            REMOVE_DUPLICATES=true \
            M={params.picard_params}
        """

rule call_snps:
    input:
        deduped_sam = rules.remove_duplicate_reads.output.deduped,
        reference = "references/{reference}.fasta"
    output:
        vcf = "process/vcfs/{reference}/{sample}.vcf"
    params:
        depth = "1000000",
        min_cov = "",
        snp_qual_threshold = "",
        snp_frequency = "",
    shell:
        """
        samtools mpileup \
            -d {params.depth} \
            {input.deduped_sam} > process/tmp.pileup \
            -f {input.reference}
        java -jar /usr/local/bin/VarScan.v2.3.9.jar mpileup2snp \
            process/tmp.pileup \
            --min-coverage {params.min_cov} \
            --min-avg-qual {params.snp_qual_threshold} \
            --min-var-freq {params.snp_frequency} \
            --strand-filter 1 \
            --output-vcf 1 > {output.vcf}
        """

rule vcf_to_consensus:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        consensus_genome = "consensus_genomes/{reference}/{sample}.consensus.fasta"
    shell:
        """
        touch {output.consensus_genome}
        """
