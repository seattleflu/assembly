"""Snakefile

Processes input fastq files to create consensus genomes and their associated
summary statistics for the Seattle Flu Study.

Install and activate the appropriate environment for the pipeline using
Conda (https://conda.io/en/latest/):
    $ conda env create -f envs/seattle-flu-environment.yaml
    $ conda activate seattle-flu

Before running perform a dry-run to test that all input files are correctly
located and that the Snakemake DAG builds correctly:
    $ snakemake -n

To run on a local machine:
    $ snakemake -k

To run on the Fred Hutch cluster (Rhino):
    $ snakemake \
        -w 60 \
        --cluster-config config/cluster.json \
        --cluster "sbatch \
            --nodes=1 \
            --tasks=1 \
            --mem={cluster.memory} \
            --cpus-per-task={cluster.cores} \
            --tmp={cluster.disk} \
            --time={cluster.time} \
            -o all_output.out" \
        -j 20 \
        -k

Basic steps:
    1. Trim raw fastq's with Trimmomatic
    2. Map trimmed reads to each reference genomes in the reference panel using
       bowtie2 # This step may change with time
    ~3. Remove duplicate reads using Picard~
    4. Call SNPs using varscan
    5. Use SNPs to generate full consensus genomes for each sample x reference
       virus combination
    6. Compute summary statistics for each sample x refernce virus combination

Adapted from Louise Moncla's illumina pipeline for influenza snp calling:
https://github.com/lmoncla/illumina_pipeline
"""
import sys, os
import json
import glob
import getpass
import requests
from requests.exceptions import HTTPError
from urllib.parse import urljoin
from datetime import datetime
from itertools import product

#### Helper functions and variable def'ns
def generate_sample_ids(cfg):
    """Return all the sample names (i.e. S04) as a list.
    Will ignore all samples listed under "ignored_samples" of config.
    """
    all_ids = set()
    for f in glob.glob("{}/*".format(cfg['fastq_directory'])):
        if f.endswith('.fastq.gz'):
            f = f.split('.')[0].split('/')[-1].split('_')[0]
            if f in cfg['ignored_samples']:
                continue
            all_ids.add(f)
    return all_ids

def generate_all_files(sample, config):
    """For a given sample, determine all fastqs as a tuple of forward
    and reverse
    """
    r1 = glob.glob('{}/*{}*R1*'.format(config['fastq_directory'], sample))
    r2 = glob.glob('{}/*{}*R2*'.format(config['fastq_directory'], sample))
    return (r1, r2)

def filter_combinator(combinator, config):
    """ Custom combinatoric function to be used with expand function.
    Only generates combination of sample and reference wildcards that are
    listed under "sample_reference_pairs" of config file.
    If no references are listed for a sample, will generate all possible
    combinations with all references.
    """
    def filtered_combinator(*args, **kwargs):
        for wc_comb in combinator(*args, **kwargs):
            if wc_comb[0][1] not in config["sample_reference_pairs"]:
                yield wc_comb
            elif len(config["sample_reference_pairs"][wc_comb[0][1]]) == 0:
                yield wc_comb
            elif wc_comb[1][1] in config["sample_reference_pairs"][wc_comb[0][1]]:
                yield wc_comb
    return filtered_combinator

# input function for the rule aggregate
def aggregate_input(wildcards):
    """
    Returns input for rule aggregate based on output from checkpoint align_rate.
    Set minimum align rate in config under "min_align_rate".
    """
    with open(checkpoints.mapped_reads.get(sample=wildcards.sample, reference=wildcards.reference).output[0]) as f:
        summary = json.load(f)
        all_segments_aligned = summary["all_segments_aligned"]
        min_reads = summary["minimum_reads_required"]
        mapped = summary["mapped_reads"]

        if not all_segments_aligned or mapped <= min_reads:
            return rules.not_mapped.output.not_mapped
        else:
            return rules.fasta_headers.output.masked_consensus

# Build static lists of reference genomes, ID's of samples, and ID -> input files
all_references = [ v for v in  config['reference_viruses'].keys() ]
all_ids = generate_sample_ids(config)
mapped = {id: generate_all_files(id, config) for id in all_ids}
filtered_product = filter_combinator(product, config)

#### Main pipeline
rule all:
    input:
        # pre_fastqc = expand("summary/pre_trim_fastqc/{fname}_fastqc.html",
        #        fname=glob.glob(config['fastq_directory']),
        post_fastqc = expand("summary/post_trim_fastqc/{sample}.trimmed_{tr}_fastqc.{ext}",
               sample=all_ids,
               tr=["1P", "1U", "2P", "2U"],
               ext=["zip", "html"]),
        bamstats = expand("summary/bamstats/{reference}/{sample}.coverage_stats.txt", filtered_product,
               sample=all_ids,
               reference=all_references),
        aggregate = expand("summary/aggregate/{reference}/{sample}.log", filtered_product,
                sample=all_ids,
                reference=all_references),

rule index_reference_genome:
    input:
        raw_reference = "references/{reference}.fasta"
    output:
        indexed_reference = "references/{reference}.1.bt2"
    # group:
    #     "pre-mapping"
    shell:
        """
        bowtie2-build {input.raw_reference} references/{wildcards.reference}
        """

rule merge_lanes:
    input:
        all_r1 = lambda wildcards: mapped[wildcards.sample][0],
        all_r2 = lambda wildcards: mapped[wildcards.sample][1]
    output:
        merged_r1 = "process/merged/{sample}_R1.fastq.gz",
        merged_r2 = "process/merged/{sample}_R2.fastq.gz"
    # group:
    #     "pre-mapping"
    shell:
        """
        cat {input.all_r1} >> {output.merged_r1}
        cat {input.all_r2} >> {output.merged_r2}
        """
# Ignored for the time being
# rule pre_trim_fastqc:
#     input:
#         all_fastq = lambda wildcards: mapped[wildcards.sample][0]+mapped[wildcards.sample][1]
#     output:
#         qc = "summary/pre_trim_fastqc/{sample}_fastqc.html",
#         zip = "summary/pre_trim_fastqc/{sample}_fastqc.zip"
#     shell:
#         """
#         fastqc {input.all_fastq} -o summary/pre_trim_fastqc
#         """

rule trim_fastqs:
    input:
        fastq_f = "process/merged/{sample}_R1.fastq.gz",
        fastq_r = "process/merged/{sample}_R2.fastq.gz"
    output:
        trimmed_fastq_1p = "process/trimmed/{sample}.trimmed_1P.fastq",
        trimmed_fastq_2p = "process/trimmed/{sample}.trimmed_2P.fastq",
        trimmed_fastq_1u = "process/trimmed/{sample}.trimmed_1U.fastq",
        trimmed_fastq_2u = "process/trimmed/{sample}.trimmed_2U.fastq"
    params:
        paired_end = config["params"]["trimmomatic"]["paired_end"],
        adapters = config["params"]["trimmomatic"]["adapters"],
        illumina_clip = config["params"]["trimmomatic"]["illumina_clip"],
        window_size = config["params"]["trimmomatic"]["window_size"],
        trim_qscore = config["params"]["trimmomatic"]["trim_qscore"],
        minimum_length = config["params"]["trimmomatic"]["minimum_length"]
    benchmark:
        "benchmarks/{sample}.trimmo"
    # group:
    #     "trim-and-map"
    shell:
        """
        trimmomatic \
            {params.paired_end} \
            -phred33 \
            {input.fastq_f} {input.fastq_r} \
            -baseout process/trimmed/{wildcards.sample}.trimmed.fastq \
            ILLUMINACLIP:{params.adapters}:{params.illumina_clip} \
            SLIDINGWINDOW:{params.window_size}:{params.trim_qscore} \
            MINLEN:{params.minimum_length}
        """

rule post_trim_fastqc:
    input:
        p1 = rules.trim_fastqs.output.trimmed_fastq_1p,
        p2 = rules.trim_fastqs.output.trimmed_fastq_2p,
        u1 = rules.trim_fastqs.output.trimmed_fastq_1u,
        u2 = rules.trim_fastqs.output.trimmed_fastq_2u
    output:
        qc_1p_html = "summary/post_trim_fastqc/{sample}.trimmed_1P_fastqc.html",
        qc_2p_html = "summary/post_trim_fastqc/{sample}.trimmed_2P_fastqc.html",
        qc_1u_html = "summary/post_trim_fastqc/{sample}.trimmed_1U_fastqc.html",
        qc_2u_html = "summary/post_trim_fastqc/{sample}.trimmed_2U_fastqc.html",
        qc_1p_zip = "summary/post_trim_fastqc/{sample}.trimmed_1P_fastqc.zip",
        qc_2p_zip = "summary/post_trim_fastqc/{sample}.trimmed_2P_fastqc.zip",
        qc_1u_zip = "summary/post_trim_fastqc/{sample}.trimmed_1U_fastqc.zip",
        qc_2u_zip = "summary/post_trim_fastqc/{sample}.trimmed_2U_fastqc.zip"
    # group:
    #     "summary-statistics"
    shell:
        """
        fastqc {input.p1} -o summary/post_trim_fastqc
        fastqc {input.p2} -o summary/post_trim_fastqc
        fastqc {input.u1} -o summary/post_trim_fastqc
        fastqc {input.u2} -o summary/post_trim_fastqc
        """

rule map:
    input:
        p1 = rules.trim_fastqs.output.trimmed_fastq_1p,
        p2 = rules.trim_fastqs.output.trimmed_fastq_2p,
        u1 = rules.trim_fastqs.output.trimmed_fastq_1u,
        u2 = rules.trim_fastqs.output.trimmed_fastq_2u,
        ref_file = rules.index_reference_genome.output.indexed_reference
    output:
        mapped_bam_file = "process/mapped/{reference}/{sample}.bam",
        bt2_log = "summary/bowtie2/{reference}/{sample}.log"
    params:
        threads = config["params"]["bowtie2"]["threads"],
        map_all = config["params"]["bowtie2"]["all"]
    benchmark:
        "benchmarks/{sample}_{reference}.bowtie2"
    # group:
    #     "trim-and-map"
    shell:
        """
        bowtie2 \
            -x references/{wildcards.reference} \
            -1 {input.p1} \
            -2 {input.p2} \
            -U {input.u1},{input.u2} \
            -P {params.threads} \
            {params.map_all} \
            --local 2> {output.bt2_log} | \
                samtools view -bSF4 - > {output.mapped_bam_file}
        """

rule sort:
    input:
        mapped_bam = rules.map.output.mapped_bam_file
    output:
        sorted_bam_file = "process/sorted/{reference}/{sample}.sorted.bam"
    benchmark:
        "benchmarks/{sample}_{reference}.sort"
    # group:
    #     "post-mapping"
    shell:
        """
        samtools view \
            -bS {input.mapped_bam} | \
            samtools sort | \
            samtools view -h > {output.sorted_bam_file}
        """

rule bamstats:
    input:
        sorted_bam = rules.sort.output.sorted_bam_file
    output:
        bamstats_file = "summary/bamstats/{reference}/{sample}.coverage_stats.txt"
    benchmark:
        "benchmarks/{sample}_{reference}.bamstats"
    # group:
    #     "summary-statistics"
    shell:
        """
        BAMStats -i {input.sorted_bam} > {output.bamstats_file}
        """

# Removed Picard because it was too computationally intensive.
# Uncomment and fix inputs of pileup to re-add
# rule remove_duplicate_reads:
#     input:
#         sorted_bam = rules.sort.output.sorted_bam_file
#     output:
#         deduped = "process/deduped/{reference}/{sample}.nodups.sam"
#     params:
#         picard_params = "summary/file.params.txt"
#     shell:
#         """
#         picard \
#             MarkDuplicates \
#             I={input.sorted_bam} \
#             O={output.deduped} \
#             REMOVE_DUPLICATES=true \
#             M={params.picard_params}
#         """

checkpoint mapped_reads:
    input:
        bt2_log = rules.map.output.bt2_log,
        bamstats = rules.bamstats.output.bamstats_file,
        reference = "references/{reference}.fasta"
    output:
        "summary/checkpoint/{reference}/{sample}.json"
    params:
        min_cov = config["params"]["varscan"]["min_cov"],
        raw_read_length = config["raw_read_length"]
    shell:
        """
        python scripts/checkpoint_mapped_reads.py \
            --bamstats {input.bamstats} \
            --bowtie2 {input.bt2_log} \
            --min-cov {params.min_cov} \
            --raw-read-length {params.raw_read_length} \
            --reference {input.reference} \
            --output {output} \
        """

rule not_mapped:
    input:
        sorted_bam = rules.sort.output.sorted_bam_file,
        reference = "references/{reference}.fasta"
    output:
        not_mapped = "summary/not_mapped/{reference}/{sample}.txt"
    shell:
        """
        touch {output.not_mapped}
        """

rule pileup:
    input:
        sorted_bam = rules.sort.output.sorted_bam_file,
        reference = "references/{reference}.fasta"
    output:
        pileup = "process/mpileup/{reference}/{sample}.pileup"
    params:
        depth = config["params"]["mpileup"]["depth"],
        min_base_qual = config["params"]["varscan"]["snp_qual_threshold"]
    benchmark:
        "benchmarks/{sample}_{reference}.mpileup"
    # group:
    #     "post-mapping"
    shell:
        """
        samtools mpileup -a -A \
            -Q {params.min_base_qual} \
            -d {params.depth} \
            {input.sorted_bam} > {output.pileup} \
            -f {input.reference}
        """

rule call_snps:
    input:
        pileup = rules.pileup.output.pileup
    output:
        vcf = "process/vcfs/{reference}/{sample}.vcf"
    params:
        min_cov = config["params"]["varscan"]["min_cov"],
        snp_qual_threshold = config["params"]["varscan"]["snp_qual_threshold"],
        snp_frequency = config["params"]["varscan"]["snp_frequency"]
    benchmark:
        "benchmarks/{sample}_{reference}.varscan"
    shell:
        """
        varscan mpileup2snp \
            {input.pileup} \
            --min-coverage {params.min_cov} \
            --min-avg-qual {params.snp_qual_threshold} \
            --min-var-freq {params.snp_frequency} \
            --strand-filter 1 \
            --output-vcf 1 > {output.vcf}
        """

rule zip_vcf:
    input:
        vcf = rules.call_snps.output.vcf
    output:
        bcf = "process/vcfs/{reference}/{sample}.vcf.gz"
    shell:
        """
        bgzip {input.vcf}
        """

rule index_bcf:
    input:
        bcf = rules.zip_vcf.output.bcf
    output:
        index = "process/vcfs/{reference}/{sample}.vcf.gz.csi"
    shell:
        """
        bcftools index {input}
        """

rule vcf_to_consensus:
    input:
        bcf = rules.zip_vcf.output.bcf,
        index = rules.index_bcf.output.index,
        ref = "references/{reference}.fasta"
    output:
        consensus_genome = "consensus_genomes/{reference}/{sample}.consensus.fasta"
    benchmark:
        "benchmarks/{sample}_{reference}.consensus"
    shell:
        """
        cat {input.ref} | \
            bcftools consensus {input.bcf} > \
            {output.consensus_genome}
        """

rule create_bed_file:
    input:
        pileup = rules.pileup.output.pileup
    output:
        bed_file = "summary/low_coverage/{reference}/{sample}.bed"
    params:
        min_cov = config["params"]["varscan"]["min_cov"],
        min_freq = config["params"]["varscan"]["snp_frequency"]
    shell:
        """
        python scripts/create_bed_file_for_masking.py \
            --pileup {input.pileup} \
            --min-cov {params.min_cov} \
            --min-freq {params.min_freq} \
            --bed-file {output.bed_file}
        """

rule mask_consensus:
    input:
        consensus_genome = rules.vcf_to_consensus.output.consensus_genome,
        low_coverage = rules.create_bed_file.output.bed_file
    output:
        masked_consensus = temp("consensus_genomes/{reference}/{sample}.temp_consensus.fasta")
    shell:
        """
        bedtools maskfasta \
            -fi {input.consensus_genome} \
            -bed {input.low_coverage} \
            -fo {output.masked_consensus}
        """


rule fasta_headers:
    input:
        masked_consensus = rules.mask_consensus.output
    output:
        masked_consensus = "consensus_genomes/{reference}/{sample}.masked_consensus.fasta"
    shell:
        """
        cat {input.masked_consensus} | \
            perl -pi -e 's/(?<=>)[^>|]*(?<=|)/{wildcards.sample}/g' > \
            {output.masked_consensus}.temp
        awk '{{split(substr($0,2),a,"|"); \
            if(a[2]) print ">"a[1]"|"a[1]"-"a[2]"-"a[3]"|"a[2]"|"a[3]; \
            else print; }}' \
            {output.masked_consensus}.temp > {output.masked_consensus}
        rm {output.masked_consensus}.temp
        """

rule combined_fasta:
    output:
        combined_fasta = "consensus_genomes/batch-{current_datetime}.fasta"
                          .format(current_datetime=datetime.now().date())
    shell:
        """
        touch {output.combined_fasta}
        """

rule aggregate:
    input:
        aggregate_input = aggregate_input
    output:
        aggregate_summary = "summary/aggregate/{reference}/{sample}.log"
    params:
        combined_fasta = rules.combined_fasta.output
    run:
        if input.aggregate_input.split(".")[-1] == "fasta":
            shell("cat {input.aggregate_input} >> {params.combined_fasta}")
        shell("echo 'Final output: {input.aggregate_input}' > {output}")

rule clean:
    message: "Removing directories: {params}"
    params:
        "summary "
        "process "
        "benchmarks "
        "references/*.bt2 "
        "references/*.fai"
    shell:
        """
        rm -rfv {params}
        """
