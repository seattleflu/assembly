"""
Snakefile specific for running assembly pipeline on Seattle Flu Study (SFS)
samples that automatically uploads resulting genomes to the SFS ID3C database.

To run:
    $ snakemake -k --snakefile Snakefile --configfile <config_file.json>
"""
import os
import json
import getpass
import requests
from requests.exceptions import HTTPError
from urllib.parse import urljoin


# Include the main assembly pipeline
include: "Snakefile-base"

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
            return rules.post_masked_consensus_and_summary_stats_to_id3c.output.successful_post


rule all:
    input:
        # pre_fastqc = expand("summary/pre_trim_fastqc/{fname}_fastqc.html",
        #        fname=glob.glob(config['fastq_directory']),
        #post_fastqc = expand("summary/post_trim_fastqc/{sample}.trimmed_{tr}_fastqc.{ext}",
        #       sample=all_ids,
        #       tr=["1P", "1U", "2P", "2U"],
        #       ext=["zip", "html"]),
        bamstats = expand("summary/bamstats/{reference}/{sample}.coverage_stats.txt", filtered_product,
               sample=all_ids,
               reference=all_references),
        aggregate = expand("summary/aggregate/{reference}/{sample}.log", filtered_product,
                sample=all_ids,
                reference=all_references),
        consensus_genomes = expand("consensus_genomes/{reference}/{sample}.consensus.fasta", filtered_product,
            sample=all_ids,
            reference=all_references
        ),


rule nwgc_sfs_map:
    params:
        mapper_file = config["barcode_match"]["mapper_filepath"],
        nwgc_column = config["barcode_match"]["nwgc_column"],
        sfs_column = config["barcode_match"]["sfs_column"],
    output:
        key_value_file = config["barcode_match"]["key_value_filepath"]
    shell:
        """
        python3 scripts/id_barcode_key_value.py \
            {params.mapper_file} \
            {params.nwgc_column} \
            {params.sfs_column} \
            {output.key_value_file}
        """

rule fasta_headers:
    input:
        key_value_file = rules.nwgc_sfs_map.output,
        masked_consensus = rules.vcf_to_consensus.output
    output:
        masked_consensus = "consensus_genomes/{reference}/{sample}.masked_consensus.fasta"
    shell:
        """
        cat {input.masked_consensus} | \
            perl -pi -e 's/(?<=>)[^>|]*(?<=|)/{wildcards.sample}/g' > \
            {output.masked_consensus}.temp
        seqkit replace -p '({wildcards.sample})' -r '{{kv}}' \
            -k {input.key_value_file} --keep-key \
            {output.masked_consensus}.temp > {output.masked_consensus}
        awk '{{split(substr($0,2),a,"|"); \
            if(a[2]) print ">"a[1]"|"a[1]"-"a[2]"-"a[3]"|"a[2]"|"a[3]; \
            else print; }}' \
            {output.masked_consensus} > {output.masked_consensus}.temp
        mv {output.masked_consensus}.temp {output.masked_consensus}

        """

rule metadata_to_json:
    input:
        all_r1 = lambda wildcards: mapped[wildcards.sample][0],
        all_r2 = lambda wildcards: mapped[wildcards.sample][1]
    output:
        temp("consensus_genomes/{reference}/{sample}.metadata.json")
    shell:
        """
        python scripts/metadata_to_json.py {input.all_r1:q} {input.all_r2:q} > {output}
        """

rule masked_consensus_to_json:
    input:
        masked_consensus = rules.fasta_headers.output.masked_consensus
    output:
        temp("consensus_genomes/{reference}/{sample}.masked_consensus.json")
    shell:
        """
        python scripts/fasta_to_json.py {input.masked_consensus} > {output}
        """

rule summary_stats_to_json:
    input:
        bam_coverage = rules.bamstats.output.bamstats_file,
        bowtie2 = rules.map.output.bt2_log
    output:
        temp("consensus_genomes/{reference}/{sample}.summary_stats.json")
    shell:
        """
        python scripts/summary_stats_to_json.py --bamstats {input.bam_coverage} \
            --bowtie2 {input.bowtie2} > {output}
        """

rule create_id3c_payload:
    input:
        metadata = rules.metadata_to_json.output,
        masked_consensus = rules.masked_consensus_to_json.output,
        summary_stats = rules.summary_stats_to_json.output
    params:
        status = 'complete'
    output:
        "consensus_genomes/{reference}/{sample}.payload.json"
    shell:
        """
        python scripts/create_id3c_payload.py \
            --masked-consensus {input.masked_consensus} \
            --summary-stats {input.summary_stats} \
            --metadata {input.metadata} \
            --status {params.status} > {output}
        """

rule post_masked_consensus_and_summary_stats_to_id3c:
    input:
        rules.create_id3c_payload.output
    output:
        successful_post = "consensus_genomes/{reference}/{sample}.successful-post.log"
    params:
        id3c_url = os.environ['ID3C_URL'],
        id3c_username = os.environ['ID3C_USERNAME'],
        id3c_password = os.environ['ID3C_PASSWORD'],
        id3c_slack_webhook = os.environ['SLACK_WEBHOOK_URL'],
    log: "consensus_genomes/{reference}/{sample}.http-response.log"
    run:
        headers = {'Content-type': 'application/json'}
        file = open(str(log), "w")

        with open(str(input)) as f:
            data = f.read()

        try:
            response = requests.post(
                urljoin(params.id3c_url, 'v1/receiving/consensus-genome'),
                data=data,
                headers=headers,
                auth=(params.id3c_username, params.id3c_password))

            response.raise_for_status()

            if response.ok:
                with open(str(output), "w"):
                    pass

        except HTTPError as http_err:
            file.write(str(http_err))

            slack_data = { "text":
                f":rotating_light: Hey {getpass.getuser()}: Assembly failed to upload to ID3C with HTTP status code: " +
                f"{http_err.response.status_code}.\nMore details at `{log}`"
            }

            try:
                slack_response = requests.post(params.id3c_slack_webhook,
                    data=json.dumps(slack_data), headers=headers)

                slack_response.raise_for_status()

            except HTTPError as slack_http_err:
                file.write(str(slack_http_err))

            raise http_err

        except Exception as err:
            file.write(str(err))

            raise Exception(f"Error: {err} in ID3C POST request.")

        finally:
            file.close()


rule aggregate:
    input:
        aggregate_input = aggregate_input
    output:
        aggregate_summary = "summary/aggregate/{reference}/{sample}.log"
    run:
        shell("echo 'Final output: {input.aggregate_input}' > {output}")
