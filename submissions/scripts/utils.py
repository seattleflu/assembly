
import conda.cli.python_api as Conda
import logging
import os
import docker
import sys
import requests
import pandas as pd
import getpass
from pathlib import Path

LOG_LEVEL = os.environ.get("LOG_LEVEL", "debug").upper()
logging.captureWarnings(True)

LOG = logging.getLogger(__name__)
LOG.setLevel(LOG_LEVEL)


def yes_no_cancel(message: str) -> bool:
    while True:
        user_input = input(f"{message} [y]es/[n]o/[c]ancel: ").lower()
        if user_input not in ['y','n','c']:
            print("Not a valid response")
            continue
        else:
            break
    if user_input == 'y':
        return True
    elif user_input == 'n':
        return False
    else:
        print("Cancelling.")
        sys.exit()


def process_with_nextclade(nextclade_dataset_dir:Path, output_tsv:Path, input_fasta:Path):
    # Get data files from NextClade
    result = Conda.run_command('run', 'nextclade', 'dataset', 'get',
        "--name=sars-cov-2",
        f"--output-dir={nextclade_dataset_dir}"
    )

    if result and len(result)==3 and result[2] == 0:
        LOG.debug("Successfully downloaded data/sars-cov-2 from NextClade.")
    else:
        raise Exception('Error: Nexclade sars-cov-2 data download failed.')


    LOG.debug("Analyzing FASTA file with NextClade.")
    result = Conda.run_command('run', 'nextclade', 'run',
        f"--input-dataset={nextclade_dataset_dir}",
        f"--output-tsv={output_tsv}",
        f"{input_fasta}"
    )

    if result[2] == 0:
        LOG.debug(f"NextClade processing complete: {output_tsv}")
    else:
        raise Exception(f"Error: NextClade processing of {input_fasta} failed:\n {result}")


def pull_previous_submissions(output_tsv:Path) -> Path:
    tsv_url = "https://raw.github.com/seattleflu/hcov19-sequence-identifiers/master/hcov19-sequence-identifiers.tsv"
    LOG.debug(f"Downloading {tsv_url}")

    token = os.environ.get('GH_ACCESS_TOKEN') or getpass.getpass('Github Personal access token:')
    try:
        r = requests.get('https://api.github.com/repos/seattleflu/hcov19-sequence-identifiers/contents/hcov19-sequence-identifiers.tsv',
            headers={
                'accept': 'application/vnd.github.v3.raw',
                'authorization': 'token {}'.format(token)
            }
        )

        if r.status_code == 200:
            with open(output_tsv, 'w') as f:
                f.write(r.text)
            return output_tsv
        else:
            LOG.error(f"Error downloading previous submissions file: {r.status_code}")
            raise Exception

    except requests.exceptions.HTTPError as e:
        LOG.error(f"Error: {e.response.text}")
        raise Exception


def calculate_excluded_vocs(nextclade_df: pd.DataFrame, output_csv:str, identifier_format:str='nwgc_id'):

    if identifier_format == 'nwgc_id':
        nextclade_df["LIMS"] = nextclade_df["seqName"].str.extract("^[^\d]*(\d+)").astype(int)
    elif identifier_format == 'sfs_sample_barcode':
        nextclade_df["LIMS"] = nextclade_df["seqName"].str.split('|').str[0]

    nextclade_df["aaSubstitutions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
    nextclade_df["aaInsertions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
    nextclade_df["aaDeletions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
    nextclade_df["aaChanges_s_total"] = nextclade_df["aaSubstitutions_s_count"] + \
                                        nextclade_df["aaInsertions_s_count"] + \
                                        nextclade_df["aaDeletions_s_count"]

    excluded_vocs = nextclade_df[
        (nextclade_df["errors"].notnull()) |
        (nextclade_df["clade"].isin(["19A","19B","recombinant"])) |
        ("S" in nextclade_df["failedGenes"].astype(str).str.split(",")) |
        (nextclade_df["aaChanges_s_total"] < 3) |
        (nextclade_df["qc.frameShifts.frameShifts"].astype(str).str.split(",").apply(count_s_occurrences) > 0)]

    # remove positive controls from excluded VOCs
    if identifier_format == 'sfs_sample_barcode':
        excluded_vocs = excluded_vocs[~excluded_vocs["LIMS"].astype(str).str.contains('-pos-con')]

    excluded_vocs[['LIMS']].to_csv(output_csv, header=False, index=False)

    LOG.debug(f'Excluded VOCs:{excluded_vocs[["LIMS", "seqName", "clade", "qc.missingData.status", "aaChanges_s_total", "qc.frameShifts.frameShifts", "failedGenes", "errors"]]}')
    LOG.debug(f"Saved {len(excluded_vocs)} NWGC ids to {output_csv}.")


def count_s_occurrences(values:str) -> int:
    count = 0
    if not isinstance(values, list):
        return count
    for item in values:
        if isinstance(item, str) and item.lower().strip().startswith("s:"):
            count += 1
    return count


def process_with_vadr(input_fasta:Path, output_batch_dir:Path, interactive:bool=False) -> str:
    # process with VADR
    process_with_vadr = None
    docker_client = None
    while True:
        process_with_vadr = (interactive==False) or yes_no_cancel("Process with VADR? (y/n)")
        if process_with_vadr:
            try:
                docker_client = docker.from_env()
                break
            except:
                LOG.warning("Docker client could not be initiated. Make sure Docker is running.")
                continue
        else:
            break

    if process_with_vadr and docker_client:
        # pull latest docker image
        docker_client.images.pull('staphb/vadr')

        # run vadr
        docker_container = docker_client.containers.run('staphb/vadr',
            command = f'/bin/bash -c \
                "/opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
                    --minlen 50 --maxlen 30000 \
                    /data/{input_fasta.name} > trimmed-genbank.fasta; \
                v-annotate.pl --split --cpu 8 --glsearch -f -s -r \
                    --noseqnamemax \
                    --nomisc --mkey sarscov2 \
                    --lowsim5seq 6 --lowsim3seq 6 \
                    --alt_fail lowscore,insertnn,deletinn \
                    --mdir /opt/vadr/vadr-models/ \
                    /data/trimmed-genbank.fasta \
                    /data/genbank"',
            detach=True,
            remove=True,
            user=f"{os.getuid()}:{os.getgid()}",
            volumes= {output_batch_dir.resolve(): {'bind': '/data', 'mode': 'rw'}}
        )
        docker_output = docker_container.attach(stdout=True, stream=True, logs=True);
        for line in docker_output:
            LOG.debug(line.decode('utf-8'))

        LOG.debug(f"Finished VADR processing, output saved to: {Path(output_batch_dir, 'genbank')}")
        return Path(output_batch_dir, 'genbank')


    raise Exception(f"Error: Something went wrong with VADR processing of {input_fasta.resolve()}")



def sequence_metrics(sequence: str, fasta_id: str) -> dict:
    """
    Calculate basic metrics from sars-cov-2 sequence
    """
    n_count = sequence.count("N")
    seq_length = len(sequence)
    ref_length = 29903  # Wuhan-Hu-1 complete genome
    metrics = {
        'CONTIG_NAME': fasta_id,
        'COUNT_A': sequence.count("A"),
        'COUNT_C': sequence.count("C"),
        'COUNT_G': sequence.count("G"),
        'COUNT_T': sequence.count("T"),
        'COUNT_N': n_count,
        'CONSENSUS_FASTA_LENGTH': len(sequence),
        'PCT_N_MAPPED': n_count / seq_length,
        'REFERENCE_LENGTH': ref_length,
        'PCT_N_REFERENCE': n_count / ref_length,
    }
    return metrics
