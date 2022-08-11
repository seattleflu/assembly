"""
Prepares files and directories for sequencing submissions.

"""
import os
import shutil
import argparse
import tarfile
from pathlib import Path
from datetime import datetime
import tempfile
import logging
import pandas as pd
import conda.cli.python_api as Conda
import docker

LOG_LEVEL = os.environ.get("LOG_LEVEL", "debug").upper()

logging.basicConfig(
    level = logging.ERROR,
    format = "[%(asctime)s] %(levelname)-8s %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S%z")

logging.captureWarnings(True)

LOG = logging.getLogger(__name__)
LOG.setLevel(LOG_LEVEL)

def count_s_occurrences(values):
    count = 0
    if not isinstance(values, list):
        return count
    for item in values:
        if isinstance(item, str) and item.lower().strip().startswith("s:"):
            count += 1
    return count

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--results-file", required=True,
        help = "File path to the assembly results compressed tar file (.tar.gz)")
    parser.add_argument("--metadata-file", required=True,
        help = "File path to metadata Excel file (.xlsx)")
    parser.add_argument("--output-dir", required=True,
        help = "Output directory")
    parser.add_argument("--batch-date", required=True,
        help = "Date in YYYYMMDD format")

    args = parser.parse_args()

    try:
        datetime.strptime(args.batch_date, '%Y%m%d')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMMDD")

    temp_dir = tempfile.mkdtemp(prefix='submissions_')
    LOG.debug(f"Created temp directory: {temp_dir}")

    batch_dir_name = 'Batch-' + args.batch_date
    fastq_dir_name = args.batch_date + '_fastq'

    # check that output directories don't already exist
    output_batch_dir = Path(args.output_dir, batch_dir_name)
    output_fastq_dir = Path(args.output_dir, fastq_dir_name)
    if os.path.isdir(output_batch_dir):
        raise FileExistsError(f"Output batch folder already exists: {output_batch_dir}. Aborting.")
    elif os.path.isdir(output_fastq_dir):
        raise FileExistsError(f"Output fastq folder already exists: {output_fastq_dir}. Aborting.")


    # create subdirectories
    batch_dir = Path(temp_dir, batch_dir_name)
    fastq_dir = Path(temp_dir, fastq_dir_name)
    LOG.debug(f"Creating subdirectories: {os.path.basename(batch_dir)}, {os.path.basename(fastq_dir)}")
    os.mkdir(batch_dir)
    os.mkdir(fastq_dir)

    # extract assembly results to temp fastq directory
    LOG.debug(f"Extracting compressed results file (.tar.gz) to: {fastq_dir}")
    assemby_results_file = tarfile.open(args.results_file)
    assemby_results_file.extractall(fastq_dir)
    assemby_results_file.close()

    LOG.debug("Copying/renaming metadata file to: external-metadata.xlsx")
    shutil.copy(args.metadata_file, temp_dir)
    original_metadata_filename = Path(args.metadata_file).name
    metadata_file_path = Path(temp_dir, original_metadata_filename)
    new_metadata_file_path = Path(batch_dir, 'external-metadata.xlsx')
    os.rename(metadata_file_path, new_metadata_file_path)

    LOG.debug(f"Copying FASTA (.fa) and metrics (.tsv) file to: {batch_dir}")
    assembly_filename_stem = Path(args.results_file).stem
    fasta_file = Path(fastq_dir, assembly_filename_stem).with_suffix('.fa')
    metrics_file = Path(fastq_dir, assembly_filename_stem).with_suffix('.metrics.tsv')
    shutil.copy(fasta_file, batch_dir)
    shutil.copy(metrics_file, batch_dir)

    # copy contents of temp_dir to output_dir
    shutil.copytree(batch_dir, output_batch_dir)
    shutil.copytree(fastq_dir, output_fastq_dir)

    LOG.debug(f"Output folders created and populated: {output_batch_dir}, {output_fastq_dir}")
    LOG.debug(f"Deleting temp directory {temp_dir}")
    shutil.rmtree(temp_dir)

    process_with_nextclade = None
    while True:
        process_with_nextclade = input("Process with NextClade? (y/n)").lower()
        if process_with_nextclade not in ['y','n']:
            print("Not a valid response")
            continue
        else:
            break

    if process_with_nextclade == 'y':
        LOG.debug("Pulling data from NextClade")

        # Get data files from NextClade
        result = Conda.run_command('run', 'nextclade', 'dataset', 'get',
            "--name=sars-cov-2",
            f"--output-dir={output_batch_dir}/data/sars-cov-2"
        )

        if result and len(result)==3 and result[2] == 0:
            LOG.debug("Successfully downloaded data/sars-cov-2 from NextClade.")
        else:
            raise Exception('Error: Nexclade sars-cov-2 data download failed.')


        LOG.debug("Analyzing FASTA file with NextClade.")
        fasta_file = Path(output_batch_dir, assembly_filename_stem).with_suffix('.fa')
        nextclade_output = Path(output_batch_dir,'nextclade.tsv')
        result = Conda.run_command('run', 'nextclade', 'run',
            f"--input-dataset={output_batch_dir}/data/sars-cov-2",
            f"--output-tsv={nextclade_output}",
            f"{fasta_file}"
        )
        if result[2] == 0:
           LOG.debug(f"NextClade processing complete: {nextclade_output}")
        else:
           raise Exception(f"Error: NextClade processing of {fasta_file} failed:\n {result}")

        #load nextclade.txv to dataframe
        nextclade_df = pd.read_csv(nextclade_output, sep="\t")

        nextclade_df["LIMS"] = nextclade_df["seqName"].str.extract("^[^\d]*(\d+)").astype(int)

        nextclade_df["aaSubstitutions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
        nextclade_df["aaInsertions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
        nextclade_df["aaDeletions_s_count"] = nextclade_df["aaSubstitutions"].str.split(",").apply(count_s_occurrences)
        nextclade_df["aaChanges_s_total"] = nextclade_df["aaSubstitutions_s_count"] + \
                                            nextclade_df["aaInsertions_s_count"] + \
                                            nextclade_df["aaDeletions_s_count"]

        excluded_vocs = nextclade_df[
                (nextclade_df["errors"].notnull()) |
                (nextclade_df["clade"].isin(["19A","19B","recombinant"])) |
                ("S" in nextclade_df["failedGenes"].str.split(",")) |
                (nextclade_df["aaChanges_s_total"] < 3) |
                (nextclade_df["qc.frameShifts.frameShifts"].str.split(",").apply(count_s_occurrences) > 0)]

        metadata = pd.read_excel(args.metadata_file)

        # remove twist positives from excluded VOCs
        twist_positives = metadata[metadata["Project"].str.lower()=="twist positive"]
        excluded_vocs=pd.merge(excluded_vocs, twist_positives, on=["LIMS"], how="outer", indicator=True)
        excluded_vocs=excluded_vocs[excluded_vocs["_merge"]=="left_only"]

        LOG.debug(f'Excluded VOCs:{excluded_vocs[["LIMS", "seqName", "clade", "qc.missingData.status", "aaChanges_s_total", "qc.frameShifts.frameShifts", "failedGenes", "errors"]]}')

        excluded_vocs[['LIMS']].to_csv(Path(output_batch_dir, 'excluded-vocs.txt'), header=False, index=False)
        LOG.debug(f"Saved {len(excluded_vocs)} NWGC ids to {Path(output_batch_dir, 'excluded-vocs.txt')}.")


    # process with VADR
    process_with_vadr = None
    docker_client = None
    while True:
        process_with_vadr = input("Process with VADR? (y/n)").lower()
        if process_with_vadr not in ['y','n']:
            print("Not a valid response")
            continue
        elif process_with_vadr == 'y':
            try:
                docker_client = docker.from_env()
                break
            except:
                print("Docker client could not be initiated. Make sure Docker is running.")
                continue
        else:
            break

    if process_with_vadr=='y' and docker_client:
        # pull latest docker image
        docker_client.images.pull('staphb/vadr')

        # run vadr
        docker_container = docker_client.containers.run('staphb/vadr',
            command = f'/bin/bash -c \
                "/opt/vadr/vadr/miniscripts/fasta-trim-terminal-ambigs.pl \
                    --minlen 50 --maxlen 30000 \
                    /data/{os.path.basename(fasta_file)} > trimmed-genbank.fasta; \
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
            volumes= {output_batch_dir: {'bind': '/data', 'mode': 'rw'}}
        )
        docker_output = docker_container.attach(stdout=True, stream=True, logs=True);
        for line in docker_output:
            print(line.decode('utf-8'))


    LOG.debug("Completed.")
