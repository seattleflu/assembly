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
import conda.cli.python_api as Conda


LOG_LEVEL = os.environ.get("LOG_LEVEL", "debug").upper()

logging.basicConfig(
    level = logging.ERROR,
    format = "[%(asctime)s] %(levelname)-8s %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S%z")

logging.captureWarnings(True)

LOG = logging.getLogger(__name__)
LOG.setLevel(LOG_LEVEL)

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

    LOG.debug("Creating empty file: excluded-vocs.txt (to be populated manually)")
    open(Path(batch_dir, 'excluded-vocs.txt'), 'w').close()

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

    LOG.debug("Completed.")
