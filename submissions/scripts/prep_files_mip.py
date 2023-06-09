
import os
import re
import sys
import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd
from Bio import SeqIO
from utils import *
import shutil

from export_lims_metadata import get_lims_sequencing_metadata

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
    parser.add_argument('--fasta', required=True,
        help='The fasta file to process.')
    parser.add_argument('--metrics', required=True,
        help='The metrics file to process.', nargs='+')
    parser.add_argument("--output-dir", required=True,
        help = "Output directory")
    parser.add_argument("--batch-date", required=True,
        help = "Date in YYYYMMDD format")
    args = parser.parse_args()

    try:
        batch_date = datetime.strptime(args.batch_date, '%Y%m%d')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMMDD")

    batch_dir_name = 'Batch-' + args.batch_date

    # check that output directory doesn't already exist
    output_batch_dir = Path(args.output_dir, batch_dir_name)
    if os.path.isdir(output_batch_dir):
        if yes_no_cancel("Output batch folder already exists: {output_batch_dir}. Do you really want to proceed?"):
            output_batch_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise FileExistsError(f"Output batch folder already exists: {output_batch_dir}. Aborting.")
    else:
        output_batch_dir.mkdir(parents=True, exist_ok=True)

    OUTPUT_PATHS = {
        'nextclade': Path(output_batch_dir,'nextclade.tsv'),
        'nextclade-sars-cov-2': Path(output_batch_dir, 'data/sars-cov-2'),

        'fasta': Path(output_batch_dir, Path(args.fasta).stem).with_suffix('.fa'),
        'metrics': Path(output_batch_dir, 'metrics.tsv'),
        'lims-metadata': Path(output_batch_dir, 'lims-metadata-with-county.csv'),

        'previous-submissions': Path(output_batch_dir, 'previous-submissions.tsv'),
        'excluded-vocs': Path(output_batch_dir, 'excluded-vocs.txt'),
        'vadr-dir': Path(output_batch_dir, 'genbank'), #this one is a folder

        'sfs-retro-sample-barcodes': Path(output_batch_dir, 'sfs-retro-sample-barcodes.csv'),
    }

    try:
        LOG.debug(f"Parsing FASTA file: {args.fasta}")
        sequences = SeqIO.parse(args.fasta, 'fasta')

        # save FASTA file to output location
        shutil.copyfile(args.fasta, OUTPUT_PATHS['fasta'])
        LOG.debug(f"FASTA file saved to: {OUTPUT_PATHS['fasta']}")
    except Exception as err:
        LOG.error(f"There was a problem parsing sequences from {args.fasta}: \n{err}")


    LOG.debug(f"Aggregating {len(args.metrics)} metrics files.")
    PICARD_METRICS_COLS = [
        'GENOME_TERRITORY',
        'MEAN_COVERAGE',
        'SD_COVERAGE',
        'MEDIAN_COVERAGE',
        'MAD_COVERAGE',
        'PCT_EXC_ADAPTER',
        'PCT_EXC_MAPQ',
        'PCT_EXC_DUPE',
        'PCT_EXC_UNPAIRED',
        'PCT_EXC_BASEQ',
        'PCT_EXC_OVERLAP',
        'PCT_EXC_CAPPED',
        'PCT_EXC_TOTAL',
        'PCT_1X',
        'PCT_5X',
        'PCT_10X',
        'PCT_15X',
        'PCT_20X',
        'PCT_25X',
        'PCT_30X',
        'PCT_40X',
        'PCT_50X',
        'PCT_60X',
        'PCT_70X',
        'PCT_80X',
        'PCT_90X',
        'PCT_100X',
        'HET_SNP_SENSITIVITY',
        'HET_SNP_Q',
    ]

    metrics_df = pd.DataFrame()
    LOG.debug(f"Processing metrics files: {', '.join(args.metrics)}")

    for metrics_file in args.metrics:
        # parse barcode from filename
        identifier = Path(metrics_file).stem.split('.')[0]

        # metrics TSV for each sample should be on lines 7 (headers) and 8 (data) of individual
        # metrics files output by Picard
        df = pd.read_csv(metrics_file, sep='\t', skiprows=6, nrows=1)

        if not all(c in df.columns for c in PICARD_METRICS_COLS):
            raise Exception(f"Metrics file {metrics_file} does not match expected columns: {', '.join(PICARD_METRICS_COLS)}")
        elif df.empty:
            raise Exception(f"Metrics file {metrics_file} does not contain expected TSV data.")

        # insert sample id as first column and append to metrics dataframe
        df.insert(0, 'SampleId', identifier)
        metrics_df = metrics_df.append(df[PICARD_METRICS_COLS + ['SampleId']], ignore_index=True)


    # Calculate additional metrics from sequence and add to metrics dataframe
    for record in sequences:
        identifier = str(record.id.split('|')[0]).lower()
        metrics = sequence_metrics(record.seq, record.id)

        # update corresponding row with additional metrics
        for k,v in metrics.items():
            metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, k] = str(v)

    metrics_df.to_csv(OUTPUT_PATHS['metrics'], sep='\t', index=False)
    LOG.debug(f"Combined metrics file saved to: {OUTPUT_PATHS['metrics']}")


    if yes_no_cancel("Process with NextClade?"):
        LOG.debug("Pulling data from NextClade")

        process_with_nextclade(OUTPUT_PATHS['nextclade-sars-cov-2'],
                            OUTPUT_PATHS['nextclade'],
                            OUTPUT_PATHS['fasta'])

        LOG.debug("NextClade processing complete")

    if os.path.isfile(OUTPUT_PATHS['nextclade']):
        nextclade_df = pd.read_csv(OUTPUT_PATHS['nextclade'], sep='\t')
    else:
        raise FileNotFoundError(f"{OUTPUT_PATHS['nextclade']}")


    identifiers = [seqname.split('|')[0] for seqname in nextclade_df['seqName']]
    non_control_identifiers = [x for x in identifiers if not str(x).lower().endswith('-pos-con')]

    missing_date = pd.DataFrame()

    if yes_no_cancel("Pull metadata from LIMS?"):
        LOG.debug("Pulling metadata from LIMS")
        OUTPUT_PATHS['lims-metadata'] = Path(output_batch_dir, 'lims-metadata-with-county.csv')

        lims_metadata = get_lims_sequencing_metadata(non_control_identifiers)

        # For Cascadia, use sfs_identifier_for_doh_reporting value as sfs_sample_identifier
        lims_metadata.loc[lims_metadata['source'].str.lower() == 'cascadia','sfs_sample_identifier'] = lims_metadata['sfs_identifier_for_doh_reporting']

        # Only retros should be captured by this filter. Barcodes are then used to pull metadata from ID3C
        metadata_from_id3c = lims_metadata[lims_metadata['source'].str.lower().isin(['sch', 'sfs'])][['sfs_sample_barcode']]

        # Only non-retros should be captured by this filter. LIMS metadata saved to CSV
        lims_metadata = lims_metadata[~lims_metadata['source'].str.lower().isin(['sch', 'sfs'])]
        lims_metadata.to_csv(OUTPUT_PATHS['lims-metadata'], index=False)

        # These additional fields are required by the `export_id3c_metadata` script
        metadata_from_id3c['nwgc_id'] = ""
        metadata_from_id3c['sample_origin'] = ""
        metadata_from_id3c[['nwgc_id', 'sfs_sample_barcode', 'sample_origin']].to_csv(OUTPUT_PATHS['sfs-retro-sample-barcodes'], index=False)
    else:
        lims_metadata = pd.read_csv(OUTPUT_PATHS['lims-metadata'])
        metadata_from_id3c = pd.read_csv(OUTPUT_PATHS['sfs-retro-sample-barcodes'])

    missing_date = lims_metadata[lims_metadata['collection_date'].isna()]

    # Check records that need metadata from ID3C (i.e. retros)
    if (len(metadata_from_id3c) > 0):
        OUTPUT_PATHS['id3c-metadata'] = Path(output_batch_dir, 'id3c-metadata-with-county.csv')

        if yes_no_cancel("Pull metadata from ID3C?"):
            LOG.debug("Pulling metadata from ID3C")

            # The `export_id3c_metadata` bash script is in the same location as current script
            export_id3c_metadata_script = Path(Path(__file__).resolve().parent, "export_id3c_metadata")

            # Temporarily point stdout to a file while running this script, then point it back
            stdout = sys.stdout
            with open(OUTPUT_PATHS['id3c-metadata'], 'w') as sys.stdout:
                result = Conda.run_command('run', f"{export_id3c_metadata_script}",
                    "--ignore-origin",
                    f"{OUTPUT_PATHS['sfs-retro-sample-barcodes']}",
                    stdout="STDOUT",
                    stderr="STDERR"
                )
            sys.stdout = stdout

            if result and len(result)==3 and result[2] == 0 and OUTPUT_PATHS['id3c-metadata'].exists() and OUTPUT_PATHS['id3c-metadata'].stat().st_size > 0:
                LOG.debug(f"Successfully saved ID3C metadata to {OUTPUT_PATHS['id3c-metadata']}")
            else:
                raise Exception(f"Error pulling metadata from ID3C.\n {result}")

        # Check for SFS samples with missing collection date
        id3c_metadata_df = pd.read_csv(OUTPUT_PATHS['id3c-metadata'])
        missing_date.append(id3c_metadata_df[id3c_metadata_df['collection_date'].isna()], ignore_index=True)


    if not missing_date.empty:
        LOG.debug(f"Warning: {len(missing_date)} samples missing collection date:\n {missing_date}")
        while True:
            user_input = input(f"Continue preparing files? [y]es/[n]o: ").lower()
            if user_input not in ['y','n']:
                print("Not a valid response")
                continue
            elif user_input == 'n':
                sys.exit()
            else:
                break

    if yes_no_cancel("Pull latest previous submissions file from Github?"):
        LOG.debug("Pulling previous submissions")

        pull_previous_submissions(OUTPUT_PATHS['previous-submissions'])
        LOG.debug("Finished pulling previous submissions")

    # calculate excluded VOCs
    calculate_excluded_vocs(nextclade_df, OUTPUT_PATHS['excluded-vocs'], identifier_format='sfs_sample_barcode')

    process_with_vadr(OUTPUT_PATHS['fasta'], output_batch_dir, interactive=True)

    print("Completed file prep.\n\n")

    previous_submissions = pd.read_csv(OUTPUT_PATHS['previous-submissions'], sep='\t')
    valid_gisaid_ids = previous_submissions[previous_submissions['strain_name'].str.contains('USA/[A-Z]{2}-S\d*/\d{4}', na=False)][['strain_name']]
    next_avail_strain_id = valid_gisaid_ids['strain_name'].apply( lambda x: x[x.find("-S")+2 : x.find("/",x.find("-S"))]).astype(int).max() + 1

    if all([x.exists() for x in OUTPUT_PATHS.values()]):
        # create_submissions script is in the same folder as current script
        create_submissions_script = Path(Path(__file__).resolve().parent, "create_submissions.py")

        print("Review outputs, then run this command to create submissions:\n\n"
            f"python3 {create_submissions_script} \\\n"
            f"--batch-name {args.batch_date} \\\n"
            f"--id3c-metadata {OUTPUT_PATHS['id3c-metadata']} \\\n"
            f"--lims-metadata {OUTPUT_PATHS['lims-metadata']} \\\n"
            f"--fasta {OUTPUT_PATHS['fasta']} \\\n"
            f"--nextclade {OUTPUT_PATHS['nextclade']} \\\n"
            f"--previous-submissions {OUTPUT_PATHS['previous-submissions']} \\\n"
            f"--excluded-vocs {OUTPUT_PATHS['excluded-vocs']} \\\n"
            f"--vadr-dir {OUTPUT_PATHS['vadr-dir']} \\\n"
            f"--output-dir {output_batch_dir} \\\n"
            f"--strain-id {next_avail_strain_id} \\\n"
            f"--metrics {OUTPUT_PATHS['metrics']} \\\n"
            f"--test-name MIPsSEQ \\\n"
            f"--gisaid-username <your username> \n"
        )
    else:
        print("Warning: Missing required inputs to create submissions:")
        print('\n'.join([str(x) for x in OUTPUT_PATHS.values() if not x.exists()]))
