
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
from datetime import datetime

from export_lims_metadata import get_lims_sequencing_metadata, add_lims_metadata
from extract_sfs_identifiers import read_all_identifiers, find_sfs_identifiers

from aggregate_metrics import combine_metrics, append_metrics_counts

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
    parser.add_argument("--metadata-file", required=True,
        help = "File path to metadata Excel file (.xlsx)")
    parser.add_argument("--pathogen", type=str, required=True,
        choices=['sars-cov-2', 'rsv-a', 'rsv-b', 'flu-a', 'flu-b'])
    parser.add_argument("--subtype", type=str, required=False,
        choices=['h1n1', 'h3n2'])
    args = parser.parse_args()

    try:
        batch_date = datetime.strptime(args.batch_date, '%Y%m%d')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMMDD")

    batch_dir_name = 'Batch-' + args.batch_date

    # check that output directory doesn't already exist
    output_batch_dir = Path(args.output_dir, batch_dir_name)
    if os.path.isdir(output_batch_dir):
        if yes_no_cancel(f"Output batch folder already exists: {output_batch_dir}. Do you really want to proceed?"):
            output_batch_dir.mkdir(parents=True, exist_ok=True)
        else:
            raise FileExistsError(f"Output batch folder already exists: {output_batch_dir}. Aborting.")
    else:
        output_batch_dir.mkdir(parents=True, exist_ok=True)

    # shared outputs for all pathogens
    OUTPUT_PATHS = {
        'fasta': Path(output_batch_dir, args.pathogen).with_suffix('.fa'),
        'lims-metadata': Path(output_batch_dir, 'lims-metadata-with-county.csv'),
        'id3c-metadata': Path(output_batch_dir, 'id3c-metadata-with-county.csv'),
        'id3c-sample-barcodes': Path(output_batch_dir, 'id3c-sample-barcodes.csv'),
        'lims-sample-barcodes': Path(output_batch_dir, 'lims-sample-barcodes.csv'),
        'metadata': Path(output_batch_dir, 'external-metadata.xlsx'),

        # pathogen-specific outputs
        'metrics': Path(output_batch_dir, f'{args.pathogen}-metrics.tsv'),
        'nextclade': Path(output_batch_dir,f'nextclade-{args.pathogen}.tsv'),
        'nextclade-data': Path(output_batch_dir, f'data/{args.pathogen}'),
        'vadr-dir': Path(output_batch_dir, f'genbank-{args.pathogen}'), #this one is a folder
        'fasta-short-ids': Path(output_batch_dir, f'{args.pathogen}-short-ids.fa') # for now, create this file for every pathogen; however, it is only needed if GenBank trimming is performed, which should only be done for RSV.
    }

    # pathogen-specific outputs
    if args.pathogen == 'sars-cov-2':
        OUTPUT_PATHS.update({
            'previous-submissions': Path(output_batch_dir, 'sars-cov-2-previous-submissions.tsv'),
            'previous-submissions-other': Path(output_batch_dir, 'rsv-flu-previous-submissions.tsv'),
            'excluded-vocs': Path(output_batch_dir, 'excluded-vocs.txt'),
        })
    else:
        OUTPUT_PATHS.update({
            'previous-submissions': Path(output_batch_dir, 'rsv-flu-previous-submissions.tsv'),
            'previous-submissions-other': Path(output_batch_dir, 'sars-cov-2-previous-submissions.tsv'),
            'fasta-short-ids': Path(output_batch_dir, f'{args.pathogen}-short-ids.fa'),
            'genbank-trimmed-fasta': Path(output_batch_dir, f'{args.pathogen}-genbank-trimmed.fa'),
            'genbank-trim-log': Path(output_batch_dir, f'{args.pathogen}-genbank-trim.xml'),
        })

    # for some reason this value is sometimes a tuple where the second value is blank
    if isinstance(OUTPUT_PATHS['previous-submissions'], tuple):
        OUTPUT_PATHS['previous-submissions'] = OUTPUT_PATHS['previous-submissions'][0]

    standardize_and_qc_external_metadata(args.metadata_file, batch_date, OUTPUT_PATHS['metadata'])

    try:
        LOG.debug(f"Parsing FASTA file: {args.fasta}")
        sequences = list(SeqIO.parse(args.fasta, 'fasta'))

        # save FASTA file to output location
        shutil.copyfile(args.fasta, OUTPUT_PATHS['fasta'])
        LOG.debug(f"FASTA file saved to: {OUTPUT_PATHS['fasta']}")
    except Exception as err:
        LOG.error(f"There was a problem parsing sequences from {args.fasta}: \n{err}")


    if yes_no_cancel("Combine metrics files?"):
        LOG.debug("Combining metrics files")

        metrics_df, original_record_ids = append_metrics_counts(combine_metrics(args.metrics, args.pathogen), args.fasta, OUTPUT_PATHS['fasta-short-ids'], args.pathogen)
        metrics_df.to_csv(OUTPUT_PATHS['metrics'], sep='\t', index=False)

        LOG.debug(f"Combined metrics file saved to: {OUTPUT_PATHS['metrics']}")


    if yes_no_cancel("Process with NextClade?"):
        LOG.debug("Pulling data from NextClade")

        process_with_nextclade(OUTPUT_PATHS['nextclade-data'],
                            OUTPUT_PATHS['nextclade'],
                            OUTPUT_PATHS['fasta'],
                            args.pathogen,
                            args.subtype or None)

        LOG.debug("NextClade processing complete")

    if os.path.isfile(OUTPUT_PATHS['nextclade']):
        nextclade_df = pd.read_csv(OUTPUT_PATHS['nextclade'], sep='\t')
    else:
        raise FileNotFoundError(f"{OUTPUT_PATHS['nextclade']}")


    # create CSV of SFS sample barcodes
    identifiers = read_all_identifiers(OUTPUT_PATHS['metadata'])
    id3c_identifiers = find_sfs_identifiers(identifiers, '^((?!cascadia).)*$')
    lims_identifiers = find_sfs_identifiers(identifiers, '^(cascadia)$')

    sfs_missing_date = pd.DataFrame()

    LOG.debug(f"LIMS identifiers: {lims_identifiers}")

    if lims_identifiers is not None and not lims_identifiers.empty:
        lims_identifiers.to_csv(OUTPUT_PATHS['lims-sample-barcodes'], index=False)
        LOG.debug(f"Saved {len(lims_identifiers)} SFS identifiers to {OUTPUT_PATHS['lims-sample-barcodes']}.")

        if yes_no_cancel("Pull metadata from LIMS?"):
            LOG.debug("Pulling metadata from LIMS")
            OUTPUT_PATHS['lims-metadata'] = Path(output_batch_dir, 'lims-metadata-with-county.csv')

            lims_metadata = add_lims_metadata(lims_identifiers)
            lims_metadata.to_csv(OUTPUT_PATHS['lims-metadata'], index=False)
            LOG.debug(f"Successfully saved LIMS metadata to {OUTPUT_PATHS['lims-metadata']}")

            sfs_missing_date = lims_metadata[(lims_metadata['source'].str.lower()=='sfs')&(lims_metadata['collection_date'].isna())]

    LOG.debug(f"ID3C identifiers: {id3c_identifiers}")

    if id3c_identifiers is not None and not id3c_identifiers.empty:
        id3c_identifiers.to_csv(OUTPUT_PATHS['id3c-sample-barcodes'], index=False)

        if yes_no_cancel("Pull metadata from ID3C?"):
            LOG.debug("Pulling metadata from ID3C")
            OUTPUT_PATHS['id3c-metadata'] = Path(output_batch_dir, 'id3c-metadata-with-county.csv')

            # The `export_id3c_metadata` bash script is in the same location as current script
            export_id3c_metadata_script = Path(Path(__file__).resolve().parent, "export_id3c_metadata")

            # Temporarily point stdout to a file while running this script, then point it back
            stdout = sys.stdout
            with open(OUTPUT_PATHS['id3c-metadata'], 'w') as sys.stdout:
                result = Conda.run_command('run', f"{export_id3c_metadata_script}",
                    "--ignore-origin",
                    f"{OUTPUT_PATHS['id3c-sample-barcodes']}",
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
            sfs_missing_date.append(id3c_metadata_df[(id3c_metadata_df['source']=='SFS')&(id3c_metadata_df['collection_date'].isna())], ignore_index=True)


        if not sfs_missing_date.empty:
            LOG.debug(f"Warning: {len(sfs_missing_date)} SFS samples missing collection date:\n {sfs_missing_date}")
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
        pull_previous_submissions(OUTPUT_PATHS['previous-submissions'], OUTPUT_PATHS['previous-submissions-other'], args.pathogen)
        LOG.debug("Finished pulling previous submissions")

    if args.pathogen == 'sars-cov-2':
        calculate_excluded_vocs(nextclade_df, OUTPUT_PATHS['excluded-vocs'], identifier_format='sfs_sample_barcode')

    if not args.pathogen.startswith('flu'): # don't trim or run vadr for flu
        if yes_no_cancel("Trim with GenBank script?"):
            # The genbank trim script has a limit of 50 characters on identifiers; save a fasta file with just NWGC ids
            #fasta_out = SeqIO.FastaIO.FastaWriter(OUTPUT_PATHS['fasta-short-ids'], wrap=None)
            #fasta_out.write_file(sequences)

            # Trim with GenBank's trimming script
            result = Conda.run_command('run',
                f'{os.path.realpath(os.path.dirname(__file__))}/fastaedit_public',
                "-in", f"{OUTPUT_PATHS['fasta-short-ids']}",
                "-trim_ambig_bases",
                "-out_seq_file", f"{OUTPUT_PATHS['genbank-trimmed-fasta']}",
                "-out", f"{OUTPUT_PATHS['genbank-trim-log']}"
            )

            # Replace the shortened IDs with the original long IDs. GenBank's trim script adds a `lcl|` prefix to each ID
            # which must be removed.
            trimmed_records = list(SeqIO.parse(OUTPUT_PATHS['genbank-trimmed-fasta'], 'fasta'))
            for r in trimmed_records:
                r.id = r.description = original_record_ids[r.id.lstrip('lcl|')]
            SeqIO.write(trimmed_records, OUTPUT_PATHS['genbank-trimmed-fasta'], 'fasta-2line')

            process_with_vadr(OUTPUT_PATHS['genbank-trimmed-fasta'], output_batch_dir, args.pathogen, interactive=True)
        else:
            process_with_vadr(OUTPUT_PATHS['fasta'], output_batch_dir, args.pathogen, interactive=True)

    print("Completed file prep.\n\n")

    previous_submissions = pd.read_csv(OUTPUT_PATHS['previous-submissions'], sep='\t')
    valid_gisaid_ids = previous_submissions[previous_submissions['strain_name'].str.contains('USA/[A-Z]{2}-S\d*/\d{4}', na=False)][['strain_name']]
    next_avail_strain_id = valid_gisaid_ids['strain_name'].apply( lambda x: x[x.find("-S")+2 : x.find("/",x.find("-S"))]).astype(int).max() + 1

    if all([x.exists() for x in OUTPUT_PATHS.values()]):
        # create_submissions scripts are in the same folder as current script
        if args.pathogen == 'sars-cov-2':
            create_submissions_script = Path(Path(__file__).resolve().parent, "create_submissions_sars_cov_2.py")
        elif args.pathogen.startswith('rsv-'):
            create_submissions_script = Path(Path(__file__).resolve().parent, "create_submissions_rsv.py")
        elif args.pathogen.startswith('flu-'):
            create_submissions_script = Path(Path(__file__).resolve().parent, "create_submissions_flu.py")


        excluded_vocs_arg = f"--excluded-vocs {OUTPUT_PATHS.get('excluded-vocs')} \\\n" if (OUTPUT_PATHS.get('excluded-vocs')) else ""
        vadr_dir_arg = f"--vadr-dir {OUTPUT_PATHS.get('vadr-dir')} \\\n" if (OUTPUT_PATHS.get('vadr-dir')) else ""
        prev_submissions_other_arg = f"--previous-submissions-rsv-flu {OUTPUT_PATHS.get('previous-submissions-other')} \\\n" if args.pathogen=='sars-cov-2' else f"--previous-submissions-sars-cov-2 {OUTPUT_PATHS.get('previous-submissions-other')} \\\n"
        subtype_arg = f"--subtype {args.pathogen.split('-')[-1]} \\\n" if (args.pathogen.startswith('rsv-') or args.pathogen.startswith('flu-')) else ""
        print("Review outputs, then run this command to create submissions:\n\n"
            f"python3 {create_submissions_script} \\\n"
            f"--batch-name {args.batch_date} \\\n"
            f"--id3c-metadata {OUTPUT_PATHS['id3c-metadata']} \\\n"
            f"--lims-metadata {OUTPUT_PATHS['lims-metadata']} \\\n"
            f"--fasta {OUTPUT_PATHS['fasta']} \\\n"
            f"--nextclade {OUTPUT_PATHS['nextclade']} \\\n"
            f"--previous-submissions {OUTPUT_PATHS['previous-submissions']} \\\n"
            f"{prev_submissions_other_arg}"
            f"--output-dir {Path(output_batch_dir, 'submissions', args.pathogen)} \\\n"
            f"--strain-id {'<not found>' if pd.isna(next_avail_strain_id) else next_avail_strain_id} \\\n"
            f"--metrics {OUTPUT_PATHS['metrics']} \\\n"
            f"--test-name MIPsSEQ \\\n"
            f"--metadata {OUTPUT_PATHS['metadata']} \\\n"
            f"{excluded_vocs_arg}"
            f"{vadr_dir_arg}"
            f"{subtype_arg}"
            f"--gisaid-username <your username> \n"
        )
    else:
        print("Warning: Missing required inputs to create submissions:")
        print('\n'.join([str(x) for x in OUTPUT_PATHS.values() if not x.exists()]))
