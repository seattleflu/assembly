"""
Prepares files and directories for sequencing submissions.

"""
import os
import re
import sys
import argparse
import tarfile
from pathlib import Path
from datetime import datetime
import logging
import pandas as pd
import conda.cli.python_api as Conda
import docker
from openpyxl import load_workbook
from extract_sfs_identifiers import read_all_identifiers, find_sfs_identifiers
import requests
import getpass

from export_lims_metadata import add_lims_metadata

LOG_LEVEL = os.environ.get("LOG_LEVEL", "debug").upper()

logging.basicConfig(
    level = logging.ERROR,
    format = "[%(asctime)s] %(levelname)-8s %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S%z")

logging.captureWarnings(True)

LOG = logging.getLogger(__name__)
LOG.setLevel(LOG_LEVEL)

def yes_no_cancel(message):
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

def count_s_occurrences(values):
    count = 0
    if not isinstance(values, list):
        return count
    for item in values:
        if isinstance(item, str) and item.lower().strip().startswith("s:"):
            count += 1
    return count

def standardize_date(date_value):
    if pd.isna(date_value) or (isinstance(date_value, str) and date_value.lower().strip() == 'na'):
        return None
    else:
        # So far, two date formats have been encountered in metadata Excel files, resulting in different formats in
        # resulting dataframes. If `date_value` does not match either of these formats, a ValueError exception will be
        # thrown. If the rejected date is in a usable format, it should be added here.
        datetime_val = pd.to_datetime(date_value, format='%Y-%m-%d %H:%M:%S', errors="coerce")
        if pd.isnull(datetime_val):
            datetime_val = pd.to_datetime(date_value, format='%m/%d/%Y')

        # The minus sign before m and d removes leading zeros, to match the format of the original Excel file.
        # Note: this may only work on unix, windows may need to use `#` instead of `-`
        return datetime_val.strftime('%-m/%-d/%Y')


def standardize_barcode(barcode):
    if pd.isna(barcode):
        return None
    elif isinstance(barcode, str):
        standardized_barcode = barcode.strip().lower()
    elif isinstance(barcode, int):
        standardized_barcode = str(barcode)
    else:
        raise TypeError(f"Barcodes must be type str or int. Invalid type {type(barcode)} on barcode: {barcode}")

    # In addition to standard format, allow barcodes ending in `_exp` which typically don't need to be submitted
    # but should be confirmed with sequencing team before dropping.
    if (bool(re.match('^[0-9a-f]{8}$', standardized_barcode)) or
        bool(re.match('^[0-9a-f]{8}_exp$', standardized_barcode)) or
        standardized_barcode in ["twist positive", "blank"]):
        return standardized_barcode
    else:
        raise ValueError(f"Identifiers must be hexidecimal barcodes, `blank`, or `twist positive`. Invalid barcode: {standardized_barcode}")

def validate_collection_dates(collection_dates: pd.Series, batch_date: datetime):
    cutoff_datetime = datetime(2020, 2, 1)

    out_of_range = collection_dates.dropna()[(pd.to_datetime(collection_dates, format='%m/%d/%Y') < cutoff_datetime)|(pd.to_datetime(collection_dates, format='%m/%d/%Y') > batch_date)]

    if not out_of_range.empty:
        raise Exception(f"Error: Collection date(s) out of range ({cutoff_datetime.strftime('%Y-%m-%d')} to {batch_date.strftime('%Y-%m-%d')}):\n {out_of_range}")

    if len(collection_dates.dropna().unique()) < 2:
        raise Exception(f"Error: All collection dates values are the same:\n {collection_dates.unique()}")



def standardize_and_qc_external_metadata(metadata_filename:str, batch_date:datetime):
    # Standardize and validate barcode and date values on external metadata Excel file. Loads 2 sheets to dataframes
    # and writes them back to Excel if data passes all checks.

    # Metadata sheet
    external_metadata_metadata = pd.read_excel(metadata_filename, sheet_name='Metadata', converters={'LIMS':int, 'lab_accession_id':str})
    external_metadata_metadata['lab_accession_id'] = external_metadata_metadata['lab_accession_id'].str.strip()
    # only standardize internal barcodes
    external_metadata_metadata['lab_accession_id'] = external_metadata_metadata.apply(lambda row: standardize_barcode(row['lab_accession_id']) if str(row['Project']).startswith(('starita_bbi_', 'SeattleChildrensDirect_')) else row['lab_accession_id'], axis=1)
    external_metadata_metadata['collection_date'] = external_metadata_metadata['collection_date'].apply(standardize_date)
    external_metadata_metadata['Seq date'] = external_metadata_metadata['Seq date'].apply(standardize_date)

    # Set submitting_lab to sentinel for records with Project = sentinel
    external_metadata_metadata.loc[external_metadata_metadata['Project'].str.lower()=='sentinel', 'submitting_lab'] = 'sentinel'

    # Samplify FC Data sheet
    external_metadata_samplify_fc_data = pd.read_excel(args.metadata_file, sheet_name='Samplify FC Data', header=1, converters={'Sample ID':int, 'Investigator\'s sample ID':str})
    external_metadata_samplify_fc_data['Investigator\'s sample ID'] = external_metadata_samplify_fc_data['Investigator\'s sample ID'].str.strip()
    # only standardize internal barcodes
    external_metadata_samplify_fc_data['Investigator\'s sample ID'] = external_metadata_samplify_fc_data.apply(lambda row: standardize_barcode(row['Investigator\'s sample ID']) if str(row['Project']).startswith(('starita_bbi_', 'SeattleChildrensDirect_')) else row['Investigator\'s sample ID'], axis=1)
    external_metadata_samplify_fc_data['Collection date'] = external_metadata_samplify_fc_data['Collection date'].apply(standardize_date)

    # Only the collection date from the Metadata sheet is used for creating submission files.
    # Validating to confirm no collection dates prior to Feb 2020 and not all dates are identical.
    validate_collection_dates(external_metadata_metadata['collection_date'], batch_date)

    sentinel_ids = ['blank', 'twist positive', 'pbs', 'lp_blank']
    # all records with unique identifiers should be present in both sheets
    non_matching_ids = pd.merge(external_metadata_metadata[~external_metadata_metadata["lab_accession_id"].str.lower().isin(sentinel_ids)],
        external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Investigator\'s sample ID'].notna()],
        how="outer",
        left_on=["LIMS", "lab_accession_id"],
        right_on=["Sample ID", "Investigator\'s sample ID"],
        indicator=True).query("_merge != 'both'")[['LIMS', 'lab_accession_id', 'Sample ID',  'Investigator\'s sample ID']]

    # for blanks and twist positives, the sample IDs should present on both sheets, but with empty `Investigator's sample ID` on Samplify FC Data sheet.
    non_matching_blank_twist_pos = pd.merge(external_metadata_metadata[external_metadata_metadata["lab_accession_id"].str.lower().isin(sentinel_ids)],
        external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Investigator\'s sample ID'].isna()],
        how="outer",
        left_on=["LIMS"],
        right_on=["Sample ID"],
        indicator=True).query("_merge != 'both'")[['LIMS', 'lab_accession_id', 'Sample ID',  'Investigator\'s sample ID']]

    assert len(non_matching_ids) == 0, f"Error: Non-matching IDs found in metadata file:\n {non_matching_ids}"
    assert len(non_matching_blank_twist_pos) == 0, f"Error: Non-matching blanks and/or twist positives found in metadata file:\n{non_matching_blank_twist_pos}"

    # Identify and optionally remove expirimental samples
    exp_samples = external_metadata_metadata[external_metadata_metadata['lab_accession_id'].str.endswith('_exp', na=False)]
    exp_samples_samplify = external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Investigator\'s sample ID'].str.endswith('_exp', na=False)]
    if not exp_samples.empty:
        print(f"Expirimental samples found: {len(exp_samples)}\n {exp_samples.to_string()}")

        if yes_no_cancel("Drop expirimental samples?"):
            external_metadata_metadata.drop(exp_samples.index, inplace=True)
            external_metadata_samplify_fc_data.drop(exp_samples_samplify.index, inplace=True)
        else:
            external_metadata_metadata = exp_samples
            external_metadata_samplify_fc_data = exp_samples_samplify

    # Identify and optionally remove cascadia samples
    cascadia_samples = external_metadata_metadata[external_metadata_metadata['county'].str.lower().str.strip() == 'cascadia']
    cascadia_samples_samplify = external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Origin'].str.lower().str.strip() == 'cascadia']
    if not exp_samples.empty or not cascadia_samples_samplify.empty:
        print(f"Cascadia samples found:\n {exp_samples.to_string()} \n {cascadia_samples_samplify.to_string()}")

        if yes_no_cancel("Drop Cascadia samples?"):
            external_metadata_metadata.drop(cascadia_samples.index, inplace=True)
            external_metadata_samplify_fc_data.drop(cascadia_samples_samplify.index, inplace=True)
        else:
            # Cascadia samples should be marked as sequence_reson: Other
            external_metadata_metadata.loc[cascadia_samples.index, 'sequence_reason'] = 'Other'

    # The Samplify FC Data sheet has an extra row before the headers. Stashing this value to reinsert the row when writing back to Excel.
    metadata_wb = load_workbook(args.metadata_file)
    samplify_fc_data_a1 = metadata_wb['Samplify FC Data']['A1'].value

    with pd.ExcelWriter(OUTPUT_PATHS['metadata'], engine='xlsxwriter', date_format='m/d/yyyy') as writer:
        external_metadata_metadata.to_excel(writer, sheet_name='Metadata', index=False)
        # put column headers on 2nd row of Samplify FC Data sheet, then populate first row with `samplify_fc_data_a1`` value
        external_metadata_samplify_fc_data.to_excel(writer, sheet_name='Samplify FC Data', startrow = 1, index=False)
        worksheet = writer.sheets['Samplify FC Data']
        worksheet.write(0, 0, samplify_fc_data_a1)
        # set format for Crt values to 2 decimals so whole numbers aren't displayed as integers (matching source file format)
        format1 = writer.book.add_format({'num_format': '0.00'})
        worksheet.set_column('O:O', None, format1)
        writer.save()

    LOG.debug(f"Saved external metadata file to {OUTPUT_PATHS['metadata']}")


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
        batch_date = datetime.strptime(args.batch_date, '%Y%m%d')
    except ValueError:
        raise ValueError("Incorrect date format, should be YYYYMMDD")

    batch_dir_name = 'Batch-' + args.batch_date

    # check that output directory doesn't already exist
    output_batch_dir = Path(args.output_dir, batch_dir_name)
    if os.path.isdir(output_batch_dir):
        raise FileExistsError(f"Output batch folder already exists: {output_batch_dir}. Aborting.")
    else:
        output_batch_dir.mkdir(parents=True, exist_ok=True)

    OUTPUT_PATHS = {
        'metadata': Path(output_batch_dir, 'external-metadata.xlsx'),
        'metrics': Path(output_batch_dir, Path(args.results_file).stem).with_suffix('.metrics.tsv'),
        'nextclade': Path(output_batch_dir,'nextclade.tsv'),
        'previous-submissions': Path(output_batch_dir, 'previous-submissions.tsv'),
        'fasta': Path(output_batch_dir, Path(args.results_file).stem).with_suffix('.fa'),
        'excluded-vocs': Path(output_batch_dir, 'excluded-vocs.txt'),
        'vadr-dir': Path(output_batch_dir, 'genbank'), #this one is a folder

        'sfs-non-retro-sample-barcodes': Path(output_batch_dir, 'sfs-non-retro-sample-barcodes.csv'),
        'sfs-retro-sample-barcodes': Path(output_batch_dir, 'sfs-retro-sample-barcodes.csv'),
        'nextclade-sars-cov-2': Path(output_batch_dir, 'data/sars-cov-2'),
    }

    # extract assembly results to temp fastq directory
    LOG.debug(f"Extracting compressed results file (.tar.gz) to: {output_batch_dir}")
    assemby_results_file = tarfile.open(args.results_file)
    assemby_results_file.extractall(output_batch_dir)
    assemby_results_file.close()

    standardize_and_qc_external_metadata(args.metadata_file, batch_date)

    # create CSV of SFS sample barcodes
    identifiers = read_all_identifiers(OUTPUT_PATHS['metadata'])
    sfs_retro_identifiers = find_sfs_identifiers(identifiers, '^.*_retro$')
    sfs_non_retro_identifiers = find_sfs_identifiers(identifiers, '^(?!.*_retro).*$')

    sfs_missing_date = pd.DataFrame()

    LOG.debug(f"sfs_non_retro_identifiers: {sfs_non_retro_identifiers}")

    if sfs_non_retro_identifiers is not None:
        sfs_non_retro_identifiers.to_csv(OUTPUT_PATHS['sfs-non-retro-sample-barcodes'], index=False)
        LOG.debug(f"Saved {len(sfs_non_retro_identifiers)} SFS identifiers to {OUTPUT_PATHS['sfs-non-retro-sample-barcodes']}.")

        #sfs_missing_date = None
        if yes_no_cancel("Pull metadata from LIMS?"):
            LOG.debug("Pulling metadata from LIMS")
            OUTPUT_PATHS['lims-metadata'] = Path(output_batch_dir, 'lims-metadata-with-county.csv')

            lims_metadata = add_lims_metadata(sfs_non_retro_identifiers)
            lims_metadata.to_csv(OUTPUT_PATHS['lims-metadata'], index=False)
            LOG.debug(f"Successfully saved LIMS metadata to {OUTPUT_PATHS['lims-metadata']}")

            sfs_missing_date = lims_metadata[(lims_metadata['source'].str.lower()=='sfs')&(lims_metadata['collection_date'].isna())]

    if sfs_retro_identifiers is not None:
        sfs_retro_identifiers.to_csv(OUTPUT_PATHS['sfs-retro-sample-barcodes'], index=False)

        if yes_no_cancel("Pull metadata from ID3C?"):
            LOG.debug("Pulling metadata from ID3C")
            OUTPUT_PATHS['id3c-metadata'] = Path(output_batch_dir, 'id3c-metadata-with-county.csv')

            # The `export_id3c_metadata` bash script is in the same location as current script
            export_id3c_metadata_script = Path(Path(__file__).resolve().parent, "export_id3c_metadata")

            # Temporarily point stdout to a file while running this script, then point it back
            stdout = sys.stdout
            with open(OUTPUT_PATHS['id3c-metadata'], 'w') as sys.stdout:
                result = Conda.run_command('run', f"{export_id3c_metadata_script}",
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
            sfs_missing_date.append(id3c_metadata_df[(id3c_metadata_df['source']=='SFS')&(id3c_metadata_df['collection_date'].isna())], ignore_index=True)


        if sfs_missing_date is not None and not sfs_missing_date.empty:
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

    while True:
        if yes_no_cancel("Pull previous submissions file from Github?"):
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
                    with open(OUTPUT_PATHS['previous-submissions'], 'w') as f:
                        f.write(r.text)
                    break
                else:
                    LOG.warning(f"Error downloading previous submissions file: {r.status_code}")
                    continue

            except requests.exceptions.HTTPError as e:
                LOG.warning(f"Error: {e.response.text}")
                continue
        else:
            break


    if yes_no_cancel("Process with NextClade?"):
        LOG.debug("Pulling data from NextClade")

        # Get data files from NextClade
        result = Conda.run_command('run', 'nextclade', 'dataset', 'get',
            "--name=sars-cov-2",
            f"--output-dir={OUTPUT_PATHS['nextclade-sars-cov-2']}"
        )

        if result and len(result)==3 and result[2] == 0:
            LOG.debug("Successfully downloaded data/sars-cov-2 from NextClade.")
        else:
            raise Exception('Error: Nexclade sars-cov-2 data download failed.')


        LOG.debug("Analyzing FASTA file with NextClade.")
        result = Conda.run_command('run', 'nextclade', 'run',
            f"--input-dataset={OUTPUT_PATHS['nextclade-sars-cov-2']}",
            f"--output-tsv={OUTPUT_PATHS['nextclade']}",
            f"{OUTPUT_PATHS['fasta']}"
        )
        if result[2] == 0:
           LOG.debug(f"NextClade processing complete: {OUTPUT_PATHS['nextclade']}")
        else:
           raise Exception(f"Error: NextClade processing of {OUTPUT_PATHS['fasta']} failed:\n {result}")

        #load nextclade.txv to dataframe
        nextclade_df = pd.read_csv(OUTPUT_PATHS['nextclade'], sep="\t")

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
                ("S" in nextclade_df["failedGenes"].astype(str).str.split(",")) |
                (nextclade_df["aaChanges_s_total"] < 3) |
                (nextclade_df["qc.frameShifts.frameShifts"].str.split(",").apply(count_s_occurrences) > 0)]

        metadata = pd.read_excel(args.metadata_file)

        # remove twist positives from excluded VOCs
        twist_positives = metadata[metadata["lab_accession_id"].str.lower()=="twist positive"]
        excluded_vocs=pd.merge(excluded_vocs, twist_positives, on=["LIMS"], how="outer", indicator=True)
        excluded_vocs=excluded_vocs[excluded_vocs["_merge"]=="left_only"]

        LOG.debug(f'Excluded VOCs:{excluded_vocs[["LIMS", "seqName", "clade", "qc.missingData.status", "aaChanges_s_total", "qc.frameShifts.frameShifts", "failedGenes", "errors"]]}')

        excluded_vocs[['LIMS']].to_csv(OUTPUT_PATHS['excluded-vocs'], header=False, index=False)
        LOG.debug(f"Saved {len(excluded_vocs)} NWGC ids to {OUTPUT_PATHS['excluded-vocs']}.")


    # process with VADR
    process_with_vadr = None
    docker_client = None
    while True:
        process_with_vadr = yes_no_cancel("Process with VADR? (y/n)")
        if process_with_vadr:
            try:
                docker_client = docker.from_env()
                break
            except:
                print("Docker client could not be initiated. Make sure Docker is running.")
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
                    /data/{OUTPUT_PATHS["fasta"].name} > trimmed-genbank.fasta; \
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
            f"--metadata {OUTPUT_PATHS['metadata']} \\\n"
            f"--id3c-metadata {OUTPUT_PATHS['id3c-metadata'] if 'id3c-metadata' in OUTPUT_PATHS else OUTPUT_PATHS['lims-metadata']} \\\n"
            f"--fasta {OUTPUT_PATHS['fasta']} \\\n"
            f"--metrics {OUTPUT_PATHS['metrics']} \\\n"
            f"--nextclade {OUTPUT_PATHS['nextclade']} \\\n"
            f"--previous-submissions {OUTPUT_PATHS['previous-submissions']} \\\n"
            f"--excluded-vocs {OUTPUT_PATHS['excluded-vocs']} \\\n"
            f"--vadr-dir {OUTPUT_PATHS['vadr-dir']} \\\n"
            f"--output-dir {output_batch_dir} \\\n"
            f"--strain-id {next_avail_strain_id} \\\n"
            f"--test-name COVSEQ"
            f"--gisaid-username <your username> \n"
        )
    else:
        print("Warning: Missing required inputs to create submissions:")
        print('\n'.join([str(x) for x in OUTPUT_PATHS.values() if not x.exists()]))
