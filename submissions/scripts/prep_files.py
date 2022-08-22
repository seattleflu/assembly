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
        # The minus sign before m and d removes leading zeros, to match the format of the original Excel file.
        # Note: this may only work on unix, windows may need to use `#` instead of `-`
        try:
            return pd.to_datetime(date_value, format='%Y-%m-%d %H:%M:%S').strftime('%-m/%-d/%Y')
        except:
            raise Exception(f"Could not convert collection_date {date_value} to valid datetime.")

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

def validate_collection_dates(collection_dates):
    cutoff_datetime = datetime(2022, 2, 1)

    too_early_dates = collection_dates.dropna()[pd.to_datetime(collection_dates, format='%m/%d/%Y') < cutoff_datetime]

    if not too_early_dates.empty:
        raise Exception(f"Error: Collection date(s) earlier than {cutoff_datetime.strftime('%Y-%m-%d')}:\n {too_early_dates}")

    if len(collection_dates.dropna().unique()) < 2:
        raise Exception(f"Error: All collection dates values are the same:\n {collection_dates.unique()}")



def standardize_and_qc_external_metadata(metadata_filename):
    # Standardize and validate barcode and date values on external metadata Excel file. Loads 2 sheets to dataframes
    # and writes them back to Excel if data passes all checks.

    # Metadata sheet
    external_metadata_metadata = pd.read_excel(metadata_filename, sheet_name='Metadata')
    external_metadata_metadata['lab_accession_id'] = external_metadata_metadata['lab_accession_id'].apply(standardize_barcode)
    external_metadata_metadata['collection_date'] = external_metadata_metadata['collection_date'].apply(standardize_date)
    external_metadata_metadata['Seq date'] = external_metadata_metadata['Seq date'].apply(standardize_date)
    # Samplify FC Data sheet
    external_metadata_samplify_fc_data = pd.read_excel(args.metadata_file, sheet_name='Samplify FC Data', header=1)
    external_metadata_samplify_fc_data['Investigator\'s sample ID'] = external_metadata_samplify_fc_data['Investigator\'s sample ID'].apply(standardize_barcode)
    external_metadata_samplify_fc_data['Collection date'] = external_metadata_samplify_fc_data['Collection date'].apply(standardize_date)

    # Only the collection date from the Metadata sheet is used for creating submission files.
    # Validating to confirm no collection dates prior to Feb 2020 and not all dates are identical.
    validate_collection_dates(external_metadata_metadata['collection_date'])

    # Identify and optionally remove expirimental samples
    exp_samples = external_metadata_metadata[external_metadata_metadata['lab_accession_id'].str.endswith('_exp', na=False)]
    exp_samples_samplify = external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Investigator\'s sample ID'].str.endswith('_exp', na=False)]
    if not exp_samples.empty:
        print(f"Expirimental samples found: {len(exp_samples)}\n {exp_samples.to_string()}")

        if yes_no_cancel("Drop expirimental samples?"):
            external_metadata_metadata.drop(exp_samples.index, inplace=True)
            external_metadata_samplify_fc_data.drop(exp_samples_samplify.index, inplace=True)

    # Identify and optionally remove cascadia samples
    cascadia_samples = external_metadata_metadata[external_metadata_metadata['county'].str.lower().str.strip() == 'cascadia']
    cascadia_samples_samplify = external_metadata_samplify_fc_data[external_metadata_samplify_fc_data['Origin'].str.lower().str.strip() == 'cascadia']
    if not exp_samples.empty or not cascadia_samples_samplify.empty:
        print(f"Cascadia samples found:\n {exp_samples.to_string()} \n {cascadia_samples_samplify.to_string()}")

        if yes_no_cancel("Drop Cascadia samples?"):
            external_metadata_metadata.drop(cascadia_samples.index, inplace=True)
            external_metadata_samplify_fc_data.drop(cascadia_samples_samplify.index, inplace=True)

    # The Samplify FC Data sheet has an extra row before the headers. Stashing this value to reinsert the row when writing back to Excel.
    metadata_wb = load_workbook(args.metadata_file)
    samplify_fc_data_a1 = metadata_wb['Samplify FC Data']['A1'].value

    with pd.ExcelWriter(Path(output_batch_dir, 'external-metadata.xlsx'), engine='xlsxwriter', date_format='m/d/yyyy') as writer:
        external_metadata_metadata.to_excel(writer, sheet_name='Metadata', index=False)
        # put column headers on 2nd row of Samplify FC Data sheet, then populate first row with `samplify_fc_data_a1`` value
        external_metadata_samplify_fc_data.to_excel(writer, sheet_name='Samplify FC Data', startrow = 1, index=False)
        worksheet = writer.sheets['Samplify FC Data']
        worksheet.write(0, 0, samplify_fc_data_a1)
        # set format for Crt values to 2 decimals so whole numbers aren't displayed as integers (matching source file format)
        format1 = writer.book.add_format({'num_format': '0.00'})
        worksheet.set_column('O:O', None, format1)
        writer.save()

    LOG.debug(f"Saved external metadata file to {Path(output_batch_dir, 'external-metadata.xlsx')}")


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

    batch_dir_name = 'Batch-' + args.batch_date

    # check that output directory doesn't already exist
    output_batch_dir = Path(args.output_dir, batch_dir_name)
    if os.path.isdir(output_batch_dir):
        raise FileExistsError(f"Output batch folder already exists: {output_batch_dir}. Aborting.")
    else:
        output_batch_dir.mkdir(parents=True, exist_ok=True)

    # extract assembly results to temp fastq directory
    LOG.debug(f"Extracting compressed results file (.tar.gz) to: {output_batch_dir}")
    assemby_results_file = tarfile.open(args.results_file)
    assemby_results_file.extractall(output_batch_dir)
    assemby_results_file.close()

    standardize_and_qc_external_metadata(args.metadata_file)

    # create CSV of SFS sample barcodes
    identifiers = read_all_identifiers(Path(output_batch_dir, 'external-metadata.xlsx'))
    sfs_identifiers = find_sfs_identifiers(identifiers)

    if sfs_identifiers is not None:
        id3c_metadata_csv = Path(output_batch_dir, 'sfs-sample-barcodes.csv')
        sfs_identifiers.to_csv(id3c_metadata_csv, index=False)
        LOG.debug(f"Saved {len(sfs_identifiers)} SFS identifiers to {id3c_metadata_csv}.")

        if yes_no_cancel("Pull metadata from ID3C?"):
            LOG.debug("Pulling metadata from ID3C")

            # Get data files from NextClade
            id3c_metadata_with_county_csv = Path(output_batch_dir, 'id3c_metadata_with_county.csv')

            # The `export_id3c_metadata` bash script is in the same location as current script
            export_id3c_metadata_script = Path(Path(__file__).resolve().parent, "export_id3c_metadata")

            # Temporarily point stdout to a file while running this script, then point it back
            stdout = sys.stdout
            with open(id3c_metadata_with_county_csv, 'w') as sys.stdout:
                result = Conda.run_command('run', f"{export_id3c_metadata_script}",
                    f"{id3c_metadata_csv}",
                    stdout="STDOUT",
                    stderr="STDERR"
                )
            sys.stdout = stdout

            if result and len(result)==3 and result[2] == 0 and id3c_metadata_with_county_csv.exists() and id3c_metadata_with_county_csv.stat().st_size > 0:
                LOG.debug(f"Successfully saved ID3C metadata to {id3c_metadata_with_county_csv}")
            else:
                raise Exception(f"Error pulling metadata from ID3C.\n {result}")

            # Check for SFS samples with missing collection date
            id3c_metadata_df = pd.read_csv(id3c_metadata_with_county_csv)
            sfs_missing_date = id3c_metadata_df[(id3c_metadata_df['source']=='SFS')&(id3c_metadata_df['collection_date'].isna())]

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

    if yes_no_cancel("Process with NextClade?"):
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
        fasta_file = Path(output_batch_dir, Path(args.results_file).stem).with_suffix('.fa')
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
                ("S" in nextclade_df["failedGenes"].astype(str).str.split(",")) |
                (nextclade_df["aaChanges_s_total"] < 3) |
                (nextclade_df["qc.frameShifts.frameShifts"].str.split(",").apply(count_s_occurrences) > 0)]

        metadata = pd.read_excel(args.metadata_file)

        # remove twist positives from excluded VOCs
        twist_positives = metadata[metadata["lab_accession_id"].str.lower()=="twist positive"]
        excluded_vocs=pd.merge(excluded_vocs, twist_positives, on=["LIMS"], how="outer", indicator=True)
        excluded_vocs=excluded_vocs[excluded_vocs["_merge"]=="left_only"]

        LOG.debug(f'Excluded VOCs:{excluded_vocs[["LIMS", "seqName", "clade", "qc.missingData.status", "aaChanges_s_total", "qc.frameShifts.frameShifts", "failedGenes", "errors"]]}')

        excluded_vocs[['LIMS']].to_csv(Path(output_batch_dir, 'excluded-vocs.txt'), header=False, index=False)
        LOG.debug(f"Saved {len(excluded_vocs)} NWGC ids to {Path(output_batch_dir, 'excluded-vocs.txt')}.")


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