
import conda.cli.python_api as Conda
import logging
import os
import re
import docker
import sys
import requests
import pandas as pd
import getpass
from pathlib import Path
from openpyxl import load_workbook
from datetime import datetime
from Bio import SeqIO

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


def process_with_nextclade(nextclade_dataset_dir:Path, output_tsv:Path, input_fasta:Path, pathogen:str, subtype:str=None, segment:str=None):
    # nextclade uses dashes in `sars-cov-2` name but underscore for RSV
    if pathogen.startswith('rsv'):
        dataset = pathogen.replace('-', '_')
    elif pathogen.startswith('flu-a'):
        assert subtype in ['h1n1', 'h3n2'], f"Error: unknown influenza subtype {subtype}"
        if subtype == 'h1n1':
            subtype += 'pdm'
        segment = segment.lower()
        assert segment in ['ha','na'], f"Error: invalid segment {segment}"
        dataset = ('flu_' + subtype + '_' + segment)
        segment_fasta_output = Path(input_fasta.parents[0], input_fasta.stem + '.' + segment + '.fasta')
        with open(segment_fasta_output, 'w') as segment_fasta:
            for record in SeqIO.parse(input_fasta, 'fasta'):
                if record.id.split('|').pop() == segment.upper():
                    SeqIO.write(record, segment_fasta, 'fasta-2line')
        input_fasta = segment_fasta_output
    elif pathogen.startswith('flu-b'):
        # TODO: DRY up this portion of the code between flu-a and flu-b
        segment = segment.lower()
        assert segment in ['ha','na'], f"Error: invalid segment {segment}"        
        dataset = 'flu_vic_' + segment
        segment_fasta_output = Path(input_fasta.parents[0], input_fasta.stem + '.' + segment + '.fasta')
        with open(segment_fasta_output, 'w') as segment_fasta:
            for record in SeqIO.parse(input_fasta, 'fasta'):
                if record.id.split('|').pop() == segment.upper():
                    SeqIO.write(record, segment_fasta, 'fasta-2line')
        input_fasta = segment_fasta_output
    else:
        dataset = pathogen

    # Get data files from NextClade
    result = Conda.run_command('run', 'nextclade', 'dataset', 'get',
        f"--name={dataset}",
        f"--output-dir={nextclade_dataset_dir}"
    )

    if result and len(result)==3 and result[2] == 0:
        LOG.debug(f"Successfully downloaded {nextclade_dataset_dir} from NextClade.")
    else:
        raise Exception(f"Error: Nexclade {nextclade_dataset_dir} download failed.")


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
    

def pull_previous_submissions(output_tsv:Path, output_tsv_other:Path, pathogen:str) -> Path:
    SARS_COV_2_TSV_URL = "https://api.github.com/repos/seattleflu/hcov19-sequence-identifiers/contents/hcov19-sequence-identifiers.tsv"
    OTHER_TSV_URL = "https://api.github.com/repos/seattleflu/rsv-flu-sequence-identifiers/contents/rsv-flu-sequence-identifiers.tsv"

    token = os.environ.get('GH_ACCESS_TOKEN') or getpass.getpass('Github Personal access token:')

    url_output_map = [
        {
            'url': SARS_COV_2_TSV_URL if pathogen=='sars-cov-2' else OTHER_TSV_URL,
            'output': output_tsv
        },
        {
            'url': OTHER_TSV_URL if pathogen=='sars-cov-2' else SARS_COV_2_TSV_URL,
            'output': output_tsv_other
        }
    ]

    result = []
    for entry in url_output_map:
        LOG.debug(f"Downloading {entry['url']}\nSaving to {entry['output']}")
        try:
            r = requests.get(entry['url'],
                headers={
                    'accept': 'application/vnd.github.v3.raw',
                    'authorization': 'token {}'.format(token)
                }
            )

            if r.status_code == 200:
                with open(entry['output'], 'w') as f:
                    f.write(r.text)
                result.append(entry['output'])
            else:
                LOG.error(f"Error downloading previous submissions file: {r}")
                raise Exception

        except requests.exceptions.HTTPError as e:
            LOG.error(f"Error: {e.response.text}")
            raise Exception

    return result

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


def process_with_vadr(input_fasta:Path, output_batch_dir:Path, pathogen:str, interactive:bool=False) -> str:
    VADR_PARAMS = {
        'sars-cov-2': {
            'mkey': 'sarscov2',
            'minlen': 50,
            'maxlen': 30000,
            'annotate_flags': '--cpu 4 --glsearch -s -r --nomisc --lowsim5seq 6 --lowsim3seq 6 --alt_fail lowscore,insertnn,deletinn',
        },
        'rsv-a': {
            'mkey': 'rsv',
            'minlen': 150,
            'maxlen': 15500,
            'annotate_flags': '-r --xnocomp',
        },
        'rsv-b': {
            'mkey': 'rsv',
            'minlen': 150,
            'maxlen': 15500,
            'annotate_flags': '-r --xnocomp',
        }
    }

    vadr_pathogen_params = VADR_PARAMS[pathogen]

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
                    --minlen {vadr_pathogen_params["minlen"]} \
                    --maxlen {vadr_pathogen_params["maxlen"]} \
                    /data/{input_fasta.name} > trimmed-genbank-{pathogen}.fasta; \
                v-annotate.pl --split -f \
                    {vadr_pathogen_params["annotate_flags"]} \
                    --noseqnamemax \
                    --mkey {vadr_pathogen_params["mkey"]}  \
                    --mdir /opt/vadr/vadr-models \
                    /data/trimmed-genbank-{pathogen}.fasta \
                    /data/genbank-{pathogen}"',
            detach=True,
            remove=True,
            user=f"{os.getuid()}:{os.getgid()}",
            volumes= {output_batch_dir.resolve(): {'bind': '/data', 'mode': 'rw'}}
        )
        docker_output = docker_container.attach(stdout=True, stream=True, logs=True);
        for line in docker_output:
            LOG.debug(line.decode('utf-8'))

        LOG.debug(f"Finished VADR processing, output saved to: {Path(output_batch_dir, 'genbank-' + pathogen)}")
        return Path(output_batch_dir, 'genbank-' + pathogen)

    # should not get this far if `process_with_vadr` is true
    if process_with_vadr:
        raise Exception(f"Error: Something went wrong with VADR processing of {input_fasta.resolve()}")



def sequence_metrics(sequence: str, fasta_id: str, pathogen: str) -> dict:
    """
    Calculate basic metrics
    """
    n_count = sequence.count("N")
    seq_length = len(sequence)
    if pathogen=='sars-cov-2':
        ref_length = 29903  # Wuhan-Hu-1 complete genome
    elif pathogen=='rsv-a':
        ref_length = 15225  # RSV-A reference genome
    elif pathogen=='rsv-b':
        ref_length = 15222  # RSV-B reference genome
    elif pathogen=='flu-b':
        ref_length = 9124  # FLU-B reference genome (EPI_ISL_983345: B/Austria/1359417/2021)
    metrics = {
        'CONTIG_NAME': fasta_id,
        'COUNT_A': sequence.count("A"),
        'COUNT_C': sequence.count("C"),
        'COUNT_G': sequence.count("G"),
        'COUNT_T': sequence.count("T"),
        'COUNT_N': n_count,
        'CONSENSUS_FASTA_LENGTH': len(sequence),
        'PCT_N_MAPPED': (n_count / seq_length) * 100,
        'REFERENCE_LENGTH': ref_length,
        'PCT_N_REFERENCE':(n_count / ref_length) * 100,
    }
    return metrics

def validate_collection_dates(collection_dates: pd.DataFrame, batch_date: datetime):
    collection_dates = collection_dates[collection_dates['collection_date'].notna()]
    collection_dates['collected'] = pd.to_datetime(collection_dates['collection_date'], format='%m/%d/%Y')

    out_of_range = collection_dates[
        ((collection_dates['Project'].str.contains('_covid_'))&(collection_dates['collected'] < datetime(2020, 2, 1)))|
        (collection_dates['collected'] > batch_date)
    ]

    if not out_of_range.empty:
        raise Exception(f"Error: Collection date(s) out of range):\n {out_of_range}")

    if len(collection_dates['collected'].unique()) < 2:
        raise Exception(f"Error: All collection dates values are the same:\n {collection_dates['collected'].unique()}")

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
        bool(re.match('^.*_exp$', standardized_barcode)) or
        standardized_barcode in ["twist positive", "blank"]):
        return standardized_barcode
    else:
        raise ValueError(f"Identifiers must be hexidecimal barcodes, `blank`, or `twist positive`. Invalid barcode: {standardized_barcode}")

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

def standardize_and_qc_external_metadata(metadata_filename:str, batch_date:datetime, metadata_output:Path):
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
    external_metadata_samplify_fc_data = pd.read_excel(metadata_filename, sheet_name='Samplify FC Info', header=1, converters={'Sample ID':int, 'Investigator\'s sample ID':str})
    external_metadata_samplify_fc_data['Investigator\'s sample ID'] = external_metadata_samplify_fc_data['Investigator\'s sample ID'].str.strip()
    # only standardize internal barcodes
    external_metadata_samplify_fc_data['Investigator\'s sample ID'] = external_metadata_samplify_fc_data.apply(lambda row: standardize_barcode(row['Investigator\'s sample ID']) if str(row['Project']).startswith(('starita_bbi_', 'SeattleChildrensDirect_')) else row['Investigator\'s sample ID'], axis=1)
    external_metadata_samplify_fc_data['Collection date'] = external_metadata_samplify_fc_data['Collection date'].apply(standardize_date)

    # Only the collection date from the Metadata sheet is used for creating submission files.
    # Validating to confirm no collection dates prior to Feb 2020 for `sars-cov-2` and not all dates are identical.
    validate_collection_dates(external_metadata_metadata[['collection_date','Project']], batch_date)

    sentinel_ids = [
        'blank',
        'twist positive',
        'pbs',
        'lp_blank',
        'sars-cov-2 control',
        'rsva control',
        'rsvb control',
        'flua control',
        'flub control',
        'h3n2 control',
        'h1n1 control'
    ]

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
        print(f"Experimental samples found: {len(exp_samples)}\n {exp_samples.to_string()}")

        if yes_no_cancel("Drop experimental samples?"):
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
        else:
            # Cascadia samples should be marked as sequence_reson: Other
            external_metadata_metadata.loc[cascadia_samples.index, 'sequence_reason'] = 'Other'

    # The Samplify FC Data sheet has an extra row before the headers. Stashing this value to reinsert the row when writing back to Excel.
    metadata_wb = load_workbook(metadata_filename)
    samplify_fc_data_a1 = metadata_wb['Samplify FC Info']['A1'].value

    with pd.ExcelWriter(metadata_output, engine='xlsxwriter', date_format='m/d/yyyy') as writer:
        external_metadata_metadata.to_excel(writer, sheet_name='Metadata', index=False)
        # put column headers on 2nd row of Samplify FC Data sheet, then populate first row with `samplify_fc_data_a1`` value
        external_metadata_samplify_fc_data.to_excel(writer, sheet_name='Samplify FC Info', startrow = 1, index=False)
        worksheet = writer.sheets['Samplify FC Info']
        worksheet.write(0, 0, samplify_fc_data_a1)
        # set format for Crt values to 2 decimals so whole numbers aren't displayed as integers (matching source file format)
        format1 = writer.book.add_format({'num_format': '0.00'})
        worksheet.set_column('O:O', None, format1)
        writer.save()

    LOG.debug(f"Saved external metadata file to {metadata_output}")


def create_biosample_submission(metadata: pd.DataFrame, output_dir: Path, batch_name: str, pathogen: str) -> pd.DataFrame:
    """
    Create the TSV for submission to BioSample through the
    NCBI webste (https://submit.ncbi.nlm.nih.gov/subs/biosample/)
    """
    def format_fields(row: pd.Series) -> pd.Series:
        """
        Add BioSample columns with specific format requirements for each *row*

        Columns include:
        - geo_loc_name
        - isolate
        """
        PATHOGEN_PREFIX = {
            'sars-cov-2': 'SARS-CoV-2',
            'rsv-a': 'RSV-A',
            'rsv-b': 'RSV-B',
            'flu-a': 'Influenza-A',
            'flu-b': 'Influenza-B',
        }
        # Follows GenBank's requirement for isolate name for easier matching
        # of BioSample accession with GenBank record
        row['isolate'] = f"{PATHOGEN_PREFIX[pathogen]}/human/{row['sample_name']}"

        # Follows NCBI's requirement of reporting location as
        # <country>:<state>,<county>
        row['geo_loc_name'] = f"USA:{row.state}"
        if not pd.isna(row['county']):
            row['geo_loc_name'] = row['geo_loc_name'] + ',' + row['county']

        return row

    biosample_columns = [
        'sample_name',
        'bioproject_accession',
        'organism',
        'collected_by',
        'collection_date',
        'geo_loc_name',
        'host',
        'host_disease',
        'isolate',
        'isolation_source'
    ]

    # Pathogen submissions (excluding SARS-CoV-2) require lat_long field
    if pathogen!='sars-cov-2':
        biosample_columns.append('lat_lon')

    column_map = {
        'originating_lab': 'collected_by',
        'strain_name': 'sample_name',
    }

    # Rename columns according to BioSample template
    metadata.rename(columns=column_map, inplace=True)

    # Apply NCBI format requirements
    metadata = metadata.apply(format_fields, axis=1)

    PATHOGEN_VALUES = {
        'sars-cov-2': {
            'bioproject_accession': 'PRJNA746979', # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA746979
            'organism': 'Severe acute respiratory syndrome coronavirus 2',
            'host_disease': 'COVID-19',
            'isolation_source': 'clinical',
        },
        'rsv-a': {
            'bioproject_accession': 'PRJNA988152', # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA988152
            'organism': 'Human respiratory syncytial virus A',
            'host_disease': 'Human respiratory syncytial virus A',
            'isolation_source': 'nasal swab',
            'lat_lon': 'missing',
        },
        'rsv-b': {
            'bioproject_accession': 'PRJNA988154', # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA988154
            'organism': 'Human respiratory syncytial virus B',
            'host_disease': 'Human respiratory syncytial virus B',
            'isolation_source': 'nasal swab',
            'lat_lon': 'missing',
        },
        'flu-a': {
            'bioproject_accession': 'PRJNA988813', # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA988813
            'organism': 'Influenza A',
            'host_disease': 'Influenza A',
            'isolation_source': 'nasal swab',
            'lat_lon': 'missing',
        },
        'flu-b': {
            'bioproject_accession': 'PRJNA988156', # https://www.ncbi.nlm.nih.gov/bioproject/PRJNA988156
            'organism': 'Influenza B',
            'host_disease': 'Influenza B',
            'isolation_source': 'nasal swab',
            'lat_lon': 'missing'
        },
    }

    # SFS BioProject accession (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA746979)
    assert pathogen in PATHOGEN_VALUES.keys(), f"Unknown pathogen, no known BioProject: {pathogen}"

    for k,v in PATHOGEN_VALUES[pathogen].items():
        metadata[k] = v

    # Hard-coded values
    metadata['host'] = 'Homo sapiens'

    metadata[biosample_columns].to_csv(output_dir / f'{batch_name}_biosample.tsv', sep='\t', index=False)
    return metadata


def create_genbank_submission(metadata: pd.DataFrame, fasta: str,
                              output_dir: Path, batch_name: str, pathogen: str = None) -> None:
    """
    Create the TSV and FASTA files necessary for submissions
    through the NCBI website (https://submit.ncbi.nlm.nih.gov/subs/genbank/)

    Creates a separate TSV + FASTA file for each project (scan, sfs, wa-doh)
    since they have different list of authors.

    Expects the provided *metadata* to contain the BioSample columns created
    by `create_biosample_submission`.
    """
    genbank_columns = [
        'Sequence_ID',
        'isolate',
        'country',
        'host',
        'collection-date',
        'isolation-source'
    ]

    column_map = {
        'geo_loc_name': 'country',
        'collection_date': 'collection-date',
        'isolation_source': 'isolation-source'
    }

    # Rename BioSample columns to match GenBank template columns
    metadata.rename(columns=column_map, inplace=True)

    # Within the sample_name USA/<state>-<strain_id>/<year>
    # the Sequence ID is the <state>-<strain_id>
    metadata['Sequence_ID'] = metadata['sample_name'].apply(
        lambda x: x.split('/')[1]
    )

    SUBMISSION_GROUPS = ['scan', 'sfs', 'wa-doh', 'cascadia']
    for group in SUBMISSION_GROUPS:
        group_metadata = metadata.loc[metadata['submission_group'] == group]
        if len(group_metadata.index) > 0:
            output_base = output_dir / f'{batch_name}_{group}_genbank'
            group_metadata[genbank_columns].to_csv(f'{output_base}_metadata.tsv', sep='\t', index=False)
            create_submission_fasta(fasta, group_metadata, 'Sequence_ID', f'{output_base}.fasta', True, pathogen)


def create_submission_fasta(fasta: str, metadata: pd.DataFrame,
                            record_id_col: str, output_fasta: str,
                            tag_baseline: bool = False, pathogen: str = None) -> None:
    """
    Create a new FASTA file *output_fasta* by filtering for sequences that have
    `nwgc_id` that are included in the provided *metadata*.
    Replaces the original FASTA record id with the corresponding
    *record_id_col* in the *metadata*.

    If *tag_baseline* is True, then add baseline tag to the sequences
    according to CDC guidelines:
        >SeqID [keyword=purposeofsampling:baselinesurveillance]
    """
    with open(fasta, 'r') as original, open(output_fasta, 'w') as output:
        for record in SeqIO.parse(original, 'fasta'):
            nwgc_id = parse_fasta_id(record.id)
            record_metadata = metadata.loc[metadata['nwgc_id'].astype(str) == nwgc_id].to_dict('records')
            # This should never happen!
            if len(record_metadata) > 1:
                sys.exit(f"Found multiple metadata records for NWGC ID «{nwgc_id}»")
            if len(record_metadata) == 1:
                record.id = record_metadata[0][record_id_col]
                record.description = ""
                if pathogen.startswith('rsv-'):
                    if 'bioproject_accession' in record_metadata[0]:
                        record.description += f" [BioProject={record_metadata[0]['bioproject_accession']}] [BioSample={record_metadata[0]['accession']}] [organism=Human respiratory syncytial virus {pathogen.split('-').pop().upper()}]"
                    else:
                        sys.exit(f"BioProject and BioSample accessions required for RSV sequences")

                # If required, add tag for baseline samples
                if (tag_baseline == True and
                    record_metadata[0]['baseline_surveillance'] == True):
                    record.description += ' [keyword=purposeofsampling:baselinesurveillance]'
                SeqIO.write(record, output, 'fasta-2line')

def parse_fasta_id(record_id: str) -> str:
    """
    Returns the nwgc_id embedded in the provided *record_id*

    Expects the *record_id* to have the format:
    `Consensus_<nwgc_id>.consensus_threshold_0.5_quality_20`
    or
    `<nwgc_id>|...`
    """
    if '|' in record_id:
        return record_id.split('|')[0]
    else:
        return record_id.split('.')[0].split('_')[1]
