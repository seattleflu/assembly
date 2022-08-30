"""
Create files for submission of SARS-CoV-2 sequences to GISAID, BioSample,
GenBank, WA DOH and reports for VoCs and assembly summary.

Creates the following files:
- <batch-name>_metadata.csv: record of all metadata for potential debugging purposes
- <batch-name>_stample_stats.csv: summary of sample status to share with NWGC in SFS #assembly
- <batch_name>_total_vocs.csv: count of all VoCs to share with SFS #sequencing
- <batch_name>_<study>_vocs.csv: summary of VoCs to share with each <study> in SFS Slack
- Batch_<batch_name>_sequencing_results.xlsx: summary of sequencing results to upload to WA DOH SFT
- Batch_<batch_name>_SCH_sequencing_results.xlsx: summary of SCH sequencing results to send to SCH
- SFS_<batch_name>_EpiCoV_BulkUpload.csv: CSV to upload to GISAID's bulk upload portal
- SFS_<batch_name>_EpiCoV_BulkUpload.fasta: FASTA to upload to GISAID's bulk upload portal
- <batch_name>_<submission_group>__genbank_metadata.tsv: TSV to upload to GenBank's submission portal
- <batch_name>_<submission_group>_genbank.fasta: FASTA to upload to GenBank's submission portal
"""
import argparse
import sys
import pandas as pd
from datetime import datetime
from typing import List, Set, Optional
from pathlib import Path
from Bio import SeqIO


base_dir = Path(__file__).resolve().parent.parent.parent
SUBMISSION_GROUPS = ['scan', 'sfs', 'wa-doh']
SFS = 'Seattle Flu Study'
IDENTIFIER_COLUMNS = [
    'nwgc_id',
    'batch',
    'strain_name',
    'phl_accession',
    'altius_sample_identifier',
    'sfs_sample_identifier',
    'sfs_sample_barcode',
    'status',
    'gisaid_accession',
    'genbank_accession',
    'clade',
    'pangolin',
]


def text_to_list(text_file: str) -> List[str]:
    """
    Convert *text_file* to a list of strings,
    where each line is a single element.
    """
    with open(text_file, 'r') as f:
        return [ line.strip() for line in f]


def parse_previous_submissions(previous_subs_file: str) -> Set[str]:
    """
    Parse identifiers of previously submitted sequences from the provided
    *previous_subs_file*

    The *previous_subs_file* is expected to be a TSV with the columns
    `phl_accession`, `sfs_sample_barcode`, and `status`, where a "submitted"
    status represents the sequence was submitted to GISAID/GenBank and a
    "pending" status represents the sequence was submitted but not yet accepted.

    Returns a list of sample ids for previously submitted sequences.
    """
    prev_subs_columns = ['phl_accession', 'sfs_sample_barcode', 'status']
    prev_subs = pd.read_csv(previous_subs_file,sep='\t', dtype='string', usecols=prev_subs_columns)
    prev_subs = prev_subs.loc[prev_subs['status'].isin(['submitted', 'pending'])]
    return set(prev_subs['phl_accession'].combine_first(prev_subs['sfs_sample_barcode']).tolist())


def parse_metadata(metadata_file: str, id3c_metadata_file: str = None) -> pd.DataFrame:
    """
    Parse the required metadata from the provided *metadata_file* and the
    optional *id3c_metadata_file*.

    The *metadata_file* is expected to be an XLSX file with a sheet named 'Metadata`
    and the optional *id3c_metadata_file* is expected to be a CSV file.
    """
    metadata_column_map = {
        'LIMS': 'nwgc_id',
        'Project': 'project',
        'lab_accession_id': 'lab_accession_id',
        'collection_date': 'collection_date',
        'submitting_lab': 'originating_lab',
        'sequence_reason': 'sequence_reason',
        'county': 'county',
    }
    metadata = pd.read_excel(metadata_file, engine='openpyxl', dtype='string',
                             sheet_name='Metadata', usecols=metadata_column_map.keys())
    metadata.rename(columns=metadata_column_map, inplace=True)
    metadata['submission_group'] = 'wa-doh'
    # Control samples are not going to be submitted
    metadata.loc[metadata['project'] == 'sentinel', 'submission_group'] = 'N/A'
    metadata['swab_type'] = 'unknown'
    # All WA DOH samples should be baseline surveillance samples
    metadata['baseline_surveillance'] = True

    if id3c_metadata_file:
        id3c_metadata = pd.read_csv(id3c_metadata_file, dtype='string')
        # Convert PostgreSQL t/f values to Python boolean
        d = {'t': True, 'f': False}
        id3c_metadata['baseline_surveillance'] = id3c_metadata['baseline_surveillance'].map(d)
        #id3c_metadata['baseline_surveillance'] = id3c_metadata['baseline_surveillance'].apply(
            #lambda x: True if x == 't' else False)

        # Find rows in external metdata that match ID3C samples
        external_metadata = metadata.loc[metadata['nwgc_id'].isin(id3c_metadata['nwgc_id'])]
        # Drop collection date, county, swab_type, baseline_surveillance columns since they are expected to come from ID3C
        columns_to_drop = ['collection_date', 'county', 'swab_type', 'baseline_surveillance']
        external_metadata = external_metadata.drop(columns=columns_to_drop)
        # Merge ID3C metadata with the external metadata
        id3c_metadata = id3c_metadata.merge(external_metadata, on=['nwgc_id'], how='left')

        # Set lab accession id to the ID3C exported sample barcode to correct
        # barcodes with any Excel formatting issues
        id3c_metadata['lab_accession_id'] = id3c_metadata['sfs_sample_barcode']
        id3c_metadata['sequence_reason'] = id3c_metadata['sequence_reason'].fillna(value='Sentinel surveillance')
        id3c_metadata['originating_lab'] = 'Seattle Flu Study'
        id3c_metadata['submission_group'] = 'sfs'
        # Label SCAN samples separately since they have different authors than SFS samples
        id3c_metadata.loc[id3c_metadata['source'] == 'SCAN', 'submission_group'] = 'scan'

        # Drop rows with nwgc_id in the ID3C metadata file
        # Ensures there are no duplicate rows in the final metadata DF
        metadata = metadata.loc[~metadata['nwgc_id'].isin(id3c_metadata['nwgc_id'])]
        metadata = pd.concat([metadata, id3c_metadata])

    else:
        metadata['sfs_sample_identifier'] = None
        metadata['sfs_sample_barcode'] = None

    # Sort by NWGC ID to ensure order of metadata
    return metadata.sort_values('nwgc_id', ascending=True, ignore_index=True)


def standardize_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize format of *metadata* collection date and location columns.
    """
    def standardize_location(metadata_row: pd.Series) -> pd.Series:
        """
        Standardize 'county' and 'state' columns within provided *metadata_row*.
        If no location data is provided, assume sequence is from Washington state
        """
        location = metadata_row['county']
        if pd.isna(location):
            metadata_row['state'] = 'Washington'
        elif location.lower() in washington_counties:
            metadata_row['state'] = 'Washington'
            metadata_row['county'] = location.title() + ' County'
        elif location in manual_annotations['provided_location'].values:
            location_correction = manual_annotations.loc[manual_annotations['provided_location'] == location]
            metadata_row['state'] = location_correction['state'].values[0]
            metadata_row['county'] = location_correction['county'].values[0]  + ' County'
        else:
            print(f"Add unknown county to manual_location_annotations.tsv: {location}")

        return metadata_row

    metadata['collection_date'] = metadata['collection_date'].apply(standardize_date, args=('%Y-%m-%d',))

    washington_counties = text_to_list(base_dir / 'submissions/source-data/washington_counties.txt')
    washington_counties = [ county.lower() for county in washington_counties ]
    manual_annotations = pd.read_csv( base_dir / 'submissions/source-data/manual_location_annotations.tsv',
                                     dtype='string', sep='\t')

    return metadata.apply(standardize_location, axis=1)


def standardize_date(date_string: str, format: str) -> Optional[str]:
    """
    Standardize provided *date_string* and return in *format*
    Returns None if unable to parse *date_string*
    """
    from dateutil import parser
    try:
        return parser.parse(date_string).strftime(format)
    except:
        return None


def add_assembly_metrics(metadata: pd.DataFrame, metrics_file: str) -> pd.DataFrame:
    """
    Parse assembly metrics from *metrics_file* and join with *metadata*.
    """
    assembly_metrics = parse_assembly_metrics(metrics_file)
    return metadata.merge(assembly_metrics, on='nwgc_id', how='left')


def parse_assembly_metrics(metrics_file: str) -> pd.DataFrame:
    """
    Parse out needed metrics data from the provided *metrics_file* TSV
    """
    # Expected metrics columns from the NWGC pipeline
    metric_column_names = [
        'SampleId',
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
        'CONTIG_NAME',
        'COUNT_A',
        'COUNT_C',
        'COUNT_G',
        'COUNT_T',
        'COUNT_N',
        'CONSENSUS_FASTA_LENGTH',
        'PCT_N_MAPPED',
        'REFERENCE_LENGTH',
        'PCT_N_REFERENCE',
        'EXTRA_COLUMN'
    ]

    # Map of needed metrics columns to final column names
    metrics_column_map = {
        'SampleId': 'nwgc_id',
        'MEAN_COVERAGE': 'coverage',
        'CONSENSUS_FASTA_LENGTH': 'length',
        'PCT_N_MAPPED': 'percent_ns',
    }

    # The TSV sometimes has an extra column at the end
    # Use header/skiprows/names options to ensure correct columns in DataFrame
    assembly_metrics = pd.read_csv(metrics_file, sep='\t', dtype='string',
                                   header=None, skiprows=1, names=metric_column_names)

    return assembly_metrics[metrics_column_map.keys()].rename(columns=metrics_column_map)


def add_clade_info(metadata: pd.DataFrame, nextclade_file: str) -> pd.DataFrame:
    """
    Parse clade info from *nextclade_file* and join it with the corresponding
    sample in *metadata*.

    Assumes that the seqName column in the *nextclade_file* follows the format:
    Consensus_<nwgc_id>.consensus_threshold_*_quality_*

    Joins with metadata using the <nwgc_id>.
    """
    nextclade = pd.read_csv(nextclade_file, sep='\t', usecols=['seqName', 'clade', 'Nextclade_pango'])
    nextclade['nwgc_id'] = nextclade['seqName'].apply(lambda s: s.split(".")[0].split('_')[1])
    nextclade.rename(columns={"Nextclade_pango": "pangolin"}, inplace=True)
    return metadata.merge(nextclade[['nwgc_id', 'clade', 'pangolin']], on='nwgc_id', how='left')


def add_sequence_status(metadata: pd.DataFrame, prev_subs: Set[str]) -> pd.DataFrame:
    """
    Add a column "status" to the provided *metadata* that reflects the final
    status of the each sequence. The provided *prev_subs* is a set of `lab_accession_id`
    for samples that have been previously submitted.

    Potential values for "status" in precedence order are:
    - failed
    - >10% Ns
    - control
    - dropped duplicate
    - missing collection date
    - submitted
    """
    # Default value is "submitted"
    metadata['status'] = 'submitted'

    # Label samples missing collection date
    metadata.loc[metadata['collection_date'].isnull(), 'status'] = 'missing collection date'

    # Find duplicate samples within current metadata set
    # Keep the sample with the lower percent Ns
    current_duplicate = metadata.sort_values('percent_ns', ascending=True).duplicated(subset='lab_accession_id', keep='first')
    # Find samples that have been previously submitted and mark as 'dropped duplicate'
    overall_duplicate = metadata['lab_accession_id'].isin(prev_subs)
    # Label duplicate samples as 'dropped duplicate'
    metadata.loc[current_duplicate | overall_duplicate, 'status'] = 'dropped duplicate'

    # Label control samples
    metadata.loc[metadata['project'] == 'sentinel', 'status'] = 'control'

    # Label samples with >10% Ns in the genome
    metadata.loc[pd.to_numeric(metadata['percent_ns'], errors='coerce') > 10, 'status'] = '>10% Ns'

    # Find samples without a genome, i.e. 'length' is null
    no_genome = pd.isna(metadata['length'])
    # Find samples with incomplete genomes, i.e. 'length' is less than 27000
    # Minimum length pulled from Nextstrain ncov: https://github.com/nextstrain/ncov/blob/master/defaults/parameters.yaml
    incomplete_genome = pd.to_numeric(metadata['length'], errors='coerce') < 27000
    # Label samples that failed to generate genome
    metadata.loc[no_genome | incomplete_genome, 'status'] = 'failed'

    return metadata


def assign_strain_identifier(metadata: pd.DataFrame, strain_id: int) -> pd.DataFrame:
    """
    For each sample that have status "submitted" in the provided *metadata*
    assign it a strain name, starting with the provided *strain_id*.

    Each strain name will have the format:
    USA/<state>-S<strain_id>/<collection_year>
    """
    def create_strain_name(row: pd.Series) -> str:
        """
        Create a strain name for the provided *row* from the metadata
        and return the strain name string.
        """
        if row.status != 'submitted':
            return 'N/A'

        collection_year = datetime.strptime(row.collection_date, '%Y-%m-%d').year
        state_abbr = state_abbreviations.get(row.state)

        if state_abbr is None:
            sys.exit(f"Could not find state «{row.state}» in us_state_abbreviations.csv")

        nonlocal strain_id
        strain_name = f"USA/{state_abbr}-S{strain_id}/{collection_year}"
        strain_id += 1

        return strain_name

    state_abbreviations = pd.read_csv(base_dir / 'submissions/source-data/us_state_abbreviations.csv')
    state_abbreviations = state_abbreviations.set_index('state')['abbreviation'].to_dict()

    metadata['strain_name'] = metadata.apply(create_strain_name, axis=1)
    return metadata


def create_identifiers_report(metadata: pd.DataFrame, output_dir: Path, batch_name: str) -> None:
    """
    Create a report that separates identifiers for easy tracking of samples.

    Follows the format of the identifiers TSV in
    https://github.com/seattleflu/hcov19-sequence-identifiers
    """
    def separate_identifiers(row: pd.Series) -> pd.Series:
        """
        Creates a separate column for `phl_accession` for the provided *row*
        and adds controls to the `sfs_sample_identifier` column.
        """
        if row['submission_group'] == 'wa-doh':
            row['phl_accession'] = row['lab_accession_id']
        else:
            row['phl_accession'] = 'N/A'

        if row['project'] == 'sentinel':
            row['sfs_sample_identifier'] = row['lab_accession_id']

        return row

    identifiers = metadata.copy(deep=True)
    identifiers['batch'] = batch_name + '_fastq'
    identifiers = identifiers.apply(separate_identifiers, axis=1)
    # We no longer sequence Altius samples
    identifiers['altius_sample_identifier'] = 'N/A'

    # These will be filled in after the sequences have been published
    identifiers['gisaid_accession'] = 'N/A'
    identifiers['genbank_accession'] = 'N/A'

    identifiers[IDENTIFIER_COLUMNS].fillna('N/A').to_csv(output_dir / 'identifiers.tsv', sep='\t', index=False)


def create_voc_reports(metadata: pd.DataFrame, excluded_vocs: str,
                       output_dir: Path, batch_name: str) -> None:
    """
    Finds variants of concern in the provided *metadata* to create a report of
    total VoCs and a CSV of VoCs for each study in SFS.

    Reports will be created in the provided *output_dir* with the *batch_name* prefix.

    Samples listed in the provided *exclude_vocs* are excluded from VoC reports
    """
    vocs = pd.read_csv(base_dir / 'submissions/source-data/variants_of_concern.tsv', sep='\t').rename(columns={'pangolin': 'clade_pangolin'})
    # Concatenate pangolin values for each clade. WHO values should be same per clade, so concatenate and remove duplicates
    vocs = vocs.groupby(['clade'], as_index = False).agg({'clade_pangolin': ', '.join, 'who': lambda x:', '.join(set(x))})
    exclude_ids = text_to_list(excluded_vocs)

    voc_samples = metadata.loc[(metadata['clade'].isin(vocs['clade'])) & (~metadata['nwgc_id'].isin(exclude_ids))]
    voc_samples = voc_samples.merge(vocs, on=['clade'], how='inner')

    all_vocs_counts = voc_samples.groupby(['clade', 'clade_pangolin', 'who']).size().reset_index(name='counts')
    all_vocs_counts.to_csv(output_dir / f'{batch_name}_total_vocs.csv', index=False)

    vocs_by_source = voc_samples.groupby(['clade', 'clade_pangolin', 'source']).size().reset_index(name='counts')
    vocs_by_source.to_csv(output_dir / f'{batch_name}_total_vocs_by_source.tsv', index=False, sep='\t')

    voc_report_columns = ['who', 'clade', 'clade_pangolin', 'pangolin', 'collection_date', 'sfs_sample_barcode', 'sfs_collection_barcode']
    sfs_vocs = voc_samples.loc[voc_samples['originating_lab'] == 'Seattle Flu Study']
    sfs_sources = sfs_vocs['source'].unique()

    for source in sfs_sources:
        source_vocs = sfs_vocs.loc[sfs_vocs['source'] == source][voc_report_columns]
        source_vocs.sort_values(by=['who'], inplace=True)
        source_vocs.to_csv(output_dir / f"{batch_name}_{source.replace('/', '_')}_vocs.csv", index=False)


def create_sample_status_report(metadata: pd.DataFrame, output_dir: Path, batch_name: str) -> None:
    """
    Create CSV from *metadata* that provides the final sample status for each
    sample in this sequencing batch. CSV to be sent to NWGC for QC purposes.
    """
    sample_status_columns = ['nwgc_id', 'status', 'percent_ns']
    metadata[sample_status_columns].to_csv(output_dir / f'{batch_name}_sample_status.csv', index=False)


def create_wa_doh_report(metadata:pd.DataFrame, pangolin: str, output_dir: Path, batch_name: str) -> None:
    """
    Create Excel report(s) of sequences for samples in *metadata*.
    Creates a separate file for SCH sequences if present because they need to
    be passed along to SCH to fill in the lab accession id.
    All other sequences are included in a single file to be submitted to
    WA DOH via SFT.

    Follows the WA DOH template with the following columns:
    - LAB_ACCESSION_ID: Accession or specimen ID
    - GISAID_ID: ID assigned to 'Virus Name' field
    - SPECIMEN_COLLECTION_DATE: Date of specimen collection. Format as MM/DD/YYYY.
    - SUBMITTING_LAB: Name of sequencing laboratory. This should be consistent across all submissions.
    - SEQUENCE_REASON: Reason for sequencing.
        Select only from the options provided:
            - Sentinel Surveillance
            - Suspected Reinfection
            - Suspected Vaccine Breakthrough
            - Outbreak
            - Other
    - SEQUENCE_STATUS: Values: Complete, Failed, Low Quality
    - PANGO_LINEAGE: Lineage written in NextClade nomenclature.
    - FIRST_NAME
    - LAST_NAME
    - MIDDLE_NAME
    - DOB
    - ALTERNATIVE_ID

    Note: PII columns are not filled in because we do not ingest PII.
    """
    doh_columns = [
        'LAB_ACCESSION_ID',
        'GISAID_ID',
        'SPECIMEN_COLLECTION_DATE',
        'SUBMITTING_LAB',
        'SEQUENCE_REASON',
        'SEQUENCE_STATUS',
        'PANGO_LINEAGE',
        'FIRST_NAME',
        'LAST_NAME',
        'MIDDLE_NAME',
        'DOB',
        'ALTERNATIVE_ID'
    ]

    status_map = {
        'submitted': 'Complete',
        '>10% Ns': 'Low Quality',
        'failed': 'Failed'
    }

    column_map = {
        'lab_accession_id': 'LAB_ACCESSION_ID',
        'strain_name': 'GISAID_ID',
        'collection_date': 'SPECIMEN_COLLECTION_DATE',
        'sequence_reason': 'SEQUENCE_REASON',
        'status': 'SEQUENCE_STATUS',
        'Nextclade_pango': 'PANGO_LINEAGE'
    }

    # PANGO lineages generated from running FASTA file through https://clades.nextstrain.org/
    pango_columns = ['seqName', 'Nextclade_pango']
    pango_lineages = pd.read_csv(pangolin, dtype='string', sep='\t', usecols=pango_columns)
    pango_lineages['nwgc_id'] = pango_lineages['seqName'].apply(parse_fasta_id)

    # Only submit sequences that are not controls, not duplicates, and not missing collection date
    not_control = metadata['originating_lab'] != 'sentinel'
    not_duplicate = metadata['status'] != 'dropped duplicate'
    not_missing_date = metadata['status'] != 'missing collection date'
    # WA DOH asked us only include sequences from Washington state as of Oct 12, 2021
    washington_state_only = metadata['state'] == 'Washington'
    doh_report = metadata.loc[(not_control) & (not_duplicate) & (not_missing_date) & (washington_state_only)].copy(deep=True)
    doh_report = doh_report.merge(pango_lineages, on=['nwgc_id'], how='left')

    # Convert date to MM/DD/YYYY format according to WA DOH template
    doh_report['collection_date'] = doh_report['collection_date'].apply(standardize_date, args=('%m/%d/%Y',))
    # Convert sequence status to WA DOH standardized values
    doh_report['status'] = doh_report['status'].apply(lambda x: status_map[x])
    # Only include PANGO lineages for completed sequences
    doh_report.loc[doh_report['status'] != 'Complete', 'Nextclade_pango'] = 'N/A'

    # Hard-coded values
    doh_report['SUBMITTING_LAB'] = 'NW Genomics'
    doh_report['FIRST_NAME'] = 'N/A'
    doh_report['LAST_NAME'] = 'N/A'
    doh_report['MIDDLE_NAME'] = 'N/A'
    doh_report['DOB'] = 'N/A'
    doh_report['ALTERNATIVE_ID'] = 'N/A'

    doh_report.rename(columns=column_map, inplace=True)

    # Get subset of SCH sequences if there were SFS sequences in the is batch
    sfs = doh_report.loc[doh_report['originating_lab'] == SFS]
    if len(sfs.index) > 0:
        sch_report = sfs.loc[sfs['source'] == 'SCH']
        if len(sch_report.index) > 0:
            # Provide the current lab accession id as the alternative id since
            # SCH will fill in the original lab accession id
            pd.set_option('mode.chained_assignment', None)
            sch_report.loc[:, 'ALTERNATIVE_ID'] = sch_report['LAB_ACCESSION_ID']
            sch_report.loc[:, 'LAB_ACCESSION_ID'] = None

            sch_report[doh_columns].to_excel(output_dir / f'Batch_{batch_name}_SCH_sequencing_results.xlsx',
                                             engine='openpyxl', index=False)

            # Remove the SCH sequences from the DOH report so we don't submit
            # the same sequences twice
            doh_report = doh_report[~doh_report.index.isin(sch_report.index)]

    doh_report[doh_columns].to_excel(output_dir / f'Batch_{batch_name}_sequencing_results.xlsx', engine='openpyxl', index=False)


def create_gisaid_submission(metadata: pd.DataFrame, fasta: str, output_dir: Path,
                             batch_name: str, submitter: str) -> None:
    """
    Create the CSV and FASTA files needed for submitting sequences to GISAID.
    Follows the bulk upload CSV format required by GISAID.
    """
    def load_authors(project: str) -> str:
        """
        Pulls author names from `submissions/source-data/authors/{project}.txt`
        and returs them in a single comma-separated string.
        """
        filename = base_dir / f'submissions/source-data/authors/{project}.txt'
        authors_list = text_to_list(filename)

        return ', '.join(authors_list)

    def add_dynamic_fields(row: pd.Series, authors: dict) -> pd.Series:
        """
        Add GISAID columns to provided *row* based on row values.
        This includes columns:
        - covv_virus_name
        - covv_location
        - covv_sampling_strategy
        - covv_coverage
        - covv_orig_lab_addr
        - covv_authors
        """
        # Create GISIAD location <region> / <country> / <state> / <county>
        row['covv_location'] = f'North America / USA / {row["state"]}'
        if not pd.isna(row['county']):
            row['covv_location'] = row['covv_location'] + ' / ' +  row['county']

        # Match origin lab address based on originating lab
        row['covv_orig_lab_addr'] = lab_addresses.get(row['covv_orig_lab'])
        if row['covv_orig_lab_addr'] is None:
            sys.exit(f"Could not find lab address for {row['covv_orig_lab']}. " \
                      "Please add the address to `submissions/source-data/lab_addresses.tsv`")

        # Assign authors based on submission group
        row['covv_authors'] = authors['wa-doh']
        if row['submission_group'] != 'wa-doh':

            if row['submission_group'] == 'scan':
                row['covv_authors'] = authors['scan']
            else:
                row['covv_authors'] = authors['sfs']

        # Add `hCoV-19/` prefix per GISAID requirements
        row['covv_virus_name'] = f'hCoV-19/{row["strain_name"]}'
        # Follow GISAIDS requirement of reporting coverage like '100x'
        row['covv_coverage'] = row['coverage'].split('.')[0] + 'x'

        # Follow CDC guidelines to label "Baseline surveillance" in `covv_sampling_strategy` field
        row['covv_sampling_strategy'] = None
        if row['baseline_surveillance'] == True:
            row['covv_sampling_strategy'] = 'Baseline surveillance'

        return row

    authors = { group: load_authors(group) for group in SUBMISSION_GROUPS }
    lab_addresses = pd.read_csv(base_dir / 'submissions/source-data/lab_addresses.tsv', sep='\t', dtype='string')
    lab_addresses = lab_addresses.set_index('lab')['address'].to_dict()
    gisaid_file_base = f'SFS_{batch_name}_EpiCoV_BulkUpload'
    gisaid_fasta = f'{gisaid_file_base}.fasta'

    gisaid_columns = [
        'submitter',
        'fn',
        'covv_virus_name',
        'covv_type',
        'covv_passage',
        'covv_collection_date',
        'covv_location',
        'covv_add_location',
        'covv_host',
        'covv_add_host_info',
        'covv_sampling_strategy',
        'covv_gender',
        'covv_patient_age',
        'covv_patient_status',
        'covv_specimen',
        'covv_outbreak',
        'covv_last_vaccinated',
        'covv_treatment',
        'covv_seq_technology',
        'covv_assembly_method',
        'covv_coverage',
        'covv_orig_lab',
        'covv_orig_lab_addr',
        'covv_provider_sample_id',
        'covv_subm_lab',
        'covv_subm_lab_addr',
        'covv_subm_sample_id',
        'covv_authors'
    ]

    column_map = {
        'collection_date': 'covv_collection_date',
        'originating_lab': 'covv_orig_lab',
    }

    # Create a deep copy so manipulations don't affect subsequent GenBank submissions
    gisaid_metadata = metadata.copy(deep=True)
    gisaid_metadata.rename(columns=column_map, inplace=True)
    # Ensure all GISAID columns are in the final DataFrame
    gisaid_metadata = gisaid_metadata.reindex(columns=list(metadata.columns) + gisaid_columns)

    # Hard-coded values
    gisaid_metadata['submitter'] = submitter
    gisaid_metadata['fn'] = gisaid_fasta
    gisaid_metadata['covv_type'] = 'betacoronavirus'
    gisaid_metadata['covv_passage'] = 'Original'
    gisaid_metadata['covv_host'] = 'Human'
    gisaid_metadata['covv_gender'] = 'unknown'
    gisaid_metadata['covv_patient_age'] = 'unknown'
    gisaid_metadata['covv_patient_status'] = 'unknown'
    gisaid_metadata['covv_seq_technology'] = 'Illumina Nextseq'
    gisaid_metadata['covv_assembly_method'] = 'Northwest Genomics Center Viral Pipeline'
    gisaid_metadata['covv_subm_lab'] = SFS
    gisaid_metadata['covv_subm_lab_addr'] = lab_addresses[SFS]

    # Dynamic values based on metadata for each sample
    gisaid_metadata = gisaid_metadata.apply(add_dynamic_fields, args=(authors,), axis=1)

    # Create the new GISAID FASTA after replacing the record ids with the covv_virus_name
    create_submission_fasta(fasta, gisaid_metadata, 'covv_virus_name', output_dir / gisaid_fasta)

    gisaid_metadata[gisaid_columns].to_csv(output_dir / f'{gisaid_file_base}.csv', index=False)


def parse_fasta_id(record_id: str) -> str:
    """
    Returns the nwgc_id embedded in the provided *record_id*

    Expects the *record_id* to have the format:
    `Consensus_<nwgc_id>.consensus_threshold_0.5_quality_20`
    """
    return record_id.split('.')[0].split('_')[1]


def create_submission_fasta(fasta: str, metadata: pd.DataFrame,
                            record_id_col: str, output_fasta: str,
                            tag_baseline: bool = False) -> None:
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
            record_metadata = metadata.loc[metadata['nwgc_id'] == nwgc_id].to_dict('records')
            # This should never happen!
            if len(record_metadata) > 1:
                sys.exit(f"Found multiple metadata records for NWGC ID «{nwgc_id}»")
            if len(record_metadata) == 1:
                record.id = record_metadata[0][record_id_col]
                record.description = ''
                # If required, add tag for baseline samples
                if (tag_baseline == True and
                    record_metadata[0]['baseline_surveillance'] == True):
                    record.description = '[keyword=purposeofsampling:baselinesurveillance]'

                SeqIO.write(record, output, 'fasta-2line')


def create_biosample_submission(metadata: pd.DataFrame, output_dir: Path, batch_name: str) -> pd.DataFrame:
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
        # Follows GenBank's requirement for isolate name for easier matching
        # of BioSample accession with GenBank record
        row['isolate'] = f"SARS-CoV-2/human/{row['sample_name']}"

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

    column_map = {
        'originating_lab': 'collected_by',
        'strain_name': 'sample_name',
    }

    # Rename columns according to BioSample template
    metadata.rename(columns=column_map, inplace=True)

    # Apply NCBI format requirements
    metadata = metadata.apply(format_fields, axis=1)

    # Hard-coded values
    # SFS BioProject accession (https://www.ncbi.nlm.nih.gov/bioproject/PRJNA746979)
    metadata['bioproject_accession'] = 'PRJNA746979'
    metadata['organism'] = 'Severe acute respiratory syndrome coronavirus 2'
    metadata['host'] = 'Homo sapiens'
    metadata['host_disease'] = 'COVID-19'
    metadata['isolation_source'] = 'clinical'

    metadata[biosample_columns].to_csv(output_dir / f'{batch_name}_biosample.tsv', sep='\t', index=False)
    return metadata


def create_genbank_submission(metadata: pd.DataFrame, fasta: str,
                              output_dir: Path, batch_name: str) -> None:
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

    for group in SUBMISSION_GROUPS:
        group_metadata = metadata.loc[metadata['submission_group'] == group]
        if len(group_metadata.index) > 0:
            output_base = output_dir / f'{batch_name}_{group}_genbank'
            group_metadata[genbank_columns].to_csv(f'{output_base}_metadata.tsv', sep='\t', index=False)
            create_submission_fasta(fasta, group_metadata, 'Sequence_ID', f'{output_base}.fasta', True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--batch-name", type=str, required=True,
        help = "The name for this batch of sequences, usually the date of sequence release, e.g. `20210701`")
    parser.add_argument("--metadata", type=str, required=True,
        help = "File path to the metadata Excel file")
    parser.add_argument("--id3c-metadata", type=str, required=False,
        help = "File path to the ID3C metadata CSV file")
    parser.add_argument("--metrics", type=str, required=True,
        help = "File path to the TSV of assembly metrics from NWGC")
    parser.add_argument("--nextclade", type=str, required=True,
        help = "File path to the NextClade TSV file")
    parser.add_argument("--previous-submissions", type=str, required=True,
        help = "File path to TSV file containing previous submissions")
    parser.add_argument("--strain-id", type=int, required=True,
        help = "The starting numerical strain ID for this batch of sequences")
    parser.add_argument("--fasta", type=str, required=True,
        help = "File path to the FASTA file")
    parser.add_argument("--gisaid-username", type=str, required=True,
        help = "Submitter's GISAID username")
    parser.add_argument("--output-dir", type=str, required=True,
        help = "Path to the output directory for all output files")
    parser.add_argument("--excluded-vocs", type=str,
        help = "File path to CSV of NWGC ids to exclude in from variants of concern")
    parser.add_argument("--vadr-dir", type=str, required=True,
        help = "Path to directory of VADR output files")

    args = parser.parse_args()

    # Verify that the provided VADR directory exists
    vadr_dir = Path(args.vadr_dir)
    if not vadr_dir.exists():
        sys.exit(f"ERROR: Provided VADR direcotry «{vadr_dir}» does not exist!")

    prev_subs = parse_previous_submissions(args.previous_submissions)

    metadata = parse_metadata(args.metadata, args.id3c_metadata)
    metadata = standardize_metadata(metadata)
    metadata = add_assembly_metrics(metadata, args.metrics)
    metadata = add_clade_info(metadata, args.nextclade)
    metadata = add_sequence_status(metadata, prev_subs)
    metadata = assign_strain_identifier(metadata, args.strain_id)

    # Create the output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_name = args.batch_name

    # Create CSV with all metadata fields for easy debugging
    metadata.to_csv(output_dir / f'{batch_name}_metadata.csv', index=False)

    create_identifiers_report(metadata, output_dir, batch_name)

    if args.id3c_metadata is not None:
        create_voc_reports(metadata, args.excluded_vocs, output_dir, batch_name)

    create_sample_status_report(metadata, output_dir, batch_name)
    create_wa_doh_report(metadata, args.nextclade, output_dir, batch_name)

    # Only create submissions for sequences that have status "submitted"
    submit_metadata = metadata.loc[metadata['status'] == 'submitted']

    if submit_metadata.empty:
        sys.exit(f"No new submissions:\n {metadata.groupby(['status'])['status'].count()}")

    create_gisaid_submission(submit_metadata, args.fasta, output_dir, batch_name, args.gisaid_username)

    # Only create NCBI submissions for sequences that passed VADR
    failed_nwgc_ids = [parse_fasta_id(id) for id in text_to_list(vadr_dir / 'genbank.vadr.fail.list')]
    ncbi_metadata = submit_metadata.loc[~submit_metadata['nwgc_id'].isin(failed_nwgc_ids)].copy(deep=True)

    biosample_metadata = create_biosample_submission(ncbi_metadata, output_dir, batch_name)
    create_genbank_submission(biosample_metadata, args.fasta, output_dir, batch_name)
