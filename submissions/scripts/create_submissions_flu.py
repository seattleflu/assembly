"""
Create files for submission of flu sequences to BioSample and reports for assembly summary.

Creates the following files:
- <batch-name>_metadata.csv: record of all metadata for potential debugging purposes
- <batch-name>_stample_stats.csv: summary of sample status to share with NWGC in SFS #assembly
- <batch-name>_biosample.tsv: TSV to upload to BioSample's submission portal
"""
import argparse
import sys
import pandas as pd
from datetime import datetime
from typing import List, Set, Optional
from pathlib import Path
from Bio import SeqIO
from utils import create_biosample_submission

base_dir = Path(__file__).resolve().parent.parent.parent
SUBMISSION_GROUPS = ['scan', 'sfs', 'wa-doh', 'cascadia']
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
    'pathogen',
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
    prev_subs_columns = ['phl_accession', 'sfs_sample_barcode', 'sfs_sample_identifier', 'status']
    prev_subs = pd.read_csv(previous_subs_file,sep='\t', dtype='string', usecols=prev_subs_columns)
    prev_subs = prev_subs.loc[prev_subs['status'].isin(['submitted', 'pending'])]
    return set(prev_subs['phl_accession'].combine_first(prev_subs['sfs_sample_barcode']).combine_first(prev_subs['sfs_sample_identifier']).tolist())



def parse_metadata(metadata_file: str, id3c_metadata_file: str = None, lims_metadata_file: str = None) -> pd.DataFrame:
    """
    Parse the required metadata from the provided *metadata_file* and the
    optional *id3c_metadata_file* and/or *lims_metadata_file*.

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

    #id3c_metadata, lims_metadata = None, None

    lims_id3c_metadata = pd.DataFrame()
    if id3c_metadata_file:
        id3c_metadata = pd.read_csv(id3c_metadata_file, dtype='string')

        # Convert PostgreSQL t/f values to Python boolean
        id3c_metadata['baseline_surveillance'] = id3c_metadata['baseline_surveillance'].map({'t': True, 'f': False})
        # add this column to match format of LIMS metadata
        id3c_metadata['sfs_identifier_for_doh_reporting'] = id3c_metadata['sfs_sample_barcode']

        lims_id3c_metadata = lims_id3c_metadata.append(id3c_metadata)

    if lims_metadata_file:
        lims_metadata = pd.read_csv(lims_metadata_file, dtype='string')

         # Convert PostgreSQL t/f values to Python boolean
        lims_metadata['baseline_surveillance'] = lims_metadata['baseline_surveillance'].map({'t': True, 'f': False})

        lims_id3c_metadata = lims_id3c_metadata.append(lims_metadata)

    if not lims_id3c_metadata.empty:
        # Drop collection date, county, swab_type, baseline_surveillance columns since they are expected to come from ID3C
        columns_to_drop = ['collection_date', 'county', 'swab_type', 'baseline_surveillance']
        external_metadata = metadata.drop(columns=columns_to_drop)

        # Merge ID3C metadata with the external metadata
        lims_id3c_metadata.drop(columns=['nwgc_id'], inplace=True)

        external_metadata['lab_accession_id'] = external_metadata['lab_accession_id'].str.lower()
        lims_id3c_metadata['sfs_sample_barcode'] = lims_id3c_metadata['sfs_sample_barcode'].str.lower()
        lims_id3c_metadata['sfs_collection_barcode'] = lims_id3c_metadata['sfs_collection_barcode'].str.lower()

        # match to metadata based on sample barcode, collection barcode, and concatenate into
        # single dataframe along with those that had no match (which may include control samples)
        metadata_1 = pd.merge(external_metadata, lims_id3c_metadata, how="inner", left_on="lab_accession_id", right_on='sfs_sample_barcode')
        metadata_2 = pd.merge(external_metadata, lims_id3c_metadata, how="inner", left_on="lab_accession_id", right_on='sfs_collection_barcode')
        metadata_missing = external_metadata[(~external_metadata['lab_accession_id'].isin(metadata_1['lab_accession_id']))&(~external_metadata['lab_accession_id'].isin(metadata_2['lab_accession_id']))]
        lims_id3c_metadata = pd.concat([metadata_1, metadata_2, metadata_missing], ignore_index=True)

        # Set lab accession id to the LIMS exported sample barcode to correct
        # barcodes with any Excel formatting issues. Use full sample identifier if barcode not available.
        lims_id3c_metadata.loc[lims_id3c_metadata['sfs_sample_barcode'].isna(), 'lab_accession_id'] = lims_id3c_metadata['sfs_sample_identifier']
        lims_id3c_metadata.loc[lims_id3c_metadata['sfs_sample_barcode'].notna(), 'lab_accession_id'] = lims_id3c_metadata['sfs_sample_barcode']

        lims_id3c_metadata['sequence_reason'] = lims_id3c_metadata['sequence_reason'].fillna(value='Sentinel surveillance')
        lims_id3c_metadata['originating_lab'] = 'Seattle Flu Study'
        lims_id3c_metadata['submission_group'] = 'sfs'
        # Label SCAN and Cascadia samples separately since they have different authors than SFS samples
        lims_id3c_metadata.loc[lims_id3c_metadata['source'].str.lower().str.strip() == 'scan', 'submission_group'] = 'scan'
        lims_id3c_metadata.loc[lims_id3c_metadata['source'].str.lower().str.strip() == 'cascadia', 'submission_group'] = 'cascadia'
    else:
        metadata['sfs_sample_identifier'] = None
        metadata['sfs_sample_barcode'] = None
        return metadata.sort_values('nwgc_id', ascending=True, ignore_index=True)

    # Sort by NWGC ID to ensure order of metadata
    return lims_id3c_metadata.sort_values('nwgc_id', ascending=True, ignore_index=True)


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
        state = metadata_row['state'] if not pd.isna(metadata_row['state']) else 'Washington'
        if pd.isna(location):
            metadata_row['state'] = 'Washington'
        elif pd.notna(state) and state == 'Washington' and location.lower() in washington_counties:
            metadata_row['county'] = location.title() + ' County'
        elif pd.notna(state) and state == 'Oregon' and location.lower() in oregon_counties:
            metadata_row['county'] = location.title() + ' County'
        elif location in manual_annotations['provided_location'].values:
            location_correction = manual_annotations.loc[manual_annotations['provided_location'] == location]
            metadata_row['state'] = location_correction['state'].values[0]
            metadata_row['county'] = location_correction['county'].values[0]  + ' County'
        else:
            print(f"Add unknown county to manual_location_annotations.tsv: {location} for state {metadata_row['state']}")

        return metadata_row

    metadata['collection_date'] = metadata['collection_date'].apply(standardize_date, args=('%Y-%m-%d',))

    washington_counties = text_to_list(base_dir / 'submissions/source-data/washington_counties.txt')
    washington_counties = [ county.lower() for county in washington_counties ]
    oregon_counties = text_to_list(base_dir / 'submissions/source-data/oregon_counties.txt')
    oregon_counties = [ county.lower() for county in oregon_counties ]

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


def add_assembly_metrics(metadata: pd.DataFrame, metrics_file: str, metadata_id: str = 'nwgc_id', metrics_id: str = 'nwgc_id') -> pd.DataFrame:
    """
    Parse assembly metrics from *metrics_file* and join with *metadata*.
    """
    assembly_metrics = parse_assembly_metrics(metrics_file)
    return metadata.merge(assembly_metrics, left_on=metadata_id, right_on=metrics_id, how='outer')


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
        #'EXTRA_COLUMN'
    ]

    # Map of needed metrics columns to final column names
    metrics_column_map = {
        'SampleId': 'segment_id',
        'MEAN_COVERAGE': 'coverage',
        'CONSENSUS_FASTA_LENGTH': 'length',
        'PCT_N_MAPPED': 'percent_ns',
    }

    assembly_metrics = pd.read_csv(metrics_file, sep='\t', dtype='string',
                                   usecols=metric_column_names)

    metrics = assembly_metrics[metrics_column_map.keys()].rename(columns=metrics_column_map)
    # get nwgc_id from segment_id
    metrics['nwgc_id'] = metrics['segment_id'].apply(lambda x: x.split('|')[0])
    return metrics


def add_clade_info(metadata: pd.DataFrame, nextclade_file: str, metadata_id: str = 'nwgc_id', nextclade_id_delimiter: str = None) -> pd.DataFrame:
    """
    Parse clade info from *nextclade_file* and join it with the corresponding
    sample in *metadata*.

    Assumes that the seqName column in the *nextclade_file* follows the format:
    Consensus_<nwgc_id>.consensus_threshold_*_quality_*

    Joins with metadata using the <nwgc_id>.
    """
    nextclade = pd.read_csv(nextclade_file, sep='\t')

    # NextClade changed the format of the `clade` column, but we'll continue using the old format
    # now in `clade_legacy`.
    if 'clade_legacy' in nextclade.columns:
        nextclade.drop(columns='clade', inplace=True)
        nextclade.rename(columns={'clade_legacy': 'clade'}, inplace=True)

    # Check format of seqName in NextClade file and parse accordingly
    if nextclade_id_delimiter is not None:
        nextclade['nwgc_id'] = nextclade['seqName'].apply(lambda s: s.split(nextclade_id_delimiter)[0])
    else:
        nextclade['nwgc_id'] = nextclade['seqName'].apply(lambda s: s.split(".")[0].split('_')[1])

    return metadata.merge(nextclade[['nwgc_id', 'clade']], left_on=metadata_id, right_on='nwgc_id', how='left')


def add_sequence_status(metadata: pd.DataFrame, prev_subs: Set[str], failed_nwgc_ids: List[str]) -> pd.DataFrame:
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

    # Find duplicate segments within current metadata set
    # Keep the segment with the lower percent Ns
    # Since flu is segmented and each segment has its own entry in the metadata df,
    # we need to consider both the lab accession id and the segment name
    # the segment name (HA, NA, etc.) needs to be parsed from the segment id first
    metadata['segment_name'] = metadata['segment_id'].apply(lambda x: x.split('|')[-1] if not pd.isnull(x) else x)
    current_duplicate = metadata.sort_values('percent_ns', ascending=True).duplicated(subset=['lab_accession_id','segment_name'], keep='first')
    
    # Find samples that have been previously submitted and mark as 'dropped duplicate'
    # this will not allow us to process segments from samples where some but not all segments have been previously submitted, but we probably won't want to do that anyway?
    overall_duplicate = metadata['lab_accession_id'].isin(prev_subs)
    
    # Label duplicate samples as 'dropped duplicate'
    metadata.loc[current_duplicate | overall_duplicate, 'status'] = 'dropped duplicate'

    # Label control samples
    metadata.loc[metadata['project'].str.lower() == 'sentinel', 'status'] = 'control'

    # Label samples that failed VADR
    metadata.loc[metadata['nwgc_id'].isin(failed_nwgc_ids), 'status'] = 'failed VADR'

    # Label samples with >10% Ns in the genome
    metadata.loc[pd.to_numeric(metadata['percent_ns'], errors='coerce') > 10, 'status'] = '>10% Ns'
    metadata.loc[pd.to_numeric(metadata['percent_ns'], errors='coerce') == 0, 'status'] = 'no Ns'

    # Find samples without a genome, i.e. 'length' is null
    no_genome = pd.isna(metadata['length'])

    # determine minimum length for flu.
    # Do we want to submit incomplete segments? Do we want to submit genomes that don't have
    # all the segments (or all complete segments)?
    # # Find samples with incomplete genomes, i.e. 'length' is less than 15000 (complete genome is ~15.2K)
    # incomplete_genome = pd.to_numeric(metadata['length'], errors='coerce') < 15000

    # Label samples that failed to generate genome
    # metadata.loc[no_genome | incomplete_genome, 'status'] = 'failed'
    metadata.loc[no_genome, 'status'] = 'failed'

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


def create_identifiers_report(metadata: pd.DataFrame, output_dir: Path, batch_name: str, pathogen: str) -> None:
    """
    Create a report that separates identifiers for easy tracking of samples.

    Follows the format of the identifiers TSV in
    https://github.com/seattleflu/sequence-identifiers
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

    identifiers['pathogen'] = pathogen
    identifiers[IDENTIFIER_COLUMNS].fillna('N/A').to_csv(output_dir / 'identifiers.tsv', sep='\t', index=False)


def create_summary_reports(metadata: pd.DataFrame, output_dir: Path, batch_name: str) -> None:
    """
    Creates summary reports of sample counts by source and batch's entire and HCT-specific date range.
    """
    all_samples = metadata
    counts_by_source = all_samples.groupby(['source']).size().reset_index(name='counts')
    counts_by_source.to_csv(output_dir / f'{batch_name}_sample_count_by_source.tsv', index=False, sep='\t')

    sample_date_range = all_samples.groupby(['source']).agg({'collection_date': ['min', 'max']})
    sample_date_range.columns = list(map(''.join, sample_date_range.columns.values))
    sample_date_range['source'] = sample_date_range.index

    sample_date_range = sample_date_range.append({'source': 'All',
        'collection_datemin': all_samples['collection_date'].dropna().min(),
        'collection_datemax': all_samples['collection_date'].dropna().max()
    }, ignore_index=True)

    sample_date_range.to_csv(output_dir / f'{batch_name}_sample_date_range.tsv', index=False, sep='\t')


def create_sample_status_report(metadata: pd.DataFrame, output_dir: Path, batch_name: str) -> None:
    """
    Create CSV from *metadata* that provides the final sample status for each
    sample in this sequencing batch. CSV to be sent to NWGC for QC purposes.
    """
    sample_status_columns = ['nwgc_id', 'status', 'percent_ns']
    metadata[sample_status_columns].to_csv(output_dir / f'{batch_name}_sample_status.csv', index=False)


# def create_gisaid_submission(metadata: pd.DataFrame, fasta: str, output_dir: Path,
#                              batch_name: str, submitter: str, test_name: str, subtype: str) -> None:
#     """
#     Create the CSV and FASTA files needed for submitting sequences to GISAID.
#     Follows the bulk upload CSV format required by GISAID.
#     """
#     def load_authors(project: str) -> str:
#         """
#         Pulls author names from `submissions/source-data/authors/{project}.txt`
#         and returs them in a single comma-separated string.
#         """
#         filename = base_dir / f'submissions/source-data/authors/{project}.txt'
#         authors_list = text_to_list(filename)

#         return ', '.join(authors_list)

#     def add_dynamic_fields(row: pd.Series, authors: dict) -> pd.Series:
#         """
#         Add GISAID columns to provided *row* based on row values.
#         This includes columns:
#         - rsv_virus_name
#         - rsv_location
#         - rsv_sampling_strategy
#         - rsv_coverage
#         - rsv_orig_lab_addr
#         - rsv_authors
#         """
#         # Create GISIAD location <region> / <country> / <state> / <county>
#         row['rsv_location'] = f'North America / USA / {row["state"]}'
#         if not pd.isna(row['county']):
#             row['rsv_location'] = row['rsv_location'] + ' / ' +  row['county']

#         # Match origin lab address based on originating lab
#         row['rsv_orig_lab_addr'] = lab_addresses.get(row['rsv_orig_lab'])
#         if row['rsv_orig_lab_addr'] is None:
#             sys.exit(f"Could not find lab address for {row['rsv_orig_lab']}. " \
#                       "Please add the address to `submissions/source-data/lab_addresses.tsv`")

#         # Assign authors based on submission group
#         assert row['submission_group'] in ['wa-doh', 'scan', 'sfs', 'cascadia', 'altius'], f'Invalid submission group: {row["submission_group"]}'
#         row['rsv_authors'] = authors[row['submission_group']]

#         # Add `hCoV-19/` prefix per GISAID requirements
#         row['rsv_virus_name'] = f'hCoV-19/{row["strain_name"]}'
#         # Follow GISAIDS requirement of reporting coverage like '100x'
#         row['rsv_coverage'] = row['coverage'].split('.')[0] + 'x'

#         # Follow CDC guidelines to label "Baseline surveillance" in `rsv_sampling_strategy` field
#         row['rsv_sampling_strategy'] = None
#         if row['baseline_surveillance'] == True:
#             row['rsv_sampling_strategy'] = 'Baseline surveillance'

#         return row

#     authors = { group: load_authors(group) for group in SUBMISSION_GROUPS }
#     lab_addresses = pd.read_csv(base_dir / 'submissions/source-data/lab_addresses.tsv', sep='\t', dtype='string')
#     lab_addresses = lab_addresses.set_index('lab')['address'].to_dict()
#     gisaid_file_base = f'SFS_{batch_name}_EpiRSV_BulkUpload'
#     gisaid_fasta = f'{gisaid_file_base}.fasta'

#     gisaid_columns = [
#         'submitter',
#         'fn',
#         'rsv_virus_name',
#         'rsv_type',
#         'rsv_passage',
#         'rsv_collection_date',
#         'rsv_location',
#         'rsv_add_location',
#         'rsv_host',
#         'rsv_add_host_info',
#         'rsv_sampling_strategy',
#         'rsv_gender',
#         'rsv_patient_age',
#         'rsv_patient_status',
#         'rsv_specimen',
#         'rsv_outbreak',
#         'rsv_last_vaccinated',
#         'rsv_treatment',
#         'rsv_seq_technology',
#         'rsv_assembly_method',
#         'rsv_coverage',
#         'rsv_orig_lab',
#         'rsv_orig_lab_addr',
#         'rsv_provider_sample_id',
#         'rsv_subm_lab',
#         'rsv_subm_lab_addr',
#         'rsv_subm_sample_id',
#         'rsv_authors'
#     ]

#     column_map = {
#         'collection_date': 'rsv_collection_date',
#         'originating_lab': 'rsv_orig_lab',
#     }

#     # Create a deep copy so manipulations don't affect subsequent GenBank submissions
#     gisaid_metadata = metadata.copy(deep=True)
#     gisaid_metadata.rename(columns=column_map, inplace=True)
#     # Ensure all GISAID columns are in the final DataFrame
#     gisaid_metadata = gisaid_metadata.reindex(columns=list(metadata.columns) + gisaid_columns)

#     # Hard-coded values
#     gisaid_metadata['submitter'] = submitter
#     gisaid_metadata['fn'] = gisaid_fasta
#     gisaid_metadata['rsv_subtype'] = subtype
#     gisaid_metadata['rsv_passage'] = 'Original'
#     gisaid_metadata['rsv_host'] = 'Human'
#     gisaid_metadata['rsv_gender'] = 'unknown'
#     gisaid_metadata['rsv_patient_age'] = 'unknown'
#     gisaid_metadata['rsv_patient_status'] = 'unknown'
#     if test_name=='MIPsSEQ':
#         gisaid_metadata['rsv_assembly_method'] = 'Brotman Baty Institute Viral MIP Sequencing'
#         gisaid_metadata['rsv_seq_technology'] = 'Illumina MiSeq'
#     else:
#         raise ValueError(f"Invalid test_name: {test_name}")
#     gisaid_metadata['rsv_subm_lab'] = SFS
#     gisaid_metadata['rsv_subm_lab_addr'] = lab_addresses[SFS]

#     # Dynamic values based on metadata for each sample
#     gisaid_metadata = gisaid_metadata.apply(add_dynamic_fields, args=(authors,), axis=1)

#     # Create the new GISAID FASTA after replacing the record ids with the rsv_virus_name
#     create_submission_fasta(fasta, gisaid_metadata, 'rsv_virus_name', output_dir / gisaid_fasta)

#     gisaid_metadata[gisaid_columns].to_csv(output_dir / f'{gisaid_file_base}.csv', index=False)


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
    parser.add_argument("--lims-metadata", type=str, required=False,
        help = "File path to the LIMS metadata CSV file" )
    parser.add_argument("--metrics", type=str, required=True,
        help = "File path to the TSV of assembly metrics from NWGC")
    parser.add_argument("--nextclade", type=str, required=True,
        help = "File path to the NextClade TSV file")
    parser.add_argument("--previous-submissions", type=str, required=True,
        help = "File path to TSV file containing previous submissions")
    parser.add_argument("--previous-submissions-sars-cov-2", type=str, required=True,
        help = "File path to TSV file containing previous SARS-Cov-2 submissions")
    parser.add_argument("--strain-id", type=int, required=True,
        help = "The starting numerical strain ID for this batch of sequences")
    parser.add_argument("--fasta", type=str, required=True,
        help = "File path to the FASTA file")
    parser.add_argument("--gisaid-username", type=str, required=True,
        help = "Submitter's GISAID username")
    parser.add_argument("--output-dir", type=str, required=True,
        help = "Path to the output directory for all output files")
    parser.add_argument("--pathogen", type=str, required=True,
        choices=['rsv-a', 'rsv-b', 'flu-a', 'flu-b'])
    parser.add_argument("--test-name", type=str, required=True,
        choices=['MIPsSEQ'],
        help = "Test name (MIPsSEQ)")

    args = parser.parse_args()

    # not running vadr for flu, so for now failed_nwgc_ids = []
    failed_nwgc_ids = []

    prev_subs = parse_previous_submissions(args.previous_submissions)

    metadata = parse_metadata(args.metadata, args.id3c_metadata, args.lims_metadata)
    metadata = standardize_metadata(metadata)

    metadata = add_assembly_metrics(metadata, args.metrics)
    metadata = add_clade_info(metadata, args.nextclade, nextclade_id_delimiter='|')

    metadata = add_sequence_status(metadata, prev_subs, failed_nwgc_ids)
    metadata = assign_strain_identifier(metadata, args.strain_id)

    metadata.to_csv(Path(args.output_dir) / 'metadata.csv', index=False)
    sys.exit()

    # Check for duplicated strain names
    previous_rsv_flu_strain_names = pd.read_csv(args.previous_submissions, sep='\t', usecols=['strain_name'])['strain_name'].dropna()
    previous_sars_cov_2_strain_names = pd.read_csv(args.previous_submissions_sars_cov_2, sep='\t', usecols=['strain_name'])['strain_name'].dropna()
    current_strain_names = metadata[metadata['strain_name']!='N/A']['strain_name'].dropna()
    all_strain_names = previous_rsv_flu_strain_names.append(previous_sars_cov_2_strain_names, ignore_index=True).append(current_strain_names, ignore_index=True)
    duplicate_strain_names = all_strain_names[all_strain_names.duplicated()]
    assert duplicate_strain_names.empty, f"Error: Overlapping identifiers\n {duplicate_strain_names}"

    # Create the output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    batch_name = args.batch_name

    # Create CSV with all metadata fields for easy debugging
    metadata.to_csv(output_dir / f'{batch_name}_metadata.csv', index=False)

    create_identifiers_report(metadata, output_dir, batch_name, args.pathogen)
    create_summary_reports(metadata, output_dir, batch_name)
    create_sample_status_report(metadata, output_dir, batch_name)

    # Only create submissions for sequences that have status "submitted"
    submit_metadata = metadata.loc[metadata['status'] == 'submitted']

    if submit_metadata.empty:
        sys.exit(f"No new submissions:\n {metadata.groupby(['status'])['status'].count()}")

    # Only create NCBI submissions for sequences that passed VADR
    ncbi_metadata = submit_metadata.loc[~submit_metadata['nwgc_id'].isin(failed_nwgc_ids)].copy(deep=True)

    biosample_metadata = create_biosample_submission(ncbi_metadata, output_dir, batch_name, args.pathogen)
    biosample_metadata.to_csv(output_dir / f'biosample_metadata.csv', index=False)
