"""
Create files for submissions to GISAID, GenBank, WA DOH and reports for VoCs
and assembly summary.

Creates the following files:
- <batch-name>_metadata.csv: record of all metadata for potential debugging purposes
- <batch-name>_stample_stats.csv: summary of sample status to share with NWGC in SFS #assembly
- <batch_name>_total_vocs.csv: count of all VoCs to share with SFS #sequencing
- <batch_name>_<study>_vocs.csv: summary of VoCs to share with each <study> in SFS Slack
- Batch_<batch_name>_sequencing_results.xlsx: summary of sequencing results to upload to WA DOH SFT
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
SFS_ADDRESS = 'University of Washington Medical Center, 1959 NE Pacific Street, Seattle, WA 98195, USA'


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
    status represents the sequence was submitted to GISAID/GenBank.

    Returns a list of sample ids for previously submitted sequences.
    """
    prev_subs_columns = ['phl_accession', 'sfs_sample_barcode', 'status']
    prev_subs = pd.read_csv(previous_subs_file,sep='\t', dtype='string', usecols=prev_subs_columns)
    prev_subs = prev_subs.loc[prev_subs['status'] == 'submitted']
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

    if id3c_metadata_file:
        id3c_metadata = pd.read_csv(id3c_metadata_file, dtype='string')

        # Find rows in external metdata that match ID3C samples
        external_metadata = metadata.loc[metadata['nwgc_id'].isin(id3c_metadata['nwgc_id'])]
        # Drop collection date and county columns since they are expected to come from ID3C
        external_metadata = external_metadata.drop(columns=['collection_date', 'county'])
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
        metadata = metadata.append(id3c_metadata)

    # Sort by NWGC ID to ensure order of metadata
    return metadata.sort_values('nwgc_id', ascending=True, ignore_index=True)


def standardize_metadata(metadata: pd.DataFrame) -> pd.DataFrame:
    """
    Standardize format of *metadata* collection date and location columns.
    """
    def standardize_date(date_string: str) -> Optional[str]:
        """
        Standardize provided *date_string* and return in format %Y-%m-%d
        Returns None if unable to parse *date_string*
        """
        from dateutil import parser
        try:
            return parser.parse(date_string).strftime('%Y-%m-%d')
        except:
            return None


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

    metadata['collection_date'] = metadata['collection_date'].apply(standardize_date)

    washington_counties = text_to_list(base_dir / 'submissions/source-data/washington_counties.txt')
    washington_counties = [ county.lower() for county in washington_counties ]
    manual_annotations = pd.read_csv( base_dir / 'submissions/source-data/manual_location_annotations.tsv', sep='\t')

    return metadata.apply(standardize_location, axis=1)


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
    nextclade = pd.read_csv(nextclade_file, sep='\t', usecols=['seqName', 'clade'])
    nextclade['nwgc_id'] = nextclade['seqName'].apply(lambda s: s.split(".")[0].split('_')[1])

    return metadata.merge(nextclade[['nwgc_id', 'clade']], on='nwgc_id', how='left')


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
    Create a repor that separates identifiers for easy tracking of samples.

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

    identifier_columns = [
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
    ]

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
    vocs = pd.read_csv(base_dir / 'submissions/source-data/variants_of_concern.tsv', sep='\t')
    exclude_ids = text_to_list(excluded_vocs)

    voc_samples = metadata.loc[(metadata['clade'].isin(vocs['clade'])) & (~metadata['nwgc_id'].isin(exclude_ids))]
    voc_samples = voc_samples.merge(vocs, on=['clade'], how='inner')

    all_vocs_counts = voc_samples.groupby(['who', 'pangolin']).size().reset_index(name='counts')
    all_vocs_counts.to_csv(output_dir / f'{batch_name}_total_vocs.csv', index=False)

    voc_report_columns = ['who', 'pangolin', 'collection_date', 'sfs_sample_barcode', 'sfs_collection_barcode']
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


def create_wa_doh_report(metadata:pd.DataFrame, output_dir: Path, batch_name: str) -> None:
    """
    Create an Excel report of sequences for samples in *metadata*.

    Follows the WA DOH template with the following columns:
    - LAB_ACCESSION_ID: Accession or specimen ID
    - GISAID_ID: ID assigned to 'Virus Name' field
    - COLLECTION_DATE: Date of specimen collection
    - SUBMITTING_LAB: Name of sequencing laboratory. This should be consistent across all submissions.
    - SEQUENCE_REASON: Reason for sequencing.
        Values: S-dropout, Suspected reinfection, Suspected vaccine breakthrough,
                Sentinel surveillance, Outbreak, Other
    - SEQUENCE_QUALITY: Values: Complete, Pending, Failed, Not Done, Low Quality
    """
    status_map = {
        'submitted': 'Complete',
        'dropped duplicate': 'Complete',
        '>10% Ns': 'Low Quality',
        'missing collection date': 'Pending',
        'failed': 'Failed'
    }

    column_map = {
        'lab_accession_id': 'LAB_ACCESSION_ID',
        'strain_name': 'GISAID_ID',
        'collection_date': 'COLLECTION_DATE',
        'originating_lab': 'SUBMITTING_LAB',
        'sequence_reason': 'SEQUENCE_REASON',
        'status': 'SEQUENCE_QUALITY'
    }

    doh_report = metadata.loc[metadata['originating_lab'] != 'sentinel'][column_map.keys()].copy(deep=True)
    doh_report.loc[doh_report['collection_date'].isnull(), 'collection_date'] = 'N/A'
    doh_report['status'] = doh_report['status'].apply(lambda x: status_map[x])
    doh_report.rename(columns=column_map, inplace=True)
    doh_report.to_excel(output_dir / f'Batch_{batch_name}_sequencing_results.xlsx', engine='openpyxl', index=False)


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
        - covv_coverage
        - covv_orig_lab_addr
        - covv_authors
        """
        # Create GISIAD location <region> / <country> / <state> / <county>
        row['covv_location'] = f'North America / USA / {row["state"]}'
        if not pd.isna(row['county']):
            row['covv_location'] = row['covv_location'] + ' / ' +  row['county']

        # Assign origin lab address and authors based on submission group
        row['covv_authors'] = authors['wa-doh']
        row['covv_orig_lab_addr'] = 'unknown'
        if row['submission_group'] != 'wa-doh':
            row['covv_orig_lab_addr'] = SFS_ADDRESS
            if row['submission_group'] == 'scan':
                row['covv_authors'] = authors['scan']
            else:
                row['covv_authors'] = authors['sfs']

        # Add `hCoV-19/` prefix per GISAID requirements
        row['covv_virus_name'] = f'hCoV-19/{row["strain_name"]}'
        # Follow GISAIDS requirement of reporting coverage like '100x'
        row['covv_coverage'] = row['coverage'].split('.')[0] + 'x'

        return row

    authors = { group: load_authors(group) for group in SUBMISSION_GROUPS }
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
        'sequence_reason': 'covv_sampling_strategy'
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
    gisaid_metadata['covv_subm_lab_addr'] = SFS_ADDRESS

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
                            record_id_col: str, output_fasta: str) -> None:
    """
    Create a new FASTA file *output_fasta* by filtering for sequences that have
    `nwgc_id` that are included in the provided *metadata*.
    Replaces the original FASTA record id with the corresponding
    *record_id_col* in the *metadata*.
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
                SeqIO.write(record, output, 'fasta-2line')


def create_genbank_submission(metadata: pd.DataFrame, fasta: str, vadr_dir: str,
                              output_dir: Path, batch_name: str) -> None:
    """
    Create the TSV and FASTA files necessary for submissions
    through the NCBI website (https://submit.ncbi.nlm.nih.gov/subs/genbank/)

    Creates a separate TSV + FASTA file for each project (scan, sfs, wa-doh)
    since they have different list of authors.

    Excludes sequences that have been failed by the VADR program listed in
    `*vadr_dir*/genbank.vadr.fail.list`
    """
    def add_dynamic_fields(row: pd.Series) -> pd.Series:
        """
        Add GenBank columns based on columns in provided *row*
        Columns include:
        - Sequence_ID
        - isolate
        - country
        """
        # Within the strain name USA/<state>-<strain_id>/<year>
        # the Sequence ID is the <state>-<strain_id>
        row['Sequence_ID'] = row['strain_name'].split('/')[1]

        # Follows GenBank's requirement for isolate name
        row['isolate'] = f"SARS-CoV-2/human/{row['strain_name']}"

        # Follows GenBank's requirement of reporting location as
        # <country>:<state>,<county>
        row['country'] = f"USA:{row.state}"
        if not pd.isna(row['county']):
            row['country'] = row['country'] + ',' + row['county']

        return row

    genbank_columns = [
        'Sequence_ID',
        'isolate',
        'country',
        'host',
        'collection-date',
        'isolation-source'
    ]

    column_map = {
        'collection_date': 'collection-date',
        'swab_type': 'isolation-source'
    }

    failed_nwgc_ids = [parse_fasta_id(id) for id in text_to_list(vadr_dir / 'genbank.vadr.fail.list')]

    genbank_metadata = metadata.loc[~metadata['nwgc_id'].isin(failed_nwgc_ids)].copy(deep=True)
    genbank_metadata.rename(columns=column_map, inplace=True)

    # Hard-coded values
    genbank_metadata['host'] = 'Homo sapiens'
    # Dynamic values based on metadata for each sample
    genbank_metadata = genbank_metadata.apply(add_dynamic_fields, axis=1)

    for group in SUBMISSION_GROUPS:
        group_metadata = genbank_metadata.loc[genbank_metadata['submission_group'] == group]
        if len(group_metadata.index) > 0:
            output_base = output_dir / f'{batch_name}_{group}_genbank'
            group_metadata[genbank_columns].to_csv(f'{output_base}_metadata.tsv', sep='\t', index=False)
            create_submission_fasta(fasta, group_metadata, 'Sequence_ID', f'{output_base}.fasta')


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
    create_voc_reports(metadata, args.excluded_vocs, output_dir, batch_name)
    create_sample_status_report(metadata, output_dir, batch_name)
    create_wa_doh_report(metadata, output_dir, batch_name)

    # Only create submissions for sequences that have status "submitted"
    submit_metadata = metadata.loc[metadata['status'] == 'submitted']

    create_gisaid_submission(submit_metadata, args.fasta, output_dir, batch_name, args.gisaid_username)
    create_genbank_submission(submit_metadata, args.fasta, vadr_dir, output_dir, batch_name)
