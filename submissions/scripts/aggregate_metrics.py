import argparse
from Bio import SeqIO
import os
import logging
import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter
import sys

LOG_LEVEL = os.environ.get("LOG_LEVEL", "debug").upper()

logging.basicConfig(
    level = logging.ERROR,
    format = "[%(asctime)s] %(levelname)-8s %(message)s",
    datefmt = "%Y-%m-%d %H:%M:%S%z")

logging.captureWarnings(True)
LOG = logging.getLogger(__name__)
LOG.setLevel(LOG_LEVEL)

def combine_metrics(metrics_files: list) -> pd.DataFrame:
    LOG.debug(f"Aggregating {len(metrics_files)} metrics files.")

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

    # flu_segment_names = [
    #     'HA',
    #     'NA',
    #     'NP',
    #     'PA',
    #     'PB1',
    #     'PB2',
    #     'MP',
    #     'NS',
    # ]

    metrics_df = pd.DataFrame()

    for metrics_file in metrics_files:
        # parse sample id from filename
        identifier = Path(metrics_file).stem.split('.')[0]

        # metrics TSV for each sample should be on lines 7 (headers) and 8 (data) of individual metrics files output by Picard
        df = pd.read_csv(metrics_file, sep='\t', skiprows=6, nrows=1)

        if not all(c in df.columns for c in PICARD_METRICS_COLS):
            raise Exception(f"Metrics file {metrics_file} does not match expected columns: {', '.join(PICARD_METRICS_COLS)}")
        elif df.empty:
            raise Exception(f"Metrics file {metrics_file} does not contain expected TSV data.")

        # for flu, create a metrics entry for each segment
        # (Note that the Picard metrics are calculated per genome, not per segment, but the whole-genome metrics will be applied to each segment's entry, which is inaccurate but the best I can do right now)
        # if pathogen.startswith('flu'):
        #     # parse sample id from filename
        #     sampleid = Path(metrics_file).stem.split('.')[0]
        #     for segment_name in flu_segment_names:
        #         # make df fresh for each segment; probably a more efficient way to do this...
        #         df = pd.read_csv(metrics_file, sep='\t', skiprows=6, nrows=1)
        #         # identifier is sample id plus segment name
        #         identifier = str(sampleid) + "|" + segment_name
        #         # insert sample id as first column and append to metrics dataframe
        #         df.insert(0, 'SampleId', identifier)
        #         metrics_df = metrics_df.append(df[PICARD_METRICS_COLS + ['SampleId']], ignore_index=True)

        # parse sample id from filename
        identifier = Path(metrics_file).stem.split('.')[0]            
        # insert sample id as first column and append to metrics dataframe
        df.insert(0, 'SampleId', identifier)
        metrics_df = metrics_df.append(df[PICARD_METRICS_COLS + ['SampleId']], ignore_index=True)
    return metrics_df


def append_metrics_counts(metrics_df:pd.DataFrame, fasta: str, fasta_short_ids: str) -> pd.DataFrame:
    try:
        LOG.debug(f"Parsing FASTA file: {fasta}")
        sequences = list(SeqIO.parse(fasta, 'fasta'))
    except Exception as err:
        LOG.error(f"There was a problem parsing sequences from {fasta}: \n{err}")

    sequence_dict = {}

    original_record_ids = {}

    # Calculate additional metrics from sequence and add to metrics dataframe
    for record in sequences:
        identifier = str(record.id.split('|')[0]).lower()
        sequence = record.seq.upper()

        # shorten identifiers to use with GenBank trimming script later
        original_record_ids[identifier] = record.id
        record.description = identifier
        record.id = identifier
        SeqIO.write(sequences, fasta_short_ids, 'fasta-2line')

        counts = {
            'COUNT_A': sequence.count("A"),
            'COUNT_C': sequence.count("C"),
            'COUNT_G': sequence.count("G"),
            'COUNT_T': sequence.count("T"),
            'COUNT_N': sequence.count("N"),
            'CONSENSUS_FASTA_LENGTH': len(sequence)
        }

        merged = dict(Counter(counts) + Counter(sequence_dict.get(identifier, {})))
        sequence_dict[identifier] = merged

    for identifier, counts in sequence_dict.items():
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier.lower(), 'CONTIG_NAME'] = identifier

        for k,v in counts.items():
            metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier.lower(), k] = str(v)

        # temporarily allow no-N sequences through so we can look at them further
        if 'COUNT_N' not in counts:
            print('Warning: ' + str(identifier) + ' does not contain Ns. Processing anyway.')
            counts['COUNT_N'] = 0

        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier.lower(), 'PCT_N_MAPPED'] = (counts['COUNT_N'] / counts['CONSENSUS_FASTA_LENGTH']) * 100
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier.lower(), 'REFERENCE_LENGTH'] = metrics_df['GENOME_TERRITORY'].astype(int)
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier.lower(), 'PCT_N_REFERENCE'] = (counts['COUNT_N'] / metrics_df['GENOME_TERRITORY']) * 100

    return metrics_df.convert_dtypes(), original_record_ids


def append_metrics_counts_flu(metrics_df:pd.DataFrame, fasta: str, fasta_short_ids: str) -> pd.DataFrame:
    """
    Appends whole genome metrics counts to metrics_df, and creates segments_metrics_df with segment metrics counts
    """
    try:
        LOG.debug(f"Parsing FASTA file: {fasta}")
        sequences = list(SeqIO.parse(fasta, 'fasta'))
    except Exception as err:
        LOG.error(f"There was a problem parsing sequences from {fasta}: \n{err}")

    sequence_dict = {}
    # keeps track of segment metrics

    original_record_ids = {}

    # Calculate additional metrics from sequence and add to metrics dataframe
    # for flu, we want to calculate all additional metrics both per genome and per segment
        # genome-wide flu metrics will be appended to metrics_df, and segment metrics will be put in a new df, segment_metrics_df
    for record in sequences:
        split_identifier = record.id.split('|')
        # for flu, create two identifiers: a segment identifier which includes the sample name and segment name,
        # and an identifier which is just sample name (nwgc id)
        # the segment identifier will be used for fasta_short_ids and segment_metrics_df
        # identifier will be used for metrics_df
        # for segment identifier, use hyphen instead of "|" when writing to fasta_short_ids due to GenBank trimming script limitations
        segment_identifier = str(split_identifier[0]).lower() + "-" + str(split_identifier[-1])
        #identifier = str(record.id.split('|')[0]).lower()
        sequence = record.seq.upper()

        # write out shortened segment ids to a fasta to use with GenBank trimming script later
        original_record_ids[segment_identifier] = record.id
        record.description = segment_identifier
        record.id = segment_identifier
        SeqIO.write(sequences, fasta_short_ids, 'fasta-2line')

        # after writing to fasta_short_ids, replace "-" with "|" for compatibility with other functions/scripts
        segment_identifier = str(split_identifier[0]).lower() + "|" + str(split_identifier[-1])

        counts = {
            'COUNT_A': sequence.count("A"),
            'COUNT_C': sequence.count("C"),
            'COUNT_G': sequence.count("G"),
            'COUNT_T': sequence.count("T"),
            'COUNT_N': sequence.count("N"),
            'CONSENSUS_FASTA_LENGTH': len(sequence)
        }

        merged = dict(Counter(counts) + Counter(sequence_dict.get(segment_identifier, {})))
        sequence_dict[segment_identifier] = merged

    # for each SampleId in metrics_df, create 8 segment entries in segment_metrics_df
    flu_segment_names = [
        'HA',
        'NA',
        'NP',
        'PA',
        'PB1',
        'PB2',
        'MP',
        'NS',
    ]

    # create list of segment ids
    segment_id_list = []
    for sampleid in metrics_df['SampleId']:
        for segment_name in flu_segment_names:
            segment_id_list.append(str(sampleid)+'|'+segment_name)

    # create segment_metrics_df
    segment_metrics_df = pd.DataFrame(data=segment_id_list,columns=['SegmentId'])
    segment_metrics_df['SampleId'] = segment_metrics_df['SegmentId'].apply(lambda x: x.split('|')[0])
    segment_metrics_df['SegmentName'] = segment_metrics_df['SegmentId'].apply(lambda x: x.split('|')[-1])

    for segment_identifier, counts in sequence_dict.items():
        # add to segment metrics df
        # contig name (note, this is redundant with SegmentId, keeping for now)
        segment_metrics_df.loc[segment_metrics_df['SegmentId'] == segment_identifier, 'CONTIG_NAME'] = segment_identifier
        # counts
        for k,v in counts.items():
            segment_metrics_df.loc[segment_metrics_df['SegmentId'] == segment_identifier, k] = str(v)
        # temporarily allow no-N sequences through so we can look at them further
        if 'COUNT_N' not in counts:
            print('Warning: ' + str(segment_identifier) + ' does not contain Ns. Processing anyway.')
            counts['COUNT_N'] = 0   
        # %N mapped (out of sequence length, how many Ns)
        segment_metrics_df.loc[segment_metrics_df['SegmentId'] == segment_identifier, 'PCT_N_MAPPED'] = (counts['COUNT_N'] / counts['CONSENSUS_FASTA_LENGTH']) * 100
        # can't calculate %N reference for segments because Picard doesn't give us reference length per segment

    # from segment_metrics_df, calculate additional whole genome metrics and append them to metrics_df
    # count_a, count_c, count_g, count_t, count_n, and consensus_fasta_length are the sums of the corresponding metrics sums for that sample
    for metric in counts:
        grouped_metric = segment_metrics_df.groupby('SampleId')[metric].apply(lambda x: sum(x.astype(int)) if x.notnull().all() else pd.NA)
        for sample in grouped_metric.index:
            metrics_df.loc[metrics_df['SampleId'] == sample, metric] = grouped_metric[sample]
        # TODO convert to whole numbers

    # TODO contig_name doesn't make sense in the context of whole genome metrics for flu, since the whole genome is actually 8 contigs
        # but for now, I am putting the identifier (nwgc_id) as the contig name for the whole genome metrics df
    metrics_df['CONTIG_NAME'] = metrics_df['SampleId']

    # calculate %N_mapped, reference length, and %N_ref
    metrics_df['PCT_N_MAPPED'] = (metrics_df['COUNT_N'] / metrics_df['CONSENSUS_FASTA_LENGTH']) * 100
    metrics_df['REFERENCE_LENGTH'] = metrics_df['GENOME_TERRITORY']
    metrics_df['PCT_N_REFERENCE'] = (metrics_df['COUNT_N'] / metrics_df['REFERENCE_LENGTH']) * 100

    return metrics_df.convert_dtypes(), original_record_ids, segment_metrics_df.convert_dtypes()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--fasta', required=True,
        help='The fasta file to process.')
    parser.add_argument('--metrics', required=True,
        help='The Picard metrics files to process.', nargs='+')
    parser.add_argument('--output', required=True,
        help='The output metrics TSV.')
    args = parser.parse_args()

    metrics_df = combine_metrics(args.metrics)
    append_metrics_counts(metrics_df, args.fasta).to_csv(args.output, sep='\t', index=False)

    LOG.debug(f"Saved {len(metrics_df)} records to: {args.output}")
