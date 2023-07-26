import argparse
from Bio import SeqIO
import os
import logging
import pandas as pd
from pathlib import Path
from collections import Counter

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

        # insert sample id as first column and append to metrics dataframe
        df.insert(0, 'SampleId', identifier)
        metrics_df = metrics_df.append(df[PICARD_METRICS_COLS + ['SampleId']], ignore_index=True)
    return metrics_df


def append_metrics_counts(metrics_df:pd.DataFrame, fasta: str) -> pd.DataFrame:
    try:
        LOG.debug(f"Parsing FASTA file: {fasta}")
        sequences = list(SeqIO.parse(fasta, 'fasta'))
    except Exception as err:
        LOG.error(f"There was a problem parsing sequences from {fasta}: \n{err}")

    sequence_dict = {}

    # Calculate additional metrics from sequence and add to metrics dataframe
    for record in sequences:
        identifier = str(record.id.split('|')[0]).lower()
        sequence = record.seq.upper()

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
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, 'CONTIG_NAME'] = identifier

        for k,v in counts.items():
            metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, k] = str(v)

        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, 'PCT_N_MAPPED'] = (counts['COUNT_N'] / counts['CONSENSUS_FASTA_LENGTH']) * 100
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, 'REFERENCE_LENGTH'] = metrics_df['GENOME_TERRITORY'].astype(int)
        metrics_df.loc[metrics_df['SampleId'].str.lower() == identifier, 'PCT_N_REFERENCE'] = (counts['COUNT_N'] / metrics_df['GENOME_TERRITORY']) * 100

    return metrics_df.convert_dtypes()


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
