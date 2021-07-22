"""
Match PUMA to county in CSV file with a 'puma' column.
"""
import argparse
import sys
import pandas as pd
from pathlib import Path

base_dir = Path(__file__).resolve().parent.parent.parent

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--metadata", type=str, required=True,
        metavar="<metadata.csv>",
        help="File path to metadata CSV")
    parser.add_argument("--puma-to-county", type=str,
        metavar="<puma-to-county.csv>",
        default=f"{base_dir}/submissions/source-data/washington_pumas_to_county.csv",
        help="File path to CSV with PUMAs to county conversion")
    parser.add_argument("--output-metadata", type=str, required=True,
        metavar="<output-metadata.csv>",
        help="File path to output metadata CSV")
    parser.add_argument("--output-missing", type=str,
        metavar="<missing-county.csv>",
        default=sys.stdout,
        help="File path to output records missing county. "
             "Prings to stdout by default.")

    args = parser.parse_args()

    pumas_to_county = pd.read_csv(args.puma_to_county, dtype='string')
    assert 'puma' in pumas_to_county.columns, \
        "<puma_to_county.csv> is missing the required `puma` column"
    assert 'county' in pumas_to_county.columns, \
        "<puma_to_county.csv> is missing the required `county` column"

    metadata = pd.read_csv(args.metadata, dtype='string')
    assert 'puma' in metadata.columns, \
        "<metadata.csv> is missing the required `puma` column"

    metadata = metadata.merge(pumas_to_county[['puma', 'county']], on='puma', how='left')

    missing_county = metadata.loc[metadata.county.isnull()]
    if not missing_county.empty:
        print(f"{len(missing_county.index)} records missing county.")
        missing_county.to_csv(args.output_missing, index=False)

    metadata.to_csv(args.output_metadata, index=False)
