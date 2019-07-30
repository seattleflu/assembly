"""
Create a JSON file with for upload to ID3C that contains details of the
sequence read set and reference combination that did not map through bowtie2.
"""
import argparse
import json
import pandas as pd


def import_arguments(metadata_file, reference, nwgc_sfs_map):
    """
    """
    with open(metadata_file) as f:
        metadata = json.load(f)

    reference = reference.split('_')[0]
    sample_id = metadata['metadata']['urls'][0].split('/')[-1].split('_')[0]
    nwgc_sfs_map = pd.read_csv(nwgc_sfs_map, sep="\t", header=None,
                               names=["nwgc_id", "sfs_uuid"],
                               dtype={'nwgc_id': 'str', 'sfs_uuid': 'str'})

    sample_uuid = nwgc_sfs_map \
                    .loc[nwgc_sfs_map['nwgc_id'] == sample_id, 'sfs_uuid']\
                    .item()

    return metadata, reference, sample_uuid


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--metadata", type=str,
        help="The metadata JSON file, output of `rule metadata_to_json`",
        required=True)
    parser.add_argument("--reference", type=str,
        help="The reference that the sample did not map to in bowtie2",
        required=True)
    parser.add_argument("--nwgc-sfs-map", type=str,
        help = "The nwgc-sfs map created by `rule nwgc_sfs_map`",
        required=True)

    args = parser.parse_args()

    metadata, reference, sample_uuid = import_arguments(args.metadata,
                                                        args.reference,
                                                        args.nwgc_sfs_map)

    print(json.dumps({ **metadata,
                       "sample_identifier": sample_uuid,
                       "reference_organism": reference,
                       "status": "notMapped"}, indent=4))

