"""
Combine the masked consensus Json fasta file and the bamstats/bowtie2 summary
stats Json file into one Json for upload to ID3C.
"""
import re
import json
import argparse

def check_sample_ids(masked_consensus: str, summary_stats: str) -> None:
    """"""
    regex = '.*/([0-9]{6,}).*.json'
    match = re.search(regex, summary_stats)
    sample_id = match.group(1)

    assert re.search(regex, masked_consensus).group(1), \
        "The NWGC sample IDs for the given files do not match."


def check_organism_ids(masked_consensus: dict, summary_stats: dict) -> None:
    """"""
    assert masked_consensus['reference_organism'].lower() == \
        summary_stats['reference_organism'].lower(), \
        "The reference organisms for the given files do not match!"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__ ,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--masked-consensus", type=str,
        help="The masked consensus Json file, output of `rule masked_consensus_to_json`",
        required=True)
    parser.add_argument("--summary-stats", type=str,
        help="The summary stats Json file, output of `rule summary_stats_to_json`",
        required=True)
    parser.add_argument("--metadata", type=str,
        help="The metadata Json file, output of `rule metadata_to_json`",
        required=True)
    parser.add_argument("--status", type=str,
        choices=['complete', 'failed', 'notMapped'],
        help="The mapping status of the consensus genome.",
        required=True)

    args = parser.parse_args()
    check_sample_ids(args.masked_consensus, args.summary_stats)

    with open(args.summary_stats) as f:
        summary_stats = json.load(f)

    with open(args.masked_consensus) as f:
        masked_consensus = json.load(f)

    with open(args.metadata) as f:
        metadata = json.load(f)

    status = { 'status': args.status }
    check_organism_ids(masked_consensus, summary_stats)

    print(json.dumps({
        **status,
        **masked_consensus,
        **summary_stats,
        **metadata }, indent=4))
