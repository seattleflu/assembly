"""
Combine consensus genomes with summary statistics and print them as
newline-delimited JSON records.

Records are output to stdout, you will likely want to redirect
stdout to a file.
"""
import argparse
import json
import pandas as pd

def import_nwgc_uuid_map(nwgc_uuid_file):
    """
    Import NWGC IDs and UUIDs as dataframe from a csv file.
    """
    nwgc_uuid = pd.read_excel(nwgc_uuid_file,
                            usecols=["sample", "uuid"],
                            dtype={"sample": "str"}) \
                  .drop_duplicates(subset="sample")
    nwgc_uuid.set_index("sample", inplace=True)
    nwgc_uuid_map = nwgc_uuid.to_dict()["uuid"]
    return nwgc_uuid_map


def combine_summary_and_sets(sequence_read_set_file, summary_file):
    """
    Combine sequence read sets with assembly generated summaries.
    """
    summaries = {}
    with open(summary_file) as f:
        for line in f:
            summary = json.loads(line)
            sample = summary.pop("sample")
            reference_organism = summary.pop("reference").lower()
            if summaries.get(sample):
                if summaries[sample].get(reference_organism):
                    summaries[sample][reference_organism]["summary_stats"] = summary
                else:
                    summaries[sample].update({reference_organism: {
                        "summary_stats": summary
                    }})
            else:
                summaries[sample]={reference_organism: {"summary_stats": summary}}
        f.close()
    with open(sequence_read_set_file) as f:
        for line in f:
            read_set = json.loads(line)
            sample = read_set.pop("sample")
            read_set["urls"] = read_set.pop("urls")
            if summaries.get(sample):
                summaries[sample].update({'metadata': read_set})
            else:
                continue
    return summaries


def convert_nwgc_to_uuid(summary, nwgc_uuid_map):
    """
    Convert the NWGC IDs to SFS UUIDs within *summary* using the
    *nwgc_uuid_map*.
    """
    for sample in summary.fromkeys(summary):
        uuid = nwgc_uuid_map[sample]
        summary[uuid] = summary.pop(sample)
    return summary

def combine_genome_sequences(consensus_file):
    """
    Combine single sequence records into one genome based on the
    sample UUID and the reference organism.
    """
    combined_genomes = {}
    with open(consensus_file) as f:
        for line in f:
            sequence = json.loads(line)
            uuid = sequence.pop("sample_identifier")
            reference_organism = sequence.pop("organism").lower()
            if combined_genomes.get(uuid):
                if combined_genomes[uuid].get(reference_organism):
                    combined_genomes[uuid][reference_organism]["masked_consensus"].append(sequence)
                else:
                    combined_genomes[uuid].update({reference_organism: {"masked_consensus": [sequence]}})
            else:
                combined_genomes.update({uuid: {reference_organism: {"masked_consensus": [sequence]}}})
        f.close()
    return combined_genomes


def print_summary_and_genomes(summary, genome):
    """
    Prints out consensus genomes with summary data and sequence read set urls
    as newline-delimited JSON records.
    """
    for sample in summary:
        for ref in summary[sample]:
            if ref == "metadata":
                continue
            final = {
                "sample_identifier": sample,
                "reference_organism": ref,
                "metadata": summary[sample]["metadata"]
            }
            if genome.get(sample) and genome[sample].get(ref):
                final.update({"status": "complete"})
                # Add summary statistics
                final.update(summary[sample][ref])
                # Add genomic sequences
                final.update(genome[sample][ref])
            else:
                final.update({"status": "notMapped"})
            print(json.dumps(final))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--sequence-read-set", type=str,
        help="ndjson file with sequence read sets",
        required = True)
    parser.add_argument("--summary", type=str,
        help="ndjson file with assembly summary statistics",
        required=True)
    parser.add_argument("--consensus", type=str,
        help="ndjson file with consensus genomes",
        required=True)
    parser.add_argument("--nwgc-uuid", type=str,
        help="Excel file containing map between NWGC ID and SFS UUID",
        required=True)

    args = parser.parse_args()

    nwgc_uuid_map = import_nwgc_uuid_map(args.nwgc_uuid)
    combined_summary = combine_summary_and_sets(args.sequence_read_set, args.summary)
    assembly_summary = convert_nwgc_to_uuid(combined_summary, nwgc_uuid_map)
    combined_genomes = combine_genome_sequences(args.consensus)
    print_summary_and_genomes(assembly_summary, combined_genomes)
