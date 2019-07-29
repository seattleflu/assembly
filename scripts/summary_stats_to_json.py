"""
Parse through Bowtie2 and BAMstats summaries given their filenames.
"""
import json
import argparse

def get_reads_and_align_rate(filename):
    """
    Parse out number of reads and align rate from bowtie2 summary log.
    """
    summary = {}
    with open(filename) as f:
        lines = f.readlines()
        number_of_reads = lines[1].split()[0]
        align_rate = lines[-1].split()[0]
        f.close()
        summary["align_rate"] = align_rate
        summary["reads"] = number_of_reads

    return summary


def get_coverage_summary(filename):
    """
    Parse out mean coverage depth for each segment of the genome and add it
    to the *bowtie_summary*.
    """
    summary = {}
    with open(filename) as f:
        lines = f.readlines()[1:]
        segments = list(map(lambda x: x.split()[0].split('|')[-1], lines))
        mean_coverage_depths = list(map(lambda x: float(x.split()[2]), lines))
        coverage_summary = dict(zip(segments, mean_coverage_depths))
        summary['mean_coverage_depths'] = coverage_summary

    return summary


def get_sample_reference_pair(filename):
    """
    Parse out sample id and reference name from *filename*
    """
    sample = filename.split('/')[-1].split('.')[0]
    reference = filename.split('/')[-2].split('_')[0]
    return sample, reference


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__ ,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bowtie2", type=str, nargs="?",
        help="bowtie2 summary file",
        required=True)
    parser.add_argument("--bamstats", type=str, nargs="?",
        help="BAMstats summary file",
        required=True)

    args = parser.parse_args()

    sample_id, reference_id = get_sample_reference_pair(args.bowtie2)
    assert (sample_id, reference_id) == (get_sample_reference_pair(args.bamstats)), \
        "The sample-reference pairs in the provided filenames do not match."

    bowtie2 = get_reads_and_align_rate(args.bowtie2)
    bamstats = get_coverage_summary(args.bamstats)

    print(json.dumps({
        'reference_organism': reference_id,
        'summary_stats': { **bowtie2, **bamstats }
    }, indent=4))
