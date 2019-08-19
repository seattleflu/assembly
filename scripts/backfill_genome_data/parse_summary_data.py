"""
Parse through Bowtie2 and BAMstats summaries given their directory paths and
print them as newline-delimited JSON records.

Summary records are output to stdout, you will likely want to redirect
stdout to a file.
"""
import glob
import json
import argparse

def get_reads_and_align_rate(directory):
    """
    Parse out number of reads and align rate from bowtie2 summary log.
    """
    summary = {}
    filepath = directory + "/*/*.log"
    for filename in glob.glob(filepath):
        sample, reference = get_sample_reference_pair(filename)
        default = {"align_rate": None, "reads": None}
        if summary.get(sample):
            summary[sample][reference] = default
        else:
            summary[sample] = { reference: default }
        with open(filename) as f:
            lines = f.readlines()
            if len(lines) > 1:
                number_of_reads = lines[1].split()[0]
                align_rate = lines[-1].split()[0]
                summary[sample][reference]["align_rate"] = align_rate
                summary[sample][reference]["reads"] = number_of_reads
            f.close()

    return summary


def get_sample_reference_pair(filename):
    """
    Parse out sample id and reference name from *filename*
    """
    sample = filename.split('/')[-1].split('.')[0]
    reference = filename.split('/')[-2].split('_')[0]
    return sample, reference


def get_coverage_summary(directory, bowtie_summary):
    """
    Parse out mean coverage depth for each segment of the genome and add it
    to the *bowtie_summary*.
    """
    summary = bowtie_summary if bowtie_summary else {}
    filepath = directory + "/*/*.coverage_stats.txt"
    for filename in glob.glob(filepath):
        sample, reference = get_sample_reference_pair(filename)
        default = {"mean_coverage_depths": None}
        if summary.get(sample):
            summary[sample][reference].update(default)
        else:
            summary[sample] = { reference: default }
        with open(filename) as f:
            lines = f.readlines()[1:]
            segments = list(map(lambda x: x.split()[0].split('|')[-1],
                                lines))
            mean_coverage_depths = list(map(lambda x: float(x.split()[2]),
                                            lines))
            coverage_summary = dict(zip(segments, mean_coverage_depths))
            summary[sample][reference]['mean_coverage_depths'] = coverage_summary
    return summary


def print_summary_as_ndjson(summary):
    """
    Print out summary in newline-delimited format
    """
    for sample in summary:
        sample_id = sample
        for target in summary[sample]:
            reference = target
            summary_stats = summary[sample][reference]
            summary_stats['sample'] = sample_id
            summary_stats['reference'] = reference
            print(json.dumps(summary_stats))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__ ,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--bowtie2", type=str,
        help="Directory of bowtie2 summaries",
        required=True)
    parser.add_argument("--bamstats", type=str,
        help="Directory of BAMstats summaries",
        required=True)

    args = parser.parse_args()

    bowtie_summary = get_reads_and_align_rate(args.bowtie2)
    bowtie_and_bamstats_summary = get_coverage_summary(args.bamstats, bowtie_summary)
    print_summary_as_ndjson(bowtie_and_bamstats_summary)
