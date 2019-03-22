"""summary_plots.py
"""
import os
import pandas as pd
import json

def count_reads(samples, references):
    """For each sample x reference combination, count how many reads there
    were total, and how many of those reads mapped to each reference genome.
    """
    out = {}
    for ref in references:
        out[ref] = {}
        for sample in samples:
            out[ref][sample] = {}
            log_file = "summary/bowtie2/{}/{}.log".format(ref, sample)
            if os.path.isfile(log_file):
                with open(log_file, 'r') as f:
                    for line in f:
                        print(line)
                        if "overall alignment rate" in line:
                            percent = float(line.split(" ")[0][:-1])
                    out[ref][sample]["total_percent_aligned"] = percent

    return out


if __name__=="__main__":
    config_path = "config/config.json"

    with open(config_path, 'r') as f:
        config = json.load(f)

    references = config["reference_viruses"].keys()
    samples = ["313913_S11_L001_R1_001", "313920_S18_L001_R1_001", "313922_S19_L001_R1_001", "313926_S23_L001_R1_001"]

    mapped_read_counts = count_reads(samples, references)
    print(mapped_read_counts)
