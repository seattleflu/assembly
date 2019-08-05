"""
Creates a config file for Snakemake for one sample/reference pair.
"""
import argparse
import gzip
import json

def check_fastq_sample(fastq_files: list, nwgc_id: str):
    """
    Check all of the *fastq_files* provided contain the *nwgc_id*.
    """
    assert all(nwgc_id in filename for filename in fastq_files), \
        f"NWGC ID {nwgc_id} was not found in all fastq filenames: {fastq_files}"


def sort_fastq_by_lane_and_reads(fastq_files: list):
    """
    Sort FASTQ files by lanes so their reads can be checked in the forward
    and reverse files.
    """
    fastq_by_lane = {}
    fastq_by_read = {
        "R1": [],
        "R2": []
    }
    for filepath in fastq_files:
        filename = filepath.split("/")[-1].split("_")
        lane = filename[-3]
        read = filename[-2]
        if fastq_by_lane.get(lane):
            fastq_by_lane[lane].update({ read: filepath })
        else:
            fastq_by_lane[lane] = { read: filepath }
        fastq_by_read[read].append(filepath)
    return fastq_by_lane, fastq_by_read


def check_r1_r2_reads(fastq_by_lane: dict):
    """
    Check the R1 and R2 reads to make sure they match.
    Prints out the difference if R1 and R2 reads don't match 100%
    """
    for lane in fastq_by_lane:
        print(f"Checking reads in lane: {lane}")
        r1_file = fastq_by_lane[lane]['R1']
        r2_file = fastq_by_lane[lane]['R2']
        r1_read_ids = find_read_ids(r1_file)
        r2_read_ids = find_read_ids(r2_file)
        different_read_ids = r1_read_ids.symmetric_difference(r2_read_ids)
        percent_different = len(different_read_ids) / len(r1_read_ids) * 100
        assert r1_read_ids == r2_read_ids, \
            f"R1 and R2 reads differ by {percent_different}. Different read ids: {different_read_ids}"
        print("100% of R1 and R2 reads are paired")


def find_read_ids(fastq_file: str) -> set:
    """
    Finds all read IDs in FASTQ file and returns them as a set.
    """
    read_ids = set()
    with gzip.open(fastq_file, "rt") as r1_file:
        for line in r1_file:
            if line.startswith("@"):
                id = line.split()[0]
                read_ids.add(id)
    return read_ids


def define_references(target: str, target_ref_map: str):
    """
    Uses the *target_ref_map* to determine which references map to the provided
    *target*
    """
    with open(target_ref_map, 'r') as f:
        reference_map = json.load(f)
        f.close()
    return reference_map[target]


def create_barcode_match(nwgc_id: str, sfs_uuid: str):
    """
    Create barcode match txt file needed for the fasta header renaming.
    """
    barcode_match_file = f"data/{sfs_uuid}_key_value.txt"
    with open(barcode_match_file, "w+") as f:
        f.write(nwgc_id + "\t" + sfs_uuid)
        f.close()
    return barcode_match_file


def create_config_file(config_template: str,
                       config_file: str,
                       fastq_by_reads: dict,
                       nwgc_id: str,
                       references: list,
                       barcode_match_file: str):
    """
    Create a new config file for the sample/reference pair based on the
    provided *config_template*
    """
    with open(config_template, "r") as template:
        config = json.load(template)
        template.close()

    config["fastq_files"] = fastq_by_reads
    config["references"] = references
    config["sample"] = nwgc_id
    config["barcode_match"] = barcode_match_file

    with open(config_file, "w") as new_config:
        new_config.seek(0)
        new_config.truncate()
        json.dump(config, new_config, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--config-template", type=str,
        metavar = "<config_template.json>",
        help = "File path to config file template",
        required = True)
    parser.add_argument("--config-file", type=str,
        metavar = "<config.json>",
        help = "File path to config file to be used for snakemake",
        required = True)
    parser.add_argument("--fastq-file", type=str, nargs="+",
        metavar = "<sample.fastq.gz>",
        help = "File path to fastq file. Expected format of fastq filename: '346757_10DB_288_S281_L001_R1_001.fastq.gz",
        required = True)
    parser.add_argument("--nwgc-id", type=str, nargs="?",
        help = "NWGC ID for the sample that matches ID in fastq filename",
        required = True)
    parser.add_argument("--sfs-uuid", type=str, nargs="?",
        help = "The SFS UUID of the sample.",
        required = True)
    parser.add_argument("--target", type=str, nargs="?",
        help = "The target to be used as reference in the assembly pipeline",
        required = True)
    parser.add_argument("--target-ref-map", type=str, nargs="?",
        help = "The JSON file that contains the target/reference map",
        default = "references/target_reference_map.json")

    args = parser.parse_args()

    # Check all fastq files are for the expected sample
    check_fastq_sample(args.fastq_file, args.nwgc_id)

    # Sort FASTQ filename by lane and by read
    fastq_by_lane, fastq_by_reads = sort_fastq_by_lane_and_reads(args.fastq_file)

    # Check all reads within one lane match in R1 and R2
    check_r1_r2_reads(fastq_by_lane)

    references = define_references(args.target, args.target_ref_map)
    barcode_match = create_barcode_match(args.nwgc_id, args.sfs_uuid)
    create_config_file(args.config_template, args.config_file,
                       fastq_by_reads, args.nwgc_id, references, barcode_match)
