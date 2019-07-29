"""
Converts a fasta file input to machine-readable, Json format. The Json document
is printed to stdout, so the user will likely want to redirect this to a new
file.
"""
import json
import argparse

def metadata_to_json(urls_r1, urls_r2):
    urls = []

    for url_set in [urls_r1, urls_r2]:
        if url_set:
            urls += url_set.split(" ")

    print(json.dumps( {
        "metadata": {
            "URLs": list(set(urls))
        }
    }, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("urls_r1", type=str, nargs="?",
        metavar = "<urls_r1>")
    parser.add_argument("urls_r2", type=str, nargs="?",
        metavar = "<urls_r2>")

    args = parser.parse_args()
    metadata_to_json(args.urls_r1, args.urls_r2)
