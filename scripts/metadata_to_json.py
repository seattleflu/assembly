"""
Converts a fasta file input to machine-readable, Json format. The Json document
is printed to stdout, so the user will likely want to redirect this to a new
file.
"""
import json
import argparse

def metadata_to_json(urls):
    print(json.dumps( {
        "metadata": {
            "urls": list(set(urls))
        }
    }, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("urls", type=str, nargs="+", metavar = "<urls>")

    args = parser.parse_args()
    metadata_to_json(args.urls)
