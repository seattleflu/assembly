"""
Converts a fasta file input to machine-readable, Json format. The Json document
is printed to stdout, so the user will likely want to redirect this to a new
file.
"""
import json
import argparse
from urllib.parse import urljoin

def metadata_to_json(urls, url_prefix):
    print(json.dumps( {
        "metadata": {
            "urls": list(set(urljoin(url_prefix, url) for url in urls))
        }
    }, indent=4))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("urls", type=str, nargs="+", metavar = "<urls>")
    parser.add_argument("--url-prefix",
        default = "file://rhino.fhcrc.org",
        help = "Base for fully-qualifying sequence read set URLs")

    args = parser.parse_args()
    metadata_to_json(args.urls, args.url_prefix)
