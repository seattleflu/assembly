"""
POST data to ID3C receiving.consensus_genome.

If POST not successful, will send HTTP Error message to Slack webhook.
"""
import os
import json
import requests
import getpass
import argparse
from requests.exceptions import HTTPError
from urllib.parse import urljoin


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument("--payload", type=str, metavar="<payload.json>",
        help = "Filepath to JSON payload")
    parser.add_argument("--output", type=str, metavar="<output.log>",
        help = "Filepath to dummy output that marks a successful POST")
    parser.add_argument("--http-response", type=str, metavar="<http-response.log>",
        help = "Filepath to http-response log to print HTTP response errors")

    args = parser.parse_args()

    headers = {'Content-type': 'application/json'}

    id3c_url = os.environ['ID3C_URL']
    id3c_username = os.environ['ID3C_USERNAME']
    id3c_password = os.environ['ID3C_PASSWORD']
    id3c_slack_webhook = os.environ['SLACK_WEBHOOK_URL']

    file = open(args.http_response, "w")

    with open(args.payload) as f:
        data = f.read()

    try:
        response = requests.post(
            urljoin(id3c_url, 'v1/receiving/consensus-genome'),
            data=data,
            headers=headers,
            auth=(id3c_username, id3c_password))

        response.raise_for_status()

        if response.ok:
            with open(args.output, "w"):
                pass

    except HTTPError as http_err:
        file.write(str(http_err))

        slack_data = { "text":
            f":rotating_light: Hey {getpass.getuser()}: Assembly failed to upload to ID3C with HTTP status code: " +
            f"{http_err.response.status_code}.\nMore details at `{args.http_response}`"
        }

        try:
            slack_response = requests.post(id3c_slack_webhook,
                data=json.dumps(slack_data), headers=headers)

            slack_response.raise_for_status()

        except HTTPError as slack_http_err:
            file.write(str(slack_http_err))

        raise http_err

    except Exception as err:
        file.write(str(err))

        raise Exception(f"Error: {err} in ID3C POST request.")

    finally:
        file.close()
