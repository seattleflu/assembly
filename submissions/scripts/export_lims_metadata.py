import base64
from datetime import datetime, timezone
import hmac
import hashlib
from http.client import HTTPException
from requests import Request, Session
import os
import pandas as pd
import json
import argparse
import sys

HMAC_KEY_ID = os.environ.get('LIMS_API_KEY_ID')
HMAC_SECRET = os.environ.get('LIMS_API_KEY_SECRET')

LIMS_SERVER = 'https://api.lims.bbi.gs.washington.edu'
LIMS_API_PATH = '/api/v1/sfs-specimens/sequencing-submission-metadata'

def get_lims_sequencing_metadata(barcodes: list):
    body = json.dumps([{"anyBarcode": b} for b in barcodes if not str(b).endswith("_exp")])

    verb = 'POST'

    # Prepare the HTTP request, with all its headers except for Authorization.
    request = Request(
        verb,
        f'{LIMS_SERVER}{LIMS_API_PATH}',
        data=body,
        headers={
            'hmac-key-id': HMAC_KEY_ID,
            'content-type': 'application/json'
        }
    )
    prepared_request = request.prepare()

    # The nonce is a value that prevents replay attacks. Currently, the server expects it to be a UNIX timestamp, i.e. the
    # number of whole milliseconds since the epoch in the UTC time zone.
    unix_timestamp_seconds_utc = int(datetime.now(timezone.utc).replace(tzinfo=timezone.utc).timestamp() * 1000)
    nonce = f'{unix_timestamp_seconds_utc}'

    # The secret provided by the LIMS application is base64-encoded, with the final "==" stripped off.
    hmac_signer = hmac.new(base64.b64decode(HMAC_SECRET + '=='), digestmod=hashlib.sha512)

    # Append the nonce, the HTTP verb, the URL (just path and query string components), and a hash digest of the body if
    # there is one to the data to be HMAC-signed.
    hmac_signer.update(nonce.encode())
    hmac_signer.update(verb.encode())
    hmac_signer.update(LIMS_API_PATH.encode())
    if body is not None:
       body_hash = hashlib.md5(prepared_request.body.encode()).hexdigest().encode()
       hmac_signer.update(body_hash)

    # Generate the signature.
    signature = hmac_signer.hexdigest()

    # Add the Authorization header.
    prepared_request.headers['Authorization'] = f'HMAC {nonce}:{signature}'

    # Make the HTTP request.
    session = Session()
    response = session.send(prepared_request)

    if response.status_code == 200 and response.content is not None:
        df = pd.json_normalize(json.loads(response.content), max_level=3)

        column_map = {
            'specimen.ids.sampleIdentifierForDohReporting': 'sfs_identifier_for_doh_reporting',
            'specimen.ids.sampleId': 'sfs_sample_identifier',
            'specimen.ids.sampleBarcode': 'sfs_sample_barcode',
            'specimen.ids.collectionBarcode': 'sfs_collection_barcode',
            'specimen.origin.collection.collectionDate': 'collection_date',
            'specimen.origin.collection.swabType': 'swab_type',
            'specimen.origin.source': 'source',
            'subject.pumaGeoId':'puma',
            'subject.county':'county',
            'subject.state': 'state',
            'sequencing.sentinelSurveillance': 'baseline_surveillance',
            'error.message': 'error_message',
        }
        df.rename(columns=column_map, inplace=True)

        df['baseline_surveillance'] = df['baseline_surveillance'].map({True:'t', False:'f'})

        return df
    else:
        print(f'Error:\n HTTP status {response.status_code}\n {response.content}')
        response.raise_for_status()


def add_lims_metadata(input_df:pd.DataFrame):

    barcodes = df['sfs_sample_barcode'].to_list()
    result = get_lims_sequencing_metadata(barcodes)
    result.to_csv('~/downloads/lims-metadata-result.csv')
    print(f"Received {len(result)} records from LIMS.")

    output_df = pd.merge(input_df, result, on=['sfs_sample_barcode'])
    remaining_records = input_df[(~input_df.nwgc_id.isin(output_df.nwgc_id))]
    remaining_records_joined = pd.merge(remaining_records, result, left_on='sfs_sample_barcode', right_on='sfs_collection_barcode')
    remaining_records_joined.to_csv('~/downloads/remaining_records_joined')

    output_df = pd.concat([output_df, remaining_records_joined], ignore_index=True).sort_values(by=['nwgc_id'])
    #remaining_records, result, left_on='sfs_sample_barcode', right_on='sfs_collection_barcode'))

    return output_df[['nwgc_id','sfs_sample_identifier','sfs_sample_barcode','sfs_collection_barcode','sfs_identifier_for_doh_reporting','collection_date','swab_type','source','puma','county','state','baseline_surveillance','error_message']]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--sfs-sample-barcodes", required=True,
        help = "File path to the sfs-sample-barcdes CSV file")

    parser.add_argument("--output-file",
        help = "Output CSV file")

    args = parser.parse_args()

    df = pd.read_csv(args.sfs_sample_barcodes)

    if args.output_file:
        add_lims_metadata(df).to_csv(args.output_file, index=False)
    else:
        print(add_lims_metadata(df))
