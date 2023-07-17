import argparse
import pandas as pd
import re
from functools import partial

def count_repl(obj, replace_str):
    global feature_count
    feature_count += 1
    return f"{replace_str}_{feature_count}"

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--identifiers", type=str, required=True,
        help = "File path to identifiers TSV")
    parser.add_argument("--vadr-annotation", type=str, required=True,
        help = "File path to VADR feature annotation")
    parser.add_argument("--output", type=str, required=True,
        help = "Output file path")

    args = parser.parse_args()

    identifiers = pd.read_csv(args.identifiers, sep='\t', usecols=['nwgc_id', 'strain_name'])
    identifiers = identifiers[~identifiers.strain_name.isna()]

    # Drop `USA/` prefix and `/YYYY` suffix from strain name to get sequence ID submitted to GenBank
    identifiers['sequence_id'] = identifiers.strain_name.str.replace('USA/', '').str.replace(r'/\d{4}$', '')

    id_dict = dict(zip(identifiers.nwgc_id, identifiers.sequence_id))

    # Replace NWGC IDs with sequence IDs and remove any feature annotations that were not assigned a sequence ID (failed QC)
    new_filedata = ''
    with open(args.vadr_annotation, 'r') as input, open(args.output, 'w') as output:
        filedata = input.read()

        for k,v in id_dict.items():
            # Replace long version IDs with sequence IDs for submission to GenBank
            pattern = r'^>Feature ' + re.escape(str(k)) + r'[|].*$'
            filedata = re.sub(pattern, '>Feature ' + v, filedata, flags=re.M)

            # Replace protein IDs with sequence IDs, including numeric incrementing suffix (e.g. WA-S###_1, WA-S###_2)
            pattern = r'protein_id\t' + re.escape(str(k)) + r'-.*$'
            feature_count = 0
            filedata = re.sub(pattern, partial(count_repl, replace_str='protein_id\t' + v), filedata, flags=re.M)

        include = False
        for line in filedata.split('\n'):
            print(line)
            if line.startswith('>Feature') and line.split(' ')[1] in id_dict.values():
                include = True
            elif line.startswith('>Feature'):
                include = False

            if include:
                new_filedata += line + '\n'

        output.write(new_filedata)
