"""
De-identifies a fasta by replacing NWGC IDs with 'seq-1', 'seq-2', etc.,
and prints out a dictionary for re-identification.
"""

import argparse
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-i", "--input-fasta", type=str, required=True,
        help = "Path to input fasta")
    parser.add_argument("-o", "--output-fasta", type=str, required=True,
        help = "Path to output fasta")

    args = parser.parse_args()

    with open(args.input_fasta, 'r') as original, open(args.output_fasta, 'w') as output:
        counter = 0
        id_dict = {}
        for record in SeqIO.parse(original, 'fasta'):
            nwgc_id = record.id.split('|')[0]
            segment = record.id.split('|')[-1]
            if nwgc_id not in id_dict:
                counter +=1
                id_dict[nwgc_id] = counter
            record.id = 'seq-'+str(counter)+'|'+segment
            record.description = ''
            SeqIO.write(record, output, 'fasta-2line')

    print(id_dict)

exit


