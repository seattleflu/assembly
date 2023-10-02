import argparse
import pandas as pd
import re
from Bio import SeqIO

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = __doc__,
        formatter_class = argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--fasta", type=str, required=True,
        help = "File path to identifiers TSV")
    parser.add_argument("--biosample-accessions", type=str, required=True,
        help = "File path to VADR feature annotation")
    parser.add_argument("--output", type=str, required=True,
        help = "Output file path")
    parser.add_argument("--pathogen", type=str, required=True,
        choices=['rsv-a', 'rsv-b'])


    PATHOGEN_ORGANISM_MAP = {
        'rsv-a': 'Human respiratory syncytial virus A',
        'rsv-b': 'Human respiratory syncytial virus B'
    }
    args = parser.parse_args()
    accessions = pd.read_csv(args.biosample_accessions, sep='\t', usecols=['accession', 'bioproject_accession', 'sample_name'])

    # Drop `USA/` prefix and `/YYYY` suffix from strain name to get sequence ID submitted to GenBank
    accessions['sequence_id'] = accessions.sample_name.str.replace('USA/', '').str.replace(r'/\d{4}$', '')

    accession_dict = accessions.groupby('sequence_id')[['sequence_id','accession','bioproject_accession']].apply(lambda x: x.set_index('sequence_id').to_dict(orient='index')).to_dict()

    with open(args.fasta) as original, open(args.output, 'w') as output:
        sequences = SeqIO.parse(original, 'fasta')
        for s in sequences:
            entry = accession_dict[s.id][s.id]
            s.description = s.id + f" [BioProject={entry['bioproject_accession']}] [BioSample={entry['accession']}] [organism={PATHOGEN_ORGANISM_MAP[args.pathogen]}]"
            SeqIO.write(s, output, 'fasta-2line')
