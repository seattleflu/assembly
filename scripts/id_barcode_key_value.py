"""
Generate NWGC sample ID and SFS barcode key-value pairs
"""
import os, sys
import argparse
import pandas as pd 

def create_key_value_pairs(filename: str, NWGC: str, SFS:str, output:str):
    """
    Find NSGC sample ID and SFS barcode within given Excel file
    and print out to a txt file as key-value pairs separated by tab
    """
    df = pd.read_excel(filename, usecols=[NWGC, SFS])
    df = df[[NWGC, SFS]]
    df.to_csv(output, header=None, index=None, sep="\t")

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Create key/value pairs of NWGC ID and SFS ID",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("match-file", type=str, nargs="?",
        metavar="NWGC-SFS-Match",
        help="File path to Excel file that matches NWGC ID to SFS ID")
    parser.add_argument("NWGC", type=str, nargs="?",
        metavar="NWGC-column",
        help="Column name of column containing NWGC sample IDs")
    parser.add_argument("SFS", type=str, nargs="?",
        metavar="SFS-column",
        help="Column name of column containing SFS UUIDs")
    parser.add_argument("output-file", type=str, nargs="?",
        metavar="Output",
        help="File path for output key/value file")
    
    try:
        args = parser.parse_args()
        create_key_value_pairs(args.match_file, 
                               args.NWGC, 
                               args.SFS, 
                               args.output_file)
    except:
        parser.print_help()
        sys.exit(0)
