import os
import sys 
import argparse
import json
import pysam
import pandas as pd

ROOT = os.path.basename(__file__)
sys.path.insert(0, ROOT)
from EpiCancer.module.methyl_qpcr_designer import BEDqPCRdesigner

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--bed', required=True)
    parser.add_argument('--bed_id_column', required=False, default='id')
    parser.add_argument('--output', required=True)
    parser.add_argument('--config', required=True)
    
    return parser

def process_arguments(parser: argparse.ArgumentParser):
    args = parser.parse_args()
    args.fasta = os.path.abspath(args.fasta)
    args.bed = os.path.abspath(args.bed)
    args.output = os.path.abspath(args.output)
    args.config = os.path.abspath(args.config)
    
    return args

def main():
    parser = parse_arguments()
    args = process_arguments(parser)

    bed_qpcr_designer = BEDqPCRdesigner(reference_fasta_path=args.fasta, bed_path=args.bed, config_path=args.config, bed_id_column=args.bed_id_column)
    bed_qpcr_designer.design_primers()
    bed_qpcr_designer.set_amplicon_df()
    bed_qpcr_designer.save_primers(args.output)

if __name__ == '__main__':
    main()

