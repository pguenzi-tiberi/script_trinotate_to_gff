#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
import pandas as pd
from sympy import im 
from lib import function
from lib import prepare
from lib import separiting_file


def run() :
    parser = argparse.ArgumentParser(
        prog="gff_builder_trinotate",
        description="\n\n A little program to make a gff from trinotate report and other files ",
        formatter_class=argparse.RawTextHelpFormatter,
        add_help=True,
    )

    mandatory_args = parser.add_argument_group("Mandatory arguments")

    # Mandatory arguments

    # Trinotate report
    
    mandatory_args.add_argument(
        "--trinotate_report",
        action="store",
        dest="trinotate_report",
        help="path of trinotate report file",
        default="",
        required=True,
    )
    '''
    mandatory_args.add_argument(
        "--targetp_report",
        action="store",
        dest="targetp_report",
        help="path of targetp file",
        default="",
        required=True,
    )
    '''
    mandatory_args.add_argument(
        "--gff_input",
        action="store",
        dest="gff_input",
        help="path of gff input file",
        default="",
        required=True,
    )

    mandatory_args.add_argument(
        "--id_genome",
        action="store",
        dest="id_genome",
        help="Name of your genome (for example : JAAAIB01",
        default="GENOMESAMPLE",
        required=True,
    )


    # Optional arguments
   
    optional_args = parser.add_argument_group("Optional arguments")
    '''
    # Recovery
    optional_args.add_argument(
        "--tRNAscan_report",
        action="store",
        dest="tRNAscan_report",
        help="path of tRNAscan gff file",
        default="",
        required=False,
    )
    '''
    # Output
    optional_args.add_argument(
        "--output",
        "-o",
        action="store",
        dest="output_dir",
        help="Output directory name",
        default="results",
        required=False,
    )

    args = parser.parse_args()

    try:
        os.mkdir(args.output_dir)
    except:
        print(
            f"\n Output directory {args.output_dir} can not be created, please erase it before launching the programm !"
        )
        sys.exit(1)

    actual_path = os.getcwd()
    os.chdir(args.output_dir)

    global_start = time.perf_counter()

    ####### Creating 1 file per gene

    list_path_file_gene_gff = separiting_file.separiting(args.gff_input)

    # Building dictionnary for trinotate report
    print(list_path_file_gene_gff)
    for sample in range (0,len(list_path_file_gene_gff)):
        line_gene=prepare.writing_gene(list_path_file_gene_gff[sample],args.trinotate_report,sample+1)
        print(line_gene)
        line_mrna=prepare.writing_mRNA(list_path_file_gene_gff[sample],args.trinotate_report,sample+1,args.id_genome)
        list_line_exon=prepare.writing_exon(list_path_file_gene_gff[sample],args.trinotate_report,sample+1,args.id_genome)
        #dic_trinotate = prepare.trinotate_dic(
        #args.trinotate_report, args.targetp_report,list_path_file_gene_gff[sample]
    #)


if __name__ == '__main__':
    run()