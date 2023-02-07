#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
import pandas as pd
from sympy import im 
from lib import merging
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
        "--uniprot_cor",
        action="store",
        dest="uniprot_cor",
        help="TSV file (first column = ID Uniprot, second column = Name of protein)",
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

    mandatory_args.add_argument(
        "--prot_accession",
        action="store",
        dest="prot_accession",
        help="For annotation. It is for protein_id in CDS part",
        default='',
        required=True,
    )

    mandatory_args.add_argument(
        "--prot_number",
        action="store",
        dest="prot_number",
        help="For annotation. It is for protein_id in CDS part",
        default='',
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
    prot_number=float(args.prot_number)
    final_table=''
    for sample in range (0,len(list_path_file_gene_gff)):
        if sample >0:
            prot_number+=1
        line_gene=prepare.writing_gene(list_path_file_gene_gff[sample],args.trinotate_report,sample+1)
        line_mrna=prepare.writing_mRNA(list_path_file_gene_gff[sample],args.trinotate_report,sample+1,args.id_genome)
        list_line_exon=prepare.writing_exon(list_path_file_gene_gff[sample],args.trinotate_report,sample+1,args.id_genome)
        list_line_CDS=prepare.writing_CDS(list_path_file_gene_gff[sample],args.trinotate_report,sample+1,args.id_genome,args.prot_accession,prot_number,args.uniprot_cor)
        if sample == 0:
            final_table=merging.merging_in_each_file(list_path_file_gene_gff[sample],args.trinotate_report,line_gene,line_mrna,list_line_exon,list_line_CDS)
        else :
            gene_table=merging.merging_in_each_file(list_path_file_gene_gff[sample],args.trinotate_report,line_gene,line_mrna,list_line_exon,list_line_CDS)
            final_table=pd.concat([final_table,gene_table],axis=0)
        print(final_table)
        if sample == 80 :
            final_table.to_csv(path_or_buf="final_table.tsv", sep="\t", header=False)
            break


if __name__ == '__main__':
    run()