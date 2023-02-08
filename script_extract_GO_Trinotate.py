#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
import pandas as pd

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
   
    args = parser.parse_args()
    actual_path = os.getcwd()

    global_start = time.perf_counter()

    table_trinotate = pd.read_csv(args.trinotate_report, sep="\t", header=None)

    file_finished = open("Go_terms_extract.txt", "w+")
    file_finished.write("gene_ID    ontology_type   nb_GO   desc\n")
    for go in range (0,table_trinotate.shape[0]):
        name_of_gene=table_trinotate.iloc[go,0]
        go_term_blastp=table_trinotate.iloc[go][13].split("`")
        go_term_pfam=table_trinotate.iloc[go][14].split("`")
        if go_term_blastp[0] != '.': 
            for x in range(0,len(go_term_blastp)):
                terms_of_go=go_term_blastp[x].split("^")
                go_ID=terms_of_go[0].split(":")
                if terms_of_go.count("cellular_component") > 0 :
                    line=name_of_gene+"\tCC\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
                if terms_of_go.count("molecular_function") > 0 :
                    line=name_of_gene+"\tMF\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
                if terms_of_go.count("biological_process") > 0 :
                    line=name_of_gene+"\tBP\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
        if go_term_pfam[0] != '.': 
            for x in range(0,len(go_term_pfam)):
                terms_of_go=go_term_pfam[x].split("^")
                go_ID=terms_of_go[0].split(":")
                if terms_of_go.count("cellular_component") > 0 :
                    line=name_of_gene+"\tCC\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
                if terms_of_go.count("molecular_function") > 0 :
                    line=name_of_gene+"\tMF\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
                if terms_of_go.count("biological_process") > 0 :
                    line=name_of_gene+"\tBP\t"+go_ID[1]+"\t"+terms_of_go[2]+"\n"
                    file_finished.write(line)
    
            
    print(
        f"\n Total running time : {float(time.perf_counter() - global_start)} seconds"
    )

    

if __name__ == '__main__':
    run()