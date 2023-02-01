#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
from isort import file
import pandas as pd 

def separiting (gff_file : file) -> list:
    table_gff = pd.read_csv(gff_file, sep="\t", header=None)
    list_path_file_gene_gff=[]
    counter=0
    line_gene=table_gff[table_gff[2]=="gene"].index.values
    print(len(line_gene))
    for sample in range(0,len(line_gene)-1):
        counter+=1
        if sample == len(line_gene)-1 :
            small_gene_df=table_gff.iloc[line_gene[sample]:-1,]
        if sample == len(line_gene)-2 :
            small_gene_df=table_gff.iloc[line_gene[sample]:line_gene[-1],]
        else :
            small_gene_df=table_gff.iloc[line_gene[sample]:line_gene[sample+1],]
        small_gene_df.to_csv(path_or_buf="mini_gff"+str(counter), sep="\t", header=False)
        list_path_file_gene_gff.append(os.path.abspath("mini_gff"+str(counter)))
    return list_path_file_gene_gff
