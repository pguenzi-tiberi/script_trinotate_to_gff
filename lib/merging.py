#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
import pandas as pd 


def merging_in_each_file(gff_report : chr ,trinotate_report : chr ,line_gene : chr ,line_mrna : chr ,list_line_exon : list ,list_line_CDS : list ):
    table_trinotate = pd.read_csv(trinotate_report, sep="\t", header=None)
    table_gff = pd.read_csv(gff_report, sep="\t", header=None)
    length_cds=len(list_line_CDS)
    length_exon=len(list_line_exon)
    for line in range (0,table_gff.shape[0]):
        if table_gff.iloc[line,3]== 'gene' :
            table_gff.loc[line,9] = line_gene
        if table_gff.iloc[line,3]== 'mRNA' :
            table_gff.loc[line,9] = line_mrna
        if table_gff.iloc[line,3]== 'exon' :
            table_gff.loc[line,9] = list_line_exon[len(list_line_exon)-length_exon]
            length_exon-=1
        if table_gff.iloc[line,3]== 'CDS' :
            table_gff.loc[line,9] = list_line_CDS[len(list_line_CDS)-length_cds]
            length_cds-=1
    return table_gff