#!/usr/bin/env python3

import argparse
import os
import sys
import time
import shutil
from numpy import product
import pandas as pd 

'''
def CDS_dic(trinotate_report: str, targetp_report: str, id_prot_uniprot: str) -> dictionnary :

    table_trinotate = pd.read_csv(trinotate_report, sep="\t", header=None)
    table_targetp = pd.read_csv(targetp_report, sep="\t", header=None)
    table_id_uniprot = pd.read_csv(id_prot_uniprot, sep="\t", header=None)

    dic_trinotate = {}
    list_exon_description = []

    for sample in range(0,table_trinotate.shape[0]):
        #name_of_gene = 'gene=' + table_trinotate.iloc[sample][6]
        print_info = 0
        name_database_and_annot = ''

        #######Number of pfam and Uniprot
        if (table_trinotate.iloc[sample][7] != '.') :
            print_info = 1
            name_database_and_annot = 'Dbxref=PFAM:'+table_trinotate.iloc[sample][7]+',UniProtKB/Swiss-Prot:'+id_swissprot

        ########THMMER
        final_thmmer=''
        information_transmembrane = ''
        if (table_trinotate.iloc[sample][9] != '.') :
            print_info = 1
            transmembrane_data=table_trinotate.iloc[sample][9].split("^")
            number_of_part=transmembrane_data[2].split("=")
            information_transmembrane=transmembrane_data[3].split("=")
            final_thmmer='Note=TransMembrane:'+number_of_part[2]+' ('+information_transmembrane[2]

        ##########SIGNALP
        final_signalp=''
        if (table_trinotate.iloc[sample][8] != '.'):
            print_info = 1
            signalp_data=table_trinotate.iloc[sample][8].split("^")
            begin_signal_peptide=signalp_data[0].split(":")[1]
            end_of_signal_peptide = signalp_data[1]
            final_signalp=')~SECRETED:SignalP('+begin_signal_peptide+'-'+end_of_signal_peptide+'):'
            
        general_go_term=''
        go_component_blastp=''
        go_molecular_function_blastp=''
        go_biological_process_blastp=''
        go_component=[]
        go_biological_process=[]
        go_molecular_function=[]
        go_term_first=[]
        ########## BLASTP
        if (table_trinotate.iloc[sample][13] != '.'):

            print_info = 1
            go_term_b=table_trinotate.iloc[sample][13].split("`")
            go_term_first=[]
            for go in range (0,len(go_term_b)):
                terms_of_go=go_term_b.split("^")
                go_term_first.append(terms_of_go[0])
                if "cellular_component" in go_term_b :
                    go_component.append(go_term_b[0])
                if "molecular_function" in go_term_b :
                    go_molecular_function.append(go_term_b[0])
                if "biological_process" in go_term_b:
                    go_biological_process.append(go_term_b[0])
            if len(go_component != 0 ):
                go_component_blastp='go_component='
                for go_c in range(0,len(go_component)):
                    terms_of_go_c=go_component[go_c].split("^")
                    if go_c != len(go_component)-1:
                        go_component_blastp=go_component_blastp+terms_of_go_c[2]+"|"terms_of_go_c[0]+"||IEA,"
                    if go_c == len(go_component)-1:
                        go_component_blastp=go_component_blastp+terms_of_go_c[2]+"|"terms_of_go_c[0]+"||IEA;"
            if len(go_molecular_function != 0 ):
                go_molecular_function_blastp="go_function="
                for go_m in range (0,len(go_molecular_function)):
                    terms_of_go_m=go_molecular_function[go_m].split("^")
                    if go_c != len(go_molecular_function)-1:
                        go_molecular_function_blastp=go_molecular_function_blastp+terms_of_go_m[2]+"|"terms_of_go_m[0]+"||IEA,"
                    if go_c == len(go_molecular_function)-1:
                        go_molecular_function_blastp=go_molecular_function_blastp+terms_of_go_m[2]+"|"terms_of_go_m[0]+"||IEA;"
            if len(go_biological_process != 0 ):
                go_biological_process_blastp="go_process="
                for go_p in range (0,len(go_biological_process)):
                    terms_of_go_p=go_biological_process[go_p].split("^")
                    if go_c != len(go_biological_process)-1:
                        go_biological_process_blastp=go_biological_process_blastp+terms_of_go_p[2]+"|"terms_of_go_p[0]+"||IEA,"
                    if go_c == len(go_biological_process)-1:
                        go_biological_process_blastp=go_biological_process_blastp+terms_of_go_p[2]+"|"terms_of_go_p[0]+"||IEA;"
            #concat_go_term=",".join(go_term_first)
            #general_go_term="Ontology_term="+concat_go_term+";gbkey=CDS;"+go_component_blastp+go_molecular_function_blastp+go_biological_process_blastp             

        #########EGGNOG
        if (table_trinotate.iloc[sample][10] != '.'):
            print_info = 1
            ###y'en aura jamais normalement

        #######KEGG
        if (table_trinotate.iloc[sample][11] != '.'):
            #KEGG
            ######KEGG ne doit pas apparaître dans le gff ########

        ########PFAM
        if (table_trinotate.iloc[sample][14] != '.'):
            ####GO term for pfam
            print_info = 1
            go_term_b=table_trinotate.iloc[sample][14].split("`")
            for go in range (0,len(go_term_b)):

                terms_of_go=go_term_b.split("^")
                #go_term_first.append(terms_of_go[0])

                if "cellular_component" in go_term_b :
                    if terms_of_go[0] in go_term_first :
                        print('in common')
                    else :
                        go_term_first.append(terms_of_go[0])
                        go_component.append(go_term_b[0])
                        
                if "molecular_function" in go_term_b :
                    if terms_of_go[0] in go_term_first :
                        print('in common')
                    else :
                        go_term_first.append(terms_of_go[0])
                        go_molecular_function.append(go_term_b[0])
                    

                if "biological_process" in go_term_b:
                    if terms_of_go[0] in go_term_first :
                        print('in common')
                    else :
                        go_term_first.append(terms_of_go[0])
                        go_biological_process.append(go_term_b[0])
                    
            if len(go_component != 0 ):
                if go_component_blastp == '':
                    go_component_blastp='go_component='
                for go_c in range(0,len(go_component)):
                    terms_of_go_c=go_component[go_c].split("^")
                    if go_c != len(go_component)-1:
                        go_component_blastp=go_component_blastp+terms_of_go_c[2]+"|"terms_of_go_c[0]+"||IEA,"
                    if go_c == len(go_component)-1:
                        go_component_blastp=go_component_blastp+terms_of_go_c[2]+"|"terms_of_go_c[0]+"||IEA;"
            if len(go_molecular_function != 0 ):
                if go_molecular_function_blastp == '':
                    go_molecular_function_blastp="go_function="
                for go_m in range (0,len(go_molecular_function)):
                    terms_of_go_m=go_molecular_function[go_m].split("^")
                    if go_c != len(go_molecular_function)-1:
                        go_molecular_function_blastp=go_molecular_function_blastp+terms_of_go_m[2]+"|"terms_of_go_m[0]+"||IEA,"
                    if go_c == len(go_molecular_function)-1:
                        go_molecular_function_blastp=go_molecular_function_blastp+terms_of_go_m[2]+"|"terms_of_go_m[0]+"||IEA;"
            if len(go_biological_process != 0 ):
                if go_biological_process_blastp == '':
                    go_biological_process_blastp="go_process="
                for go_p in range (0,len(go_biological_process)):
                    terms_of_go_p=go_biological_process[go_p].split("^")
                    if go_c != len(go_biological_process)-1:
                        go_biological_process_blastp=go_biological_process_blastp+terms_of_go_p[2]+"|"terms_of_go_p[0]+"||IEA,"
                    if go_c == len(go_biological_process)-1:
                        go_biological_process_blastp=go_biological_process_blastp+terms_of_go_p[2]+"|"terms_of_go_p[0]+"||IEA;"

            concat_go_term=",".join(go_term_first)
            general_go_term="Ontology_term="+concat_go_term+";gbkey=CDS;"+go_component_blastp+go_molecular_function_blastp+go_biological_process_blastp
        
        

        if (print_info == 1):
            note_gff = +';gbkey=CDS;;locus_tag=BGZ75_000005;orig_transcript_id=gnl|WGS:JAAAIB|BGZ75_000006-T1_mrna;'+

'''

def writing_gene(gff_file : chr , trinotate_report : chr , number_of_gene : int) -> list:
    table_trinotate = pd.read_csv(trinotate_report, sep="\t", header=None)
    #table_id_uniprot = pd.read_csv(id_prot_uniprot, sep="\t", header=None), id_prot_uniprot : chr 
    table_gff = pd.read_csv(gff_file, sep="\t", header=None)
    line_gene = table_gff.loc[table_gff[3]=='gene'].index.values
    list_gene_gff = table_gff.iloc[line_gene , 9][0].split("=")
    ID_gene_gff = list_gene_gff[1]
    gene_element=''
    name_of_gene=''
    line_trinotate_gene=table_trinotate.iloc[number_of_gene,6].split('^')
    trinotate_gene_name=line_trinotate_gene[0].split('_')
    if table_trinotate.iloc[number_of_gene,6] != '.':
        name_of_gene=trinotate_gene_name[0]
        gene_element="gene="+name_of_gene+";"
    else :
        name_of_gene=ID_gene_gff
    text_gene = "ID=gene-"+ID_gene_gff+";Name="+name_of_gene+";gbkey=Gene;"+gene_element+"gene_biotype=protein_coding;locus_tag="+ID_gene_gff
    return text_gene

def writing_mRNA ( gff_file : chr , trinotate_report : chr , number_of_gene : int , name_of_genome : chr) -> chr :
    table_trinotate = pd.read_csv(trinotate_report, sep="\t", header=None)
    #table_id_uniprot = pd.read_csv(id_prot_uniprot, sep="\t", header=None), id_prot_uniprot : chr 
    table_gff = pd.read_csv(gff_file, sep="\t", header=None)
    line_gene = table_gff.loc[table_gff[3]=='gene'].index.values
    list_gene_gff = table_gff.iloc[line_gene , 9][0].split("=")
    ID_transcript_gff = list_gene_gff[1]
    ID_gene_gff = list_gene_gff[1]
    mRNA_element=''
    gene_element=''
    name_of_gene=''
    line_trinotate_gene = table_trinotate.iloc[number_of_gene,6].split('^')
    line_trinotate_gene_2 = table_trinotate.iloc[number_of_gene,6].split('^')
    line_trinotate_gene_3 = table_trinotate.iloc[number_of_gene,8].split('^')
    line_trinotate_gene_4 = table_trinotate.iloc[number_of_gene,9].split('^')
    trinotate_gene_name = line_trinotate_gene[0].split('_')
    ####emplacement numero 1
    if table_trinotate.iloc[number_of_gene,6] != '.':
        complete_name_of_gene = line_trinotate_gene[5].split('=')
        name_of_gene=trinotate_gene_name[0]
        gene_element="gene="+name_of_gene+";"
        product_mrna="product="+complete_name_of_gene[1]
    else :
        name_of_gene=ID_gene_gff
        product_mrna="product=hypothetical protein"
    text_mrna="ID=rna-gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1-mrna;Parent="+ID_gene_gff+";gbkey=mRNA;"+gene_element+"locus_tag="+ID_gene_gff+";orig_protein_id=gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1;orig_transcript_id:gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1-mrna;"+product_mrna
    print(text_mrna)
    return(text_mrna)
    #########function OK############

def writing_exon( gff_file : chr , trinotate_report : chr , number_of_gene : int , name_of_genome : chr) -> list :
    table_trinotate = pd.read_csv(trinotate_report, sep="\t", header=None)
    #table_id_uniprot = pd.read_csv(id_prot_uniprot, sep="\t", header=None), id_prot_uniprot : chr 
    list_text_exon=[]
    table_gff = pd.read_csv(gff_file, sep="\t", header=None)
    line_gene = table_gff.loc[table_gff[3]=='exon'].index.values
    for sample_exon in range(0,len(line_gene)):
        list_gene_gff = table_gff.iloc[line_gene[sample_exon], 9].split("=")
        list_ID_transcript_gff = list_gene_gff[2].split('.')
        ID_gene_gff = list_ID_transcript_gff[0]
        gene_element=''
        name_of_gene=''
        line_trinotate_gene = table_trinotate.iloc[number_of_gene,6].split('^')
        trinotate_gene_name = line_trinotate_gene[0].split('_')
        if table_trinotate.iloc[number_of_gene,6] != '.':
            complete_name_of_gene = line_trinotate_gene[5].split('=')
            name_of_gene=trinotate_gene_name[0]
            gene_element="gene="+name_of_gene+";"
            product_mrna="product="+complete_name_of_gene[1]
        else :
            name_of_gene=ID_gene_gff
            product_mrna="product=hypothetical protein"
        text_exon="ID=exon-gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1-mrna-"+str(sample_exon+1)+";Parent=rna-gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1-mrna;gbkey=mRNA;"+gene_element+"locus_tag="+ID_gene_gff+";orig_protein_id=gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1;orig_transcript_id:gnl|WGS:"+name_of_genome+"|"+ID_gene_gff+"-T1-mrna;"+product_mrna
        list_text_exon.append(text_exon)
    return list_text_exon
    ####### FUNCTION DONE###########
'''
def gene_dic ( trinotate_report: str, targetp_report: str, id_prot_uniprot: str) -> dictionnary :


def mRNA_dic ( trinotate_report: str, targetp_report: str, id_prot_uniprot: str) -> dictionnary :


def exon_dic ( trinotate_report: str, targetp_report: str, id_prot_uniprot: str) -> dictionnary :
'''  
    