# script_trinotate_to_gff
A little script to print a gff3 file 

Script to extract name of the protein in  Trinotate report :

`` 
awk -F "\t" '{print $7}' report_trino.txt | awk -F "^" '{print $1}' | sort | uniq > prot_id 
``

``
(Delete '.' and the name of the column in prot_id) 
`` 

``
aria2c https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping.dat.gz
``

``
grep -f "prot_id"   idmapping.dat > prot_id_final
`` 


I will publish test version and test file after my publication but you can try to understand how to used it with this command line :

``
script_gff_builder_Trinotate.py --gff_input augustus_annotation_corrected_with_agat.gff3 --trinotate_report report_trino.txt --output results_test --id_genome AAAAAAQ --uniprot_cor output_ID_correspondance_2 --prot_accession TESTF --prot_number 000555555771
`` 
