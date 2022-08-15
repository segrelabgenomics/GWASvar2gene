#!/bin/bash

## INPUT FILE
# Table with variant ID and chr position of all variants analyzed in GTEx v8. File can be downloaded from: https://zenodo.org/deposit/6994426  
GTEx='../data/GTEx_v8_HG38_all_variants.tsv.gz'
#GWAS=sys.argv[2]
GWAS = '../example/gwas_catalog_sample_file.tsv'

## OUTPUT FILE
#out_file = sys.argv[3]
out_file = '../output/out_GWAS_variant_IDs_added.tsv'

module purge
module load pandas/2.5.0

./add_GTEx_variant_ID_to_GWAS_variant_list.py $GTEx $GWAS $out_file
