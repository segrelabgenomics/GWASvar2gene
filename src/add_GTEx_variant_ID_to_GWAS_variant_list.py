#!/usr/bin/env python3

# Author: John Rouhana
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School
# Date: 2020

import sys
import pandas as pd

## INPUT FILES:
# List of variants (chr number, position, variant ID (CHR_POS_REF_ALT_BUILD) from QTL database, in this case GTEx v8:
# Required columns: [1]CHROM     [2]POS  [3]REF  [4]ALT  [5]ID
GTEx = sys.argv[1]
#GTEx = '../data/GTEx_v8_HG38_all_variants.tsv.gz'

# List of variants that will be annnotated with eQTLs and/or sQTLs, e.g., variants from GWAS catalog
# Required columns:  'CHR_ID', 'CHR_POS', 'STRONGEST SNP-RISK ALLELE' (e.g., rs1925953-T)
GWAS = sys.argv[2]
#GWAS = '../example/gwas_catalog_sample_file.tsv'

## OUTPUT FILE
out_file = sys.argv[3]
#out_file = '../output/out_GWAS_variant_IDs_added.tsv'

GTEx_df = pd.read_csv(GTEx, sep='\t', compression='gzip')
GWAS_df = pd.read_csv(GWAS, sep='\t', low_memory=False)

GTEx_df.columns = ['CHR_ID', 'CHR_POS', 'REF', 'ALT', 'GWAS_variant']
# Assuming genome build 38
GWAS_df["CHR_ID"] = GWAS_df["CHR_ID"].apply(lambda x: "chr"+str(x) if not str(x).startswith("chr") else str(x))

def ALT_splitter(x):
    if '-' in x:
        ret = x.split('-')[1]
    else:
        ret = None
    if ret == '?':
        ret = None
    return(ret)

GWAS_df['RISK_ALLELE'] = GWAS_df['STRONGEST SNP-RISK ALLELE'].apply(lambda x: ALT_splitter(x))

GWAS_df['CHR_ID'] = GWAS_df['CHR_ID'].astype(str)
GWAS_df['CHR_POS'] = GWAS_df['CHR_POS'].astype(str)
GTEx_df['CHR_POS'] = GTEx_df['CHR_POS'].astype(str)
GTEx_df['CHR_ID'] = GTEx_df['CHR_ID'].astype(str)


dropped_indels = GTEx_df.loc[GTEx_df.REF.apply(lambda x: len(x) == 1)].copy()
new_df = GWAS_df.merge(dropped_indels[['CHR_ID', 'CHR_POS', 'REF']], how='left', on=['CHR_ID', 'CHR_POS']).copy()


ref_is_risk = new_df.loc[new_df['RISK_ALLELE'] == new_df['REF']].copy()
alt_is_risk = new_df.loc[new_df['RISK_ALLELE'] != new_df['REF']].copy()

alt_is_known = alt_is_risk.loc[alt_is_risk.RISK_ALLELE.apply(lambda x: x != None)].copy()
alt_is_unknown = alt_is_risk.loc[alt_is_risk.RISK_ALLELE.apply(lambda x: x == None)].copy()


final_df = []
#For alt is known, we already have REF and ALT
alt_is_known['ALT'] = alt_is_known['RISK_ALLELE']
final_df.append(alt_is_known)

#For REF is risk, simply grab all alternatives as variants of interest
ref_is_risk = ref_is_risk.merge(dropped_indels[['CHR_ID', 'CHR_POS', 'REF', 'ALT']], how='left', on=['CHR_ID', 'CHR_POS', 'REF']).copy()
final_df.append(ref_is_risk)

#for alt is unknown, grab all possible variants at that position
alt_is_unknown.drop('REF', axis=1, inplace=True)
alt_is_unknown = alt_is_unknown.merge(GTEx_df[['CHR_ID', 'CHR_POS', 'REF', 'ALT']], how='left', on=['CHR_ID', 'CHR_POS']).copy()
final_df.append(alt_is_unknown)

final_df = pd.concat(final_df)
final_df = final_df.drop_duplicates().copy()

final_df['GWAS_variant'] = final_df['CHR_ID'] + '_' + final_df['CHR_POS'] + '_' + final_df['REF'] + '_' + final_df['ALT'] + '_' + 'b38'
final_df.to_csv(out_file, sep='\t', index=False)
