# Author: John Rouhana
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School
# Date: 2020

import argparse
import sys
import gzip
import os
import tempfile
import pandas as pd
import numpy as np
from new_tools import (str2bool, ERRORPARSING, convert_variant_format, call_plink,
                       dataframe, clean_plink, deep_clean_plink, g_open, 
                       parse_gtf, gene)

import re
import subprocess
import glob
import dask.dataframe as dd


def parse_eqtl_file(sig_file, mapping_df, #LD_variants_dict,
                    gene_dict, tissue, eGene_dict, type='eQTL'):
    """
    Reads in an eQTL file
    Returns rows that for which the variant are in LD table as LD variant

    sig_file: file with significant variant-gene pairs
    mapping_df: a df with GWAS_variant, proxy_pairing_id, and LD_variant
    LD_variants_dict: dict pointing to sets by chromosome
    gene_dict: gene_id to gene_name dict
    tissue: name of tissue
    eGene_dict: dict indicating best eVariant for each eGene by tissue
    """
    #Define DataFrame to fill
    eVariants={}
    gene_ids={}
    tissues={}
    tss_distances={}
    pval_nominals={}
    slopes={}
    slope_ses={}
    pval_nominal_thresholds={}
    min_pval_nominals={}
    pval_betas={}
    chrs={}
    eVariant_bp_poss={}
    gene_names={}
    gene_id_truncateds={}
    best_eqtls={}

    if type == 'sQTL': #If we're working with sQTLs
        cluster_ids = {}

    #Parse through sig_file, get what we need
    with g_open(sig_file) as f:
        for i, line in enumerate(f):
            if i != 0: #skip header
                row = line.strip().split('\t')
                eVariant = row[0]
                chr = eVariant.split('_')[0].replace('chr', '')
#                if chr not in LD_variants_dict.keys():
#                    continue
                #TODO: Attempt this if it's not fast enough
#                elif eVariant not in LD_variants_dict[chr]:
#                    continue
                gene_id = row[1]
                if type == 'sQTL':
                    cluster_id = gene_id
                    gene_id = cluster_id.split(':')[4]

                eVariants[str(i)] = eVariant
                gene_ids[str(i)] = gene_id
                tissues[str(i)] = tissue
                tss_distances[str(i)] = row[2]
                pval_nominals[str(i)] = row[6]
                slopes[str(i)] = row[7]
                slope_ses[str(i)] = row[8]
                pval_nominal_thresholds[str(i)] = row[9]
                min_pval_nominals[str(i)] = row[10]
                pval_betas[str(i)] = row[11]
                chrs[str(i)] = eVariant.split('_')[0]
                eVariant_bp_poss[str(i)] = eVariant.split('_')[1]
                try:
                    gene_names[str(i)] = gene_dict[gene_id.split('.')[0]].gene_name
                except:
                    gene_names[str(i)] = None
                try:
                    gene_id_truncateds[str(i)] = gene_id.split('.')[0]
                except:
                    gene_id_truncateds[str(i)] = None
                try:
                    best_eqtls[str(i)] = determine_best_eqtl(eVariant, gene_id, tissue, eGene_dict) #Returns True/False
                except:
                    best_eqtls[str(i)] = None
                   
                if type == 'sQTL':
                    cluster_ids[str(i)] =  cluster_id

    ret_df = dd.from_pandas(pd.DataFrame({'eVariant':list(eVariants.values()), 
                                          'gene_id':list(gene_ids.values()),
                                          'tissue':list(tissues.values()),
                                          'tss_distance':list(tss_distances.values()),
                                          'pval_nominal':list(pval_nominals.values()),
                                          'slope':list(slopes.values()),
                                          'slope_se':list(slope_ses.values()),
                                          'pval_nominal_threshold':list(pval_nominal_thresholds.values()),
                                          'min_pval_nominal':list(min_pval_nominals.values()),
                                          'pval_beta':list(pval_betas.values()),
                                          'chr':list(chrs.values()),
                                          'eVariant_bp_pos':list(eVariant_bp_poss.values()),
                                          'gene_name':list(gene_names.values()),
                                          'gene_id_truncated':list(gene_id_truncateds.values()),
                                          'best_eqtl':list(best_eqtls.values())}), npartitions=20)
                
    if type == 'sQTL': #Rearrange index
        ret_df['cluster_id'] = dd.from_array(np.array(list(cluster_ids.values())))
        ret_df = ret_df[['eVariant', 'gene_id', 'cluster_id', 'tissue', 'tss_distance',
                         'pval_nominal', 'slope', 'slope_se', 'pval_nominal_threshold',
                         'min_pval_nominal', 'pval_beta', 'chr', 'eVariant_bp_pos',
                         'gene_name', 'gene_id_truncated', 'best_eqtl']].copy()

        ret_df = ret_df.rename(columns={'eVariant':'sVariant',
                                        'eVariant_bp_pos':'sVariant_bp_pos',
                                        'best_eqtl':'best_sqtl'}).copy()
    ret_df.reset_index(drop=True)
    #Get rows that are LD variants to our GWAS
    if type == 'eQTL':
        ret_df = ret_df.merge(mapping_df[['LD_variant', 'GWAS_variant', 'proxy_pairing_ID']],
                              how='inner', left_on='eVariant', right_on='LD_variant').copy()
        ret_df = ret_df.drop('LD_variant', axis=1).copy()
    elif type == 'sQTL':
        ret_df = ret_df.merge(mapping_df[['LD_variant', 'GWAS_variant', 'proxy_pairing_ID']],
                                          how='inner', left_on='sVariant', right_on='LD_variant').copy()

        ret_df = ret_df.drop('LD_variant', axis=1).copy()

    #Turn to 1-n-n returns 
    if type == 'eQTL':
        ret_1 = ret_df[['eVariant', 'gene_id', 'tss_distance', 'chr', 'eVariant_bp_pos', 
                        'gene_name', 'gene_id_truncated', 'GWAS_variant',
                        'proxy_pairing_ID']].drop_duplicates().copy()
        ret_2 = ret_df[['eVariant', 'gene_id', 'tss_distance', 'chr', 'eVariant_bp_pos',
                        'gene_name', 'gene_id_truncated', 'tissue', 'best_eqtl', 
                        'pval_nominal', 'slope', 'slope_se', 'pval_nominal_threshold',
                        'min_pval_nominal', 'pval_beta']].drop_duplicates().copy()
    elif type == 'sQTL':
        ret_1 = ret_df[['sVariant', 'gene_id', 'cluster_id', 'tss_distance', 'chr', 'sVariant_bp_pos',
                        'gene_name', 'gene_id_truncated', 'GWAS_variant',
                        'proxy_pairing_ID']].drop_duplicates().copy()
        ret_2 = ret_df[['sVariant', 'gene_id', 'cluster_id', 'tss_distance', 'chr', 'sVariant_bp_pos',
                        'gene_name', 'gene_id_truncated', 'tissue', 'best_sqtl',
                        'pval_nominal', 'slope', 'slope_se', 'pval_nominal_threshold',
                        'min_pval_nominal', 'pval_beta']].drop_duplicates().copy()

    return(ret_1, ret_2)

def determine_best_eqtl(eVariant, gene_id, tissue, eGene_dict):
    """
    Returns True/False indicating whether this is
    the best eVariant-eGene-tissue combination
    for a particular gene. Works for sGenes as well

    eGene_dict: dictionary containing desired results
    """
    tissue_dict = eGene_dict[tissue]
    if tissue_dict[gene_id] == eVariant:
        return(True)
    else:
        return(False)
               
def build_best_eGene_dict(tissue_file, type='eQTL'):
    """
    tissue_file: file mapping tissues to significant eQTL files

    returns a dict of 1 dict per tissue indicating best 
    eVariant for each eGene in tissue.
    """
    eGene_dict = {}

    with g_open(tissue_file) as f:
        for i, line in enumerate(f):
            row = line.strip().split('\t')
            tissue = row[0]
            eGene_dict[tissue] = {tissue:{}} #Initiate dict of dicts
            try:
                if type == 'eQTL':
                    egene_file = row[2]
                elif type == 'sQTL':
                    egene_file = row[4]
                with g_open(egene_file) as z:
                    for idx, l in enumerate(z):
                        if idx == 0:
                            continue
                        l_row = l.strip().split('\t')
                        gene_id = l_row[0]
                        if type == 'sQTL':
                             gene_id = gene_id.split(':')[4]
                        eVariant = l_row[11]
                        eGene_dict[tissue][gene_id] = eVariant
            except:
                print("No best eGene file for "+row[0]+". Continuing.")
    return(eGene_dict)





parser = argparse.ArgumentParser(description="In-depths annotation for GWAS variants from the GWAS Catalog")
parser.add_argument('--GWAS_variants', type=str, required=True,
                    help='NHGRI-EBI Catalog without ontology annotations and with a "GWAS_variant" column added to annotate variant in "chr1_12345_A_T_b38" format, or the likes.')
parser.add_argument('--primary_plink_directory', type=str, required=True,
                    help='Directory where the primary set of PLINK files to be used is stored.')
parser.add_argument('-primary_plink_prefix', type=str, default='plink1',
                    help="Name of the study that the primary PLINK files originated from.")
parser.add_argument('-primary_subset', type=str, default=None, 
                    help='Subset of samples from the primary PLINK files to use for LD analysis.')
parser.add_argument('-primary_ldr', type=float, default=0.5, 
                    help='LD r^2 value to use for annotation cutoff from primary PLINK files.')
parser.add_argument('-secondary_plink_directory', type=str, default=None,
                    help='If variants are not in the primary PLINK files, a second set of PLINK files from a separate study can be included to find proxies.')
parser.add_argument('-secondary_plink_prefix', type=str, default='plink2', 
                    help="Name of the study that the secondary PLINK files originated from.")
parser.add_argument('-secondary_subset', type=str, default=None, 
                    help='Subset of samples from the secondary PLINK files to use for finding proxy variants')
parser.add_argument('-secondary_ldr', type=float, default=0.8,
                    help='LD r^2 value to use for finding proxies to variants with no LD variants in the primary PLINK fileset.')
parser.add_argument('-use_db', type=str2bool, default=False,
                    help='Whether to write outputs to a mySQL database.')
parser.add_argument('-db_host', type=str, default=None,
                    help='Host of the database to write to.')
parser.add_argument('-db_user', type=str, default=None,
                    help='Username in database.')
parser.add_argument('-db_pass', type=str, default=None,
                    help='Password in database.')
parser.add_argument('-db_name', type=str, default='GWASVar2Gene_output',
                    help='Name of database to connect to.')
parser.add_argument('-db_overwrite', type=str, default=False,
                    help='Whether to overwrite database if it already exists.')
parser.add_argument('-plink_memory_allocation', type=int, default=4096,
                    help='Amount of memory to allocate toward PLINK running. Note: this restriction does not apply to the rest of the program.')
parser.add_argument('-ld_MAF_cut', type=str2bool, default=True,
                    help='Whether to exclude candidate LD variants with MAF < 0.01.')
parser.add_argument('-proxy_MAF_cut', type=str2bool, default=True,
                    help='Whether to exclude candidate proxies with MAF < 0.01.') 
parser.add_argument('-GENCODE', type=str, default=None,
                    help='GENCODE gtf or gtf.gz file for gene annotation.')
parser.add_argument('-all_QTL_analysis_variants', type=str, default=None, 
                    help='All variants used in computing QTLs.')
parser.add_argument('-tissue_file', type=str, default=None,
                    help='A tab-delimited file with appropriate tissue-to-eQTL/sQTL mapping, as indicated by documentation..')
parser.add_argument('-output_dir', type=str, default='',
                    help='Desired output directory.')
parser.add_argument('-output_prefix', type=str, default='GWASVar2Gene_output',
                    help='Prefix to assign to output. Do not include desired directory output.')
parser.add_argument('-variants_only', type=str2bool, default=False, 
                    help="If true, will run a minimalist run of GWASVar2Gene using only PLINK and stopping after the variants have been annotated for LD. Expects a file with GWAS_variants column.")
parser.add_argument('-ldkb', type=int, default=5000,
                    help="Kilobase window to search for LD.")
parser.add_argument('-development_skip', type=str2bool, default=False,
                    help="If True, assumes this is a development run where the proxy/trait table have already been generated under output name. Loads the files and skips LD calculations.")
#Parse arguments
args = parser.parse_args()

#Check that necessary files and folders exist
ERRORPARSING(args)

output_dir = args.output_dir
gwas_file = args.GWAS_variants
gtf_file = args.GENCODE
QTL_variants = args.all_QTL_analysis_variants
primary_dir = args.primary_plink_directory
secondary_dir = args.secondary_plink_directory
primary_prefix = args.primary_plink_prefix
secondary_prefix = args.secondary_plink_prefix
primary_subset = args.primary_subset
secondary_subset = args.secondary_subset
memory = args.plink_memory_allocation
primary_ldr = args.primary_ldr
secondary_ldr = args.secondary_ldr
ldkb = args.ldkb
variants_only = args.variants_only
output_prefix = args.output_prefix
tissue_file = args.tissue_file
d_skip = args.development_skip

if not d_skip:
    #Verify output folder exists
    if output_dir[-1] != "/":
        output_dir = output_dir+"/"

    if not primary_dir.endswith('/'):
        primary_dir = primary_dir + '/'

    if secondary_dir is not None:
        if not secondary_dir.endswith('/'):
            secondary_dir = secondary_dir + '/'

    if not ((os.path.exists(output_dir)) or (output_dir =='')):
        os.mkdir(output_dir)

    #Read in variant data. Assume variants_only is False for now.
    full_GWAS_table = pd.read_csv(gwas_file, sep='\t',
                                  dtype={'PUBMEDID':str, 'REPLICATION SAMPLE SIZE':str,
                                  'CHR_ID':str, 'UPSTREAM_GENE_ID':str, 'DOWNSTREAM_GENE_ID':str,
                                  'CHR_POS':str, 'SNP_ID_CURRENT':str, 'RISK ALLELE FREQUENCY':str,
                                  'P-VALUE':str})#, encoding='iso-8859-1')

    #Get rid of lines where GWAS_variant is not given
    full_GWAS_table.drop_duplicates(inplace=True)
    full_GWAS_table = full_GWAS_table.loc[~(full_GWAS_table.GWAS_variant.isnull())].copy()
    variants = list(full_GWAS_table.GWAS_variant)

    #Set up traits output table
    traits = full_GWAS_table[['GWAS_variant', 'DISEASE/TRAIT', 'PUBMEDID', 
                              'SNPS', 'REPORTED GENE(S)', 'MAPPED_GENE', 'UPSTREAM_GENE_ID', 
                              'DOWNSTREAM_GENE_ID', 'UPSTREAM_GENE_DISTANCE', 'DOWNSTREAM_GENE_DISTANCE',
                              'STRONGEST SNP-RISK ALLELE', 'CONTEXT', 'INTERGENIC', 'P-VALUE',
                              'PVALUE_MLOG', 'OR or BETA', 'INITIAL SAMPLE SIZE'
                             ]].copy()
    traits.drop_duplicates(inplace=True)
    traits.reset_index(drop=True, inplace=True)

    traits.columns = ['GWAS_variant', 'trait', 'pubmed_id', 'rs_id', 'GWAS_reported_genes', 'GWAS_mapped_genes',
                      'GWAS_variant_upstream_gene_id', 'GWAS_variant_downstream_gene_id', 'GWAS_variant_upstream_gene_distance',
                      'GWAS_variant_downstream_gene_distance', 'GWAS_strongest_SNP_risk_allele', 'GWAS_variant_context',
                      'GWAS_variant_intergenic', 'GWAS_p_value', 'GWAS_p_value_mlog', 'GWAS_OR_or_BETA', 'Sample_Ancestries']


    #Determine if we're going to be using a second study
    second_study = False
    if secondary_dir is not None:
        second_study = True

    #Determine desired variant format
    #Pick a file, then sample a variant from it
    primary_files = [os.path.join(primary_dir, f) for f in os.listdir(primary_dir) if os.path.isfile(os.path.join(primary_dir, f))]
    bimone = [f for f in primary_files if f.endswith('bim')][0]
    with open(bimone) as f:
        first_line = f.readline()
    example_primary_variant = first_line.split('\t')[1]
    
    #Fix GWAS variants in traits
    traits['GWAS_variant'] = [convert_variant_format(x, example_variant=example_primary_variant) for x in list(traits['GWAS_variant'])]
    
    #Repeat for secondary variants
    if second_study:
        secondary_files = [os.path.join(secondary_dir, f) for f in os.listdir(secondary_dir) if os.path.isfile(os.path.join(secondary_dir, f))]
        bimtwo = [f for f in secondary_files if f.endswith('bim')][0]
        with open(bimtwo) as f:
            first_line = f.readline()
        example_secondary_variant = first_line.split('\t')[1]
        #Convert variants to second format needed for PLINK
        secondary_variants = list(set([convert_variant_format(x, example_variant=example_secondary_variant) for x in variants]))



    #convert variants to first format needed for PLINK
    print("Converting variants to format needed for PLINK. NOTE THAT THIS TOOL DOES NOT LIFTOVER VARIANTS.")
    variants = [convert_variant_format(x, example_variant=example_primary_variant) for x in variants]

    #Start setting up table. Get chr and pos_bp from variants
    GWAS_variant_table = pd.DataFrame({'GWAS_variant':list(set(variants))})
    GWAS_variant_table['chr'] = GWAS_variant_table['GWAS_variant'].apply(lambda x: re.split(':|_', x)[0].replace('chr', ''))
    GWAS_variant_table['pos_bp'] = GWAS_variant_table['GWAS_variant'].apply(lambda x: re.split(':|_', x)[1])
    #Lines needed if second study is in use
    if second_study:
        GWAS_variant_table[secondary_prefix+'_proxies_needed'] = False
        GWAS_variant_table[secondary_prefix+'_proxies_found'] = np.nan

    #Run first round of PLINK through primary study
    #Find files with chr in name in primary study dir
    primary_chr_files = [os.path.join(primary_dir, f) for f in os.listdir(primary_dir) if os.path.isfile(os.path.join(primary_dir, f))]
    #Get rid of suffixes
    primary_chr_files = [os.path.splitext(f)[0] for f in primary_chr_files]
    primary_chr_files = list(set(primary_chr_files))
    #Do the same for second study if available
    if second_study:
        secondary_chr_files = [os.path.join(secondary_dir, f) for f in os.listdir(secondary_dir) if os.path.isfile(os.path.join(secondary_dir, f))]
        secondary_chr_files = [os.path.splitext(f)[0] for f in secondary_chr_files]
        secondary_chr_files = list(set(secondary_chr_files))




    #Get unique values to avoid PLINK redundancy
    primary_variants = list(set(variants))
    all_variants = primary_variants #backup

    #Get unique temp filename, write variants for PLINK to use
    tf = tempfile.NamedTemporaryFile()
    tf_name = output_dir+tf.name.split('/')[2]
    with open(tf_name, mode='wt') as variant_file:
        variant_file.write('\n'.join(primary_variants))

    #Figure out which variant files we need to query
    var_chrs = list(GWAS_variant_table['chr'].drop_duplicates())

    #Run PLINK on each chr for which we have variants
    prim_lines = []
    sec_lines = []
    for chr in var_chrs:
        try:
            primary_query_file = [s for s in primary_chr_files if s.endswith('chr'+str(chr))][0]
        except: 
            print("No files to query for chromosome "+chr+", skipping.")
            continue
        temp_out=tf_name+'chr'+str(chr)
        call_plink(query_file=primary_query_file,
                   subset_file=primary_subset, snp_file=tf_name,
                   ld_kb=ldkb, ld_r2=primary_ldr, out_name=temp_out, memory=str(memory))
        try:
            #TODO
            add_df = pd.read_csv(temp_out+'.ld', header=0, sep='\s+')
            print("Len of chr"+str(chr)+" proxy variants is "+str(len(add_df)))
            prim_lines.append(add_df)
        except:
            print('No LD variants for our GWAS variants recovered from chr'+str(chr)+'.')
        #recycle temp file; PLINK doesn't overwrite, it uses different names
        clean_plink(tf_name)

        #If second study is presented, run PLINK now
        if second_study:
            secondary_query_file = [s for s in secondary_chr_files if s.endswith('chr'+str(chr))][0]
            #Write variants in format that second PLINK files can find
            ntf = tempfile.NamedTemporaryFile()
            ntf_name = output_dir + ntf.name.split('/')[2]
            with open(ntf_name, mode='wt') as variant_file:
                variant_file.write('\n'.join(secondary_variants))

            n_temp_out = ntf_name+'chr'+str(chr)
            call_plink(query_file=secondary_query_file,
                       subset_file=secondary_subset, snp_file=ntf_name,
                       ld_kb=ldkb, ld_r2=secondary_ldr, out_name=n_temp_out, memory=str(memory))
            try:
                sec_lines.append(pd.read_csv(n_temp_out+'.ld', header=0, sep='\s+'))
            except:
                print("No LD variants for our GWAS variants recovered from chr"+str(chr)+" in "+secondary_prefix+".")
            ntf.close()
            deep_clean_plink(ntf_name)

    #cleanup
    print("Cleaning up temp files...")
    tf.close()
    deep_clean_plink(tf_name)
    print("Done cleaning up.")

    #cleanup



    #Turn into one df
    if len(prim_lines) > 1:
        prim_lines=pd.concat(prim_lines)
    elif len(prim_lines) == 1:
        prim_lines = prim_lines[0]
    #TODO
    print("Len of primary lines is "+str(len(prim_lines)))

#else:
#    raise Exception("None of your variants appear in your primary study. Please re-check inputs.")

    if second_study:
        if len(sec_lines) > 1:
            sec_lines = pd.concat(sec_lines)
        elif len(sec_lines) == 1:
            sec_lines = sec_lines[0]


    if ((second_study) & (len(sec_lines) >=1)): 
        #Find variants that don't have proxies in primary VCFs
        if len(prim_lines) > 0:
            primary_series = pd.Series(primary_variants)
            need_proxies = list(set(primary_series.loc[~(primary_series.isin(prim_lines['SNP_A']))]))
            #TODO
            print("Number of variants that need proxies "+str(len(need_proxies)))
        else: 
            need_proxies = primary_variants #If none of the variants were found in primary study
    
        if len(need_proxies) > 0:
            GWAS_variant_table.loc[GWAS_variant_table['GWAS_variant'].isin(need_proxies), secondary_prefix+'_proxies_needed'] = True
            GWAS_variant_table.loc[GWAS_variant_table['GWAS_variant'].isin(need_proxies), secondary_prefix+'_proxies_found'] = False
            print("Some of the given variants were not found in the primary VCF. Trying to find viable proxies from the secondary VCF.") 
            need_proxies = [convert_variant_format(x, example_variant=example_secondary_variant) for x in need_proxies]
            sec_candidates = list(sec_lines.loc[sec_lines['SNP_A'].isin(need_proxies), 'SNP_B'].drop_duplicates().copy())

            #Write temporary file, then query for it 
            sec_query_first = [convert_variant_format(x, example_variant=example_primary_variant) for x in sec_candidates]
        
            #Get unique temp filename, write variants for PLINK to use
            tf = tempfile.NamedTemporaryFile()
            tf_name = output_dir+tf.name.split('/')[2]
            with open(tf_name, mode='wt') as variant_file:
                variant_file.write('\n'.join(sec_query_first))

            #Figure out which variant files we need to query
            var_chrs = list(sec_lines['CHR_A'].drop_duplicates())

            #Run PLINK on each chr for which we have variants
            prim_ext_lines = []
            for chr in var_chrs:
                try:
                    primary_query_file = [s for s in primary_chr_files if s.endswith('chr'+str(chr))][0]
                except:
                    print("No files to query for chromosome "+str(chr)+", skipping.")
                    continue
                temp_out=tf_name+'chr'+str(chr)
                call_plink(query_file=primary_query_file,
                           subset_file=primary_subset, snp_file=tf_name,
                           ld_kb=ldkb, ld_r2=primary_ldr, out_name=temp_out, memory=str(memory))
                try:
                    prim_ext_lines.append(pd.read_csv(temp_out+'.ld', header=0, sep='\s+'))
                except:
                    print('No LD variants for our GWAS variants recovered from chr'+str(chr)+'.')
                #recycle temp file; PLINK doesn't overwrite, it uses different names
            deep_clean_plink(tf_name)

            #Turn into one df
            if len(prim_ext_lines) > 1:
                prim_ext_lines=pd.concat(prim_ext_lines)
            elif len(prim_ext_lines) == 1:
                prim_ext_lines = prim_ext_lines[0]
            else:  
                print("No suitable proxies were found.")

            #TODO
            print("Number of proxies is "+str(len(prim_ext_lines)))
   
            print('Cleaning up temp files...')
            tf.close()
            deep_clean_plink(tf_name)
            print("Done cleaning up.")

            #Find variants for which we found proxies
            prim_ext_lines = prim_ext_lines.loc[prim_ext_lines.SNP_A != prim_ext_lines.SNP_B].copy()
            subs_found = prim_ext_lines.SNP_A.drop_duplicates()
            GWAS_variant_table.loc[GWAS_variant_table['GWAS_variant'].isin(subs_found), secondary_prefix+'_proxies_found'] = True

    ### Proxy Searching completed ###
    print("Proxy queries completed.")

    #['GWAS_variant', 'chr', 'pos_bp', 'plink2_proxies_needed', 'plink2_proxies_found']
    #['CHR_A', 'BP_A', 'SNP_A', 'MAF_A', 'CHR_B', 'BP_B', 'SNP_B', 'MAF_B','R2', 'DP']
    #prim_lines ##variants with LD in GTEx, SNP A in GWAS_variant_table.GWAS_variant
    #sec_lines ##data with LD variant information in 1000 genomes
    #prim_ext_lines ##LD variants in GTEx for 1000 genomes proxies

    print("Determining where proxies are needed...")

    #Convert variants from secondary format to primary format
    sec_lines['SNP_A'] = sec_lines.SNP_A.apply(lambda x: convert_variant_format(x, example_variant=example_primary_variant))
    sec_lines['SNP_B'] = sec_lines.SNP_B.apply(lambda x: convert_variant_format(x, example_variant=example_primary_variant))
    sec_lines = sec_lines.loc[sec_lines.SNP_A != sec_lines.SNP_B].copy()

    #grab proxies
    new_table = GWAS_variant_table.merge(prim_lines[['SNP_A', 'SNP_B', 'BP_B',  'MAF_A', 'MAF_B', 'R2', 'DP']], how='left', left_on='GWAS_variant', right_on='SNP_A').copy()
    print("Len of primary lines table after merge "+str(len(new_table)))
    new_table.drop('SNP_A', axis=1,inplace=True)
    new_table.columns = ['GWAS_variant', 'chr', 'pos_bp', new_table.columns[3], new_table.columns[4], 'LD_variant', 'LD_pos_bp',
                         'GWAS_variant_AF_in_'+primary_prefix, 'LD_variant_AF_in_'+primary_prefix, 'LD_R2', 'LD_DP']

    #Figure out variants that found proxies
    fp = new_table.loc[new_table[new_table.columns[4]] == True, 'GWAS_variant'].copy()
    sl = sec_lines.loc[sec_lines.SNP_A.isin(list(fp))][['SNP_A','CHR_A', 'BP_A','SNP_B', 'BP_B', 'MAF_A','MAF_B','R2','DP']].copy() #proxies
    sl = sl.loc[sl.R2.astype(float) >= secondary_ldr].copy()
    sl.columns = ['GWAS_variant', 'chr', 'pos_bp', 'proxy_variant', 'proxy_pos_bp', 
                  'GWAS_AF_in_'+secondary_prefix,'proxy_AF_in_'+secondary_prefix,'proxy_R2','proxy_DP']

    try:
        print("Secondary lines for which proxies were found "+str(len(sl)))
        print("Lines by which primary lines will be extended "+str(len(prim_ext_lines)))
        slm = sl.merge(prim_ext_lines[['SNP_A','SNP_B', 'BP_B', 'MAF_A','MAF_B','R2','DP']], how='right', left_on='proxy_variant', right_on='SNP_A')
    except:
        print("No proxies for variants that needed proxies were found.")
        slm=sl.copy()
        slm['SNP_A'] = np.nan
        slm['SNP_B'] = "No proxies found."
        slm['BP_B'] = "No proxies found."
        slm['MAF_A'] = "No proxies found."
        slm['MAF_B'] = "No proxies found."
        slm['R2'] = "No proxies found."
        slm['DP'] = "No proxies found."
    print("Lines by which primary lines will be extended "  + str(len(slm)))
    slm.columns = ['GWAS_variant','chr', 'pos_bp','proxy_variant', 'proxy_pos_bp', slm.columns[5],
                   slm.columns[6], 'proxy_R2', 'proxy_DP', 'SNP_A', 'LD_variant', 'LD_pos_bp',
                   'proxy_AF_in_'+primary_prefix, 'LD_variant_AF_in_'+primary_prefix, 'LD_R2', 'LD_DP']
    slm.drop('SNP_A', axis=1, inplace=True)

    #Create place-holder columns for all values
    slm[secondary_prefix+'_proxies_needed'] = True
    slm[secondary_prefix+'_proxies_found'] = True
    slm['GWAS_variant_AF_in_'+primary_prefix] = 0

    new_table['proxy_variant'] = 'None'
    new_table['proxy_pos_bp'] = 'No_proxies_taken'
    new_table['GWAS_AF_in_'+secondary_prefix] = 'Not_checked'
    new_table['proxy_AF_in_'+secondary_prefix] = 'No_proxies_taken'
    new_table['proxy_R2'] = 'No_proxies_taken'
    new_table['proxy_DP'] = 'No_proxies_taken'
    new_table['proxy_AF_in_'+primary_prefix] = 'No_proxies_taken'
    new_table = new_table[slm.columns]

    proxy_table = pd.concat([new_table, slm])
    print("Final proxy table length "+str(len(proxy_table)))
    proxy_table = proxy_table.drop_duplicates()

    #Determine "tested" variant
    proxy_table['tested_variant'] = proxy_table['GWAS_variant']
    proxy_table.loc[((proxy_table[secondary_prefix+'_proxies_needed'] ==True) & (proxy_table[secondary_prefix+'_proxies_found']==True)), 
                    'tested_variant'] = proxy_table.loc[((proxy_table[secondary_prefix+'_proxies_needed'] ==True) & (proxy_table[secondary_prefix+'_proxies_found']==True)), 
                    'proxy_variant']

    proxy_table['tested_variant_pos_bp'] = proxy_table['pos_bp']
    proxy_table.loc[((proxy_table[secondary_prefix+'_proxies_needed'] ==True) & (proxy_table[secondary_prefix+'_proxies_found']==True)), 
                    'tested_variant_pos_bp'] = proxy_table.loc[((proxy_table[secondary_prefix+'_proxies_needed'] ==True) & (proxy_table[secondary_prefix+'_proxies_found']==True)), 
                    'proxy_pos_bp']

    proxy_table = proxy_table[['chr', 'GWAS_variant', proxy_table.columns[17], 'tested_variant', 'tested_variant_pos_bp',
                               'proxy_DP', 'proxy_R2', proxy_table.columns[5], proxy_table.columns[6], 
                               'LD_variant', 'LD_pos_bp', proxy_table.columns[11], proxy_table.columns[12],
                               'LD_R2', 'LD_DP', proxy_table.columns[15], proxy_table.columns[16]]].copy()

    #Drop values where no proxies were filled
    proxy_table = proxy_table.loc[~((proxy_table[proxy_table.columns[-2]] == True) & (proxy_table[proxy_table.columns[-1]]==True) & 
                                    (proxy_table[proxy_table.columns[2]].isnull()))].copy()

    print("Double checeking final proxy table length "+str(len(proxy_table)))

    print("Adding proxy pairing ID...")
    #Add a unique key
    proxy_groups = proxy_table.groupby('GWAS_variant')
    new_proxy_table = []
    for group in proxy_groups.groups:
        proxy_table = proxy_groups.get_group(group).copy()
        ranger = [str(x) for x in range(0, len(proxy_table))]
        proxy_table['proxy_pairing_ID'] = proxy_table['GWAS_variant'] +'_'+ ranger
        new_proxy_table.append(proxy_table)

    print("Writing proxy table.")
    proxy_table = pd.concat(new_proxy_table)

    print("Triple checking proxy table length "+str(len(proxy_table)))
    print("Length of proxy table is "+str(len(proxy_table)))
    proxy_table.to_csv(output_dir+output_prefix+'_proxy_table.tsv', sep='\t', index=False)
    traits.to_csv(output_dir+output_prefix+'_trait_table.tsv', sep='\t', index=False)
    print("Done.")        

    if variants_only:
        print("User has indicated that this run is intended only to retrieve LD information. Ending now.")
        exit()

else: 
    print("Development build selected. Loading proxy file...")
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    proxy_table = dd.read_csv(output_dir+output_prefix+'_proxy_table.tsv', sep='\t', low_memory=False,
                              dtype={'chr':str})
    traits_table = dd.read_csv(output_dir+output_prefix+'_trait_table.tsv', sep='\t', low_memory=False,
                              dtype={'chr':str})
    print("Done.")

#Create gene dict for gene names
gtf_dict = parse_gtf(gtf_file)

#Build best eGene dict
print("Building best eGene dict for annotation...")
best_eGene_dict = build_best_eGene_dict(tissue_file)
print("Building best sGene dict for annotation...")
best_sGene_dict = build_best_eGene_dict(tissue_file, type='sQTL')
print("Done.")

#print("Setting up LD variants dict for eQTL/sQTL parsing...")
mapping_df = proxy_table[['LD_variant', 'GWAS_variant', 'proxy_pairing_ID', 'chr']].copy()
#chrs = set(mapping_df.chr)
#LD_variants_dict = {chr:set(mapping_df.loc[mapping_df.chr == chr, 'LD_variant']) for chr in chrs}
print("Done.")
#
#Parse through each tissue 
eQTL_table = []
eQTL_tissue_table = []
sQTL_table = []
sQTL_tissue_table = []
#eQTL_DAPG_table = []
#sQTL_DAPG_table = []

print("Reading in eQTL/sQTL files")
with g_open(tissue_file) as f:
    for i, line in enumerate(f):
        row = line.strip().split('\t')
        tissue = row[0]
        print(tissue)
        esig_file = row[1]
        eqtl_df, eqtl_tiss_df = parse_eqtl_file(esig_file, mapping_df, #LD_variants_dict,
                                gtf_dict, tissue, best_eGene_dict)
        
        eQTL_table.append(eqtl_df)
        eQTL_tissue_table.append(eqtl_tiss_df)

        if len(eQTL_table) > 1:
            eQTL_table = [dd.concat(eQTL_table)]
            eQTL_tissue_table = [dd.concat(eQTL_tissue_table)]
        try:
            ssig_file = row[3]    
            sqtl_df, sqtl_tiss_df = parse_eqtl_file(ssig_file, mapping_df, #LD_variants_dict,
                                    gtf_dict, tissue, best_sGene_dict, type='sQTL')

            sQTL_table.append(sqtl_df)
            sQTL_tissue_table.append(sqtl_tiss_df)


            if len(sQTL_table) > 1:
                sQTL_table = [dd.concat(sQTL_table).drop_duplicates()]
                sQTL_tissue_table = [dd.concat(sQTL_tissue_table).drop_duplicates()]
        except:
            print("No sQTLs for "+row[0]+". Continuing.")


print("Done processing eQTL/sQTL files.")

eQTL_table = dd.concat(eQTL_table).drop_duplicates()
eQTL_tissue_table = dd.concat(eQTL_tissue_table).drop_duplicates()
sQTL_table = dd.concat(sQTL_table).drop_duplicates()
sQTL_tissue_table = dd.concat(sQTL_tissue_table).drop_duplicates()
#eQTL_table.reset_index(drop=True)
#sQTL_table.reset_index(drop=True)

print("Writing eQTL and sQTL tables...")
eQTL_table.to_csv(output_dir+output_prefix+'_eQTL_table.tsv', sep='\t', index=False, single_file=True)
print('1')
eQTL_tissue_table.to_csv(output_dir+output_prefix+'_eQTL_tissue_table.tsv', sep='\t', index=False, single_file=True)
print('2')
sQTL_table.to_csv(output_dir+output_prefix+'_sQTL_table.tsv', sep='\t', index=False, single_file=True)
print('3')
sQTL_tissue_table.to_csv(output_dir+output_prefix+'_sQTL_tissue_table.tsv', sep='\t', index=False, single_file=True)
print("Done.")



