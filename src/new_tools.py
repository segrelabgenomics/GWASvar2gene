# Author: John Rouhana
# Segre lab, Massachusetts Eye and Ear, Harvard Medical School
# Date: 2020

from collections import defaultdict
import argparse
import numpy as np
import pandas as pd
import re
import os
from argparse import ArgumentParser as eparser
import subprocess
import glob
import gzip

def g_open(file):
    """
    Determines if a file is gzipped
    opens as appropriate
    """
    if file.endswith('gz'):
        return(gzip.open(file, 'rt'))
    else:
        return(open(file, 'r'))

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def ERRORPARSING(args):
    #Verify secondary plink files given if subset file received
    if (not (args.secondary_subset is None)) and (args.secondary_plink_directory is None):
        raise Exception("-secondary_subset requires -secondary_plink_directory.")

    #Verify host exists if writing database
    if args.use_db and (args.db_host is None):
        raise Exception( "-mySQL write requested but host not given.")

    if not args.variants_only:
        #Verify inputs needed given
        if args.all_QTL_analysis_variants is None:
            raise Exception("-all_QTL_analysis_variants required when variants_only is False.")

        if args.GENCODE is None:
            raise Exception("-GENCODE required when variants_only is False.")

        if args.tissue_file is None:
            raise Exception("-tissue_file required when variants_only is False.")

        #Verify GENCODE file exists
        if not os.path.isfile(args.GENCODE):
            raise Exception("-GENCODE should point to a GENCODE gtf file. Current file does not exist.")

        #Verify QTL variants exists
        if not os.path.isfile(args.all_QTL_analysis_variants):
            raise Exception("-all_QTL_analysis_variants should point to a file containing all variants used in QTL analysis. Current file does not exist.")

        #Verify eQTL_dir exists
        if not os.path.isfile(args.tissue_file):
            raise Exception("-tissue_file should point to a tab-delimited file associating tissues with eQTL/sQTL files. Selected file does not exist.")

    #Verify secondary file exists
    if args.secondary_plink_directory is not None:
        if not os.path.isdir(args.secondary_plink_directory):
            raise Exception("-secondary_plink_directory should point to a folder containing PLINK files from a study. Selected directory does not exist.")
        secondary_dir = args.secondary_plink_directory
        secondary_files = [os.path.join(secondary_dir, f) for f in os.listdir(secondary_dir) if os.path.isfile(os.path.join(secondary_dir, f))]
        try:
            bimtwo = [f for f in secondary_files if f.endswith('bim')][0]
        except:
            raise Exception("-secondary_plink_directory does not appear to point to a directory containing PLINK files for a second study.")

def call_plink(query_file, subset_file, snp_file, ld_kb, ld_r2, out_name,
               memory='None', threads='1', ld_window='99999', silent=True):
    """
    Function to call PLINK from PYTHON
    must import subprocess

    chr_dir: directory where PLINK files are located
    query_file: Name of the file being queried 
    subset_file: File containing subset of individuals  to be kept from PLINK files
    snp_file: name of file with variants of interest
    ld_kb: kilobase window range to look for LD (suggest: 5000)
    ld_r2: r2 cutoff to consider two variants to be in LD (suggest: 0.5)
    ld_window: number of variants allowed between variants a and b to still be considered in LD


    memory: how much memory (in MB) to partition
    threads: how many threads to allot for PLINK
    
    """
    new_list = ['plink', '--bfile', query_file,
                '--r2', 'dprime', 'with-freqs', '--ld-snp-list', 
                snp_file, '--ld-window-kb', str(ld_kb),
                '--ld-window-r2', str(ld_r2), '--ld-window', str(ld_window)]
    #If subset file given, include that parameter
    if subset_file is not None:
        new_list = new_list + ['--keep', subset_file]
    #If memory specified, pass parameter
    if memory is not None:
        new_list = new_list + ['--memory', str(memory)]  
    #same for threads 
    if threads is not None:
        new_list = new_list + ['--threads', str(threads)]
    if silent == True:
        new_list = new_list + ['--silent']
    #add out parameter at end
    new_list = new_list + ['--out', out_name]
    subprocess.call(new_list)



def convert_variant_format(in_variant, example_variant=None, ref_format=None):
    """
    Function to convert variant between different formats.
    Either an example variant or the desired output format
    are required to run. Available output formats are
    GTEx_v7, GTEx_v8, G1000 (for Thousand Genomes).

    If both example_variant and format are provided, 
    example_variant is used.

    Function can read the above three listed variant formats
    for examples.
    """
    gtex_v7_regex = re.compile('[0-9,X,Y]+?_[0-9]*?_[ATCG]*?_[ATCG]*?_b37$')
    gtex_v8_regex = re.compile('chr[0-9,X,Y]+?_[0-9]*?_[ATCG]*?_[ATCG]*?_b38$')
    G1000_regex = re.compile('[0-9,X,Y]+?:[0-9]*?:[ATCG]*?:[ATCG]*$')
    G1000_alt_regex = re.compile('chr[0-9,X,Y]+?:[0-9]*?:[ATCG]*?:[ATCG]*$')
    alt_regex = re.compile('[0-9,X,Y]+?_[0-9]*?_[ATCG]*?_[ATCG]*$')
    if example_variant is not None:
        if gtex_v7_regex.match(example_variant) is not None:
            ref_format = 'gtex_v7'
            if gtex_v7_regex.match(in_variant) is not None:
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                raise ValueError("Reference variant and conversion variant "
                                 "appear to be from difference genome builds.")
                sys.exit()
            elif G1000_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':','_')+'_b37'
                return in_variant
            elif G1000_alt_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':', '_').replace('chr', '')+'_b37'
            elif alt_regex.match(in_variant) is not None:
                in_variant = in_variant+'_b37'
                return(in_variant)
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif gtex_v8_regex.match(example_variant) is not None:
            ref_format = 'gtex_v8'
            if gtex_v7_regex.match(in_variant) is not None:
                raise ValueError("Reference variant and conversion variant "
                                 "appear to be from difference genome builds.")
                sys.exit()
            elif gtex_v8_regex.match(in_variant) is not None:
                return in_variant
            elif G1000_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant.replace(':','_')+'_b38'
                return in_variant
            elif G1000_alt_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':', '_')+'_b38'
                return(in_variant)
            elif alt_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant+'_b38'
                return(in_variant)
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif G1000_regex.match(example_variant) is not None:
            ref_format = 'G1000'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':')[:-4]
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':')[3:-4]
                return in_variant
            elif G1000_regex.match(in_variant):
                return in_variant
            elif G1000_alt_regex.match(in_variant):
                in_variant = in_variant.replace('chr', '')
                return in_variant
            elif alt_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':')
                return in_variant
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return in_variant

        elif G1000_alt_regex.match(example_variant) is not None:
            ref_format = 'G1000_alt'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant.replace('_', ':')[:-4]
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':')[:-4]
                return in_variant
            elif G1000_regex.match(in_variant):
                in_variant = 'chr'+in_variant
                return in_variant
            elif G1000_alt_regex.match(in_variant):
                return in_variant
            elif alt_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant.replace('_', ':')
                return in_variant
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif alt_regex.match(example_variant) is not None:
            ref_format = 'alt'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = in_variant[:-4]
                return(variant)
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant[3:-4]
                return(variant)
            elif G1000_regex.match(in_variant):
                in_variant = in_variant.replace(':', '_')
                return(in_variant)
            elif G1000_alt_regex.match(in_variant):
                in_variant = in_variant.replace(':', '_').replace('chr', '')
            elif alt_regex.match(in_variant) is not None:
                return(in_variant)
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)


        else: # If example is not recognized
            raise ValueError("Example variant format not recognized.")
            sys.exit()
    elif ref_format is not None: #If example is not provided
        if ref_format == 'GTEx_v7':
            ref_format = 'gtex_v7'
            if gtex_v7_regex.match(in_variant) is not None:
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                raise ValueError("Reference variant and conversion variant "
                                 "appear to be from different genome builds.")
                sys.exit()
            elif G1000_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':','_')+'_b37'
                return in_variant
            elif G1000_alt_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':', '_').replace('chr', '')+'_b37'
            elif alt_regex.match(in_variant) is not None:
                in_variant = in_variant+'_b37'
                return in_variant
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif ref_format == 'GTEx_v8':
            ref_format = 'gtex_v8'
            if gtex_v7_regex.match(in_variant) is not None:
                raise ValueError("Reference variant and conversion variant "
                                 "appear to be from difference genome builds.")
                sys.exit()
            elif gtex_v8_regex.match(in_variant) is not None:
                return in_variant
            elif G1000_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant.replace(':','_')+'_b38'
                return in_variant
            elif G1000_alt_regex.match(in_variant) is not None:
                in_variant = in_variant.replace(':', '_')+'_b38'
                return in_variant
            elif alt_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant+'_b38'
                return in_variant
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif ref_format == 'G1000':
            ref_format = 'G1000'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':').replace(":b37", '')
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':').replace(':b38', '').replace('chr', '')
                return in_variant
            elif G1000_regex.match(in_variant):
                return in_variant
            elif G1000_alt_regex.match(in_variant):
                in_variant = in_variant.replace('chr', '')
                return(in_variant)
            elif alt_regex.match(in_variant):
                in_variant = in_variant.replace('_', ':')
                return(in_variant)
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif ref_format == 'G1000_alt':
            ref_format = 'G1000_alt'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = 'chr'+in_variant.replace('_', ':').replace(":b37", '')
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant.replace('_', ':').replace(':b38', '')
                return in_variant
            elif G1000_regex.match(in_variant):
                in_variant = 'chr'+in_variant
                return in_variant
            elif G1000_alt_regex.match(in_variant):
                return in_variant
            elif alt_regex.match(in_variant):
                in_variant = 'chr'+in_variant.replace('_', ':')
                return in_variant
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        elif ref_format == 'alt':
            ref_format = 'alt'
            if gtex_v7_regex.match(in_variant) is not None:
                in_variant = in_variant[:-4]
                return in_variant
            elif gtex_v8_regex.match(in_variant) is not None:
                in_variant = in_variant[3:-4]
                return in_variant
            elif G1000_regex.match(in_variant):
                in_variant = in_variant.replace(':', '_')
                return in_variant
            elif G1000_alt_regex.match(in_variant):
                in_variant = in_variant.replace(':', '_').replace('chr', '')
            elif alt_regex.match(in_variant):
                return(in_variant)
            else:
                print("Variant format for "+in_variant+" not recognized. Returning as-is.")
                return(in_variant)

        else:
            raise ValueError("Reference format not recognized. Options are "
                             "GTEx_v7, GTEx_v8, G1000 (for Thousand Genomes),"
                             "and alt. See documentation. Case sensitive.")
            sys.exit()

    else: #If neither example nor in_format provided
        raise ValueError("No example_variant or desired format provided "
              "for the conversion of variant format.")
        sys.exit()



#"""
#GTF parser written and open-sourced by 
#Kamil Slowikowski. Modified for readability.
#See more here:
#https://gist.github.com/slowkow/8101481#file-gtf-py-L52
#"""
GTF_HEADER  = ['seqname', 'source', 'feature', 'start', 'end', 'score',
               'strand', 'frame']
R_SEMICOLON = re.compile(r'\s*;\s*')
R_COMMA     = re.compile(r'\s*,\s*')
R_KEYVALUE  = re.compile(r'(\s+|\s*=\s*)')


def dataframe(filename):
    """
    Inputs:
    filename= name of gtf file, optionally gzipped
    
    Outputs:
    pandas dataframe of gzipped file
    """
    # Each column is a list stored as a value in this dict.
    result = defaultdict(list)

    for i, line in enumerate(lines(filename)):
        for key in line.keys():
            # This key has not been seen yet, so set it to None for all
            # previous lines.
            if key not in result:
                result[key] = [None] * i

        # Ensure this row has some value for each column.
        for key in result.keys():
            result[key].append(line.get(key, None))

    return pd.DataFrame(result)


def lines(filename):
    """
    Inputs: 
    filename= name of gtf file, optionally gzipped
    Outputs:
    generator yielding line by line, excluding header
    """
    if filename.endswith('.gz'):
        fn_open = gzip.open
    else:
        fn_open = open
    with fn_open(filename, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            else:
                yield parse(line)

def parse(line):
    """
    Inputs:
    line: A single line from a GTF file

    Outputs:
    A Dict pointing header to value
    """
    result = {}

    #check if bytes loaded
    fields = line.rstrip().split('\t')
    for i, col in enumerate(GTF_HEADER):
        result[col] = _get_value(fields[i])

    # INFO field consists of "key1=value;key2=value;...".
    infos = [x for x in re.split(R_SEMICOLON, fields[8]) if x.strip()]

    for i, info in enumerate(infos, 1):
        # It should be key="value".
        try:
            key, _, value = re.split(R_KEYVALUE, info, 1)
        # But sometimes it is just "value".
        except ValueError:
            key = 'INFO{}'.format(i)
            value = info
        # Ignore the field if there is no value.
        if value:
            result[key] = _get_value(value)

    return result


def _get_value(value):
    if not value:
        return None

    # Strip double and single quotes.
    value = value.strip('"\'')

    # Return a list if the value has a comma.
    if ',' in value:
        value = re.split(R_COMMA, value)
    # These values are equivalent to None.
    elif value in ['', '.', 'NA']:
        return None

    return value
#"""
#End GTF parser
#"""

def clean_plink(tf_name):
    """
    Deletes PLINK outputs
    """
    for f1 in glob.glob(tf_name+'chr*'):
        try:
            os.remove(f1)
        except:
            continue

def deep_clean_plink(tf_name):
    """
    Deletes all temporary files under a temporary ID
    """
    for f1 in glob.glob(tf_name+'*'):
        try:
            os.remove(f1)
        except:
            continue


###
#GENCODE PARSER FOR DICT
###
#Object containing basic information for each gene
class gene:
    def __init__(self, chr, gene_name, gene_id, gene_start, gene_finish):
        self.gene_name = gene_name
        self.gene_id = gene_id
        self.start = int(gene_start)
        self.end = int(gene_finish)
        self.chr = chr
        self.gene_range = range(self.start, self.end)

#Object for intervals
class interval:
    def __init__(self, chr, interval_start, interval_end):
        self.chr = chr
        self.start = int(interval_start)
        self.end = int(interval_end)
        self.interval_range = range(self.start, self.end)

#Parses gtf file, returns 'gene' objects
def parse_gtf(gtf_file, compression=None):
    """
    gtf_file: path to a GTF file, compressed or not
    """
    gene_dict = {}
    if compression == None:
        gtf = open(gtf_file)
    elif compression == 'gz':
        gtf = gzip.open(gtf_file, 'rb')
    else:
        raise Exception("Compression of gtf file not recognized. Please used "
                        ".gz or no compression.")
    with gtf:
        for row in gtf:
            if isinstance(row, bytes):
                row = row.decode('utf-8')
            row = row.strip().split('\t')
            if row[0][0]=='#' or row[2]!='gene': continue
            attr = dict([i.split() for i in row[8].replace('"','').split(';') if i!=''])
            gene_dict[attr['gene_id'].split('.')[0]]=(gene(row[0],
                      attr['gene_name'],attr['gene_id'].split('.')[0], row[3], row[4]))
    return(gene_dict)

