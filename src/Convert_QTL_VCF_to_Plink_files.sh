#!/bin/bash

# Script converts VCF files split by chromosome to Plink files per chromosome.
# This example is run on the GTEx v8 WGS analysis freeze VCF available via dbGaP (accession no. phs000424.v8)

in_dir='VCF_chrs/'
for chr in {{1..22},X}
do
    file_name="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.chr$chr"
    file="GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.chr$chr.vcf.gz"
    file=$in_dir$file
    plink --vcf $file --memory 20480 --threads 1 --out $file_name
done

