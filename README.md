## GWASvar2gene

GWASvar2gene is a tool that annotates a set of GWAS variants for complex diseases or trait or any variants of interest, with putative causal genes and regulatory mechanisms, based on eQTLs and sQTLs from GTEx v8 or other molecular QTL studies that are in linkage disequilibrium with each variant.

**Author:** John Rouhana, written in the Segre lab under the guidance of Ayellet Segre, Massachusetts Eye and Ear, Harvard Medical School, Boston, MA

Contact: ayellet_segre@meei.harvard.edu

## Repository structure

`src`: directory contains scripts for the tool 

`data`: directory contains input files required to run GWASvar2gene using cis-eQTLs and cis-sQTLs from GTEx v8 

`example`: directory contains sample data and paths to source files for running GWASvar2gene on a sample list of GWAS variants

`output`: directory contains sample output files from the different steps of the GWASvar2gene package.


## Dependencies

GWASvar2gene was written in Python (at least 3.8) and requires the following libraries and modules:

• pandas 2.5.0 

• numpy 1.16.2 

• scipy 0.17.1 

• matplotlib 2.2.2 

• ggplot2 3.2.1 

• dplyr 0.8.3 

• readr 1.3.1 

• cowplot 1.0.0 

• miniconda/4.7.12

• PLINK/1.9b-6.10


## Steps for running GWASvar2gene

### STEP 1. Prepare table with list of GWAS variants or variants of interest, and run the following script to add variant IDs according to GTEx format (CHR_POS_REF_ALT_GENOMEBUILD) if this column does not exist. 

GWAS variants for complex diseases or traits of interest can be taken from the NHGRI-NBI GWAS catalog (https://www.ebi.ac.uk/gwas/docs/file-downloads) or Open Targets Genetics (https://genetics-docs.opentargets.org/data-access/data-download).

Minimum required columns: 'CHR_ID', 'CHR_POS', 'STRONGEST SNP-RISK ALLELE'

#### Example shell script to run step 1:  
`/src/run_add_GTEx_variant_IDs_to_GWAS_variant_list.sh`

Runs: ./add_GTEx_variant_ID_to_GWAS_variant_list.py

### INPUT FILES:
Three input files need to be defined on the top of `add_GTEx_variant_ID_to_GWAS_variant_list.py` 

1. List of variants tested in QTL database. Required columns: [1]CHROM	[2]POS  [3]REF  [4]ALT  [5]ID
(i.e., chr number, position, REF allele, ALT allele, variant ID (CHR_POS_REF_ALT_GENOMEBUILD) )

We provide table of variants analyzed in GTEx release v8 for eQTLs and sQTLs:
GTEx = '/data/GTEx_v8_HG38_all_variants.tsv.gz'

2.  Table of variants that GWASvar2gene will annnotate with eQTLs and sQTLs, e.g., variants downloaded from GWAS catalog (https://www.ebi.ac.uk/gwas/docs/file-downloads).
Required columns in input GWAS variant table: 'CHR_ID', 'CHR_POS', 'STRONGEST SNP-RISK ALLELE' (e.g., rs1925953-T)

Example file:
GWAS = '/example/gwas_catalog_sample_file.tsv'

3. Output file name:
out_file = '/output/out_GWAS_variant_IDs_added.tsv'

The GTEx formatted variant IDs are in column 'GWAS_variant'.

### STEP 2. Run GWASvar2gene to map genes to GWAS variants based on eQTLs and sQTLs

#### Example shell script to run step 2:
`./src/run_GWASVar2Gene.sh`

Runs: GWASVar2Gene.py

### Input arguments:

`--GWAS_variants`: Table with list of GWAS variants to be annotated. Default table is the GWAS catalog table. Minimum required columns: 'CHR_ID', 'CHR_POS', 'STRONGEST SNP-RISK ALLELE'

`--primary_plink_directory`: Directory with genotype files of the primary reference panel in Plink format. Best to use genotype study that is most compatible with QTL study. As default we use the GTEx WGS plink files from release v8 that requires dbGaP access (e.g., GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.chr*.bed, *.bim, *.fam) or one can use another reference panel that matches the e/sQTL study used.

The script `/src/Convert_QTL_VCF_to_Plink_files.sh` can be used to convert the QTL study VCF files to Plink files by chromosome.

`-primary_plink_prefix`: Name of primary genotype reference panel (Default: GTEx)

`-primary_ldr`: r2 cutoff for determining if an e/sQTL is in linkage disequillibrium (LD) with a GWAS variant using the primary reference panel Plink files (Default: 0.5; one can use 0.8 to be more stringent).

`-primary_subset`: Sample IDs of samples in primary reference panel to be considered for the LD calculation (useful if you want to analyze only a subset of samples from a given ancestral background). This tab-delimied file with no header should contain two columns with the sample IDs (column 1: family ID if available, otherwise sample ID; column 2: sample ID; both columns can be sample IDs (identical))

`-secondary_plink_directory`: Directory with genotype files in Plink format for secondary reference panel. This reference panel is used if GWAS variant is not found in primary reference panel.

`-secondary_plink_prefix`: Name of secondary genotype reference panel (Default: 1000_genomes)

'-secondary_ldr': r2 cutoff for finding LD proxy variants in secondary panel for GWAS variants with no LD variants in the primary PLINK file set (Default: 0.8).

`-secondary_subset`: Sample IDs of samples in secondary reference panel to be considered in the LD calculation (See `-primary_subset` for file format)

`-output_dir`: Name of output directory

`-output_prefix`: Prefix of output file name

`-GENCODE`: Gencode GTF file with gene annotations. As default we use gencode.v26.annotation.gtf used in GTEx v8, which can be downloaded here:
https://www.gtexportal.org/home/datasets. For other gencode versions see: https://www.gencodegenes.org/human

`-tissue_file`: A tab-delimited file with a list of the eQTL and sQTL results file names per tissue (one row per tissue). Should not contain header. Columns: (1) Tissue name as appears in e/sQTL file names, (2) significant eQTL variant-gene pair file (all significant or fine-mapped variants), (3) Best eQTL per eGene file (contains all genes tested with q-value column), (4) significant sQTL variant-gene pair file (all significant or fine-mapped variants), (5) Best sQTL per sGene file (contains all genes tested with q-value column). The GTEx v8 eQTL and sQTL files were donwloaded from: https://gtexportal.org/home/datasets (e.g., GTEx_TISSUE_NAME.v8.signif_variant_gene_pairs.txt.gz, GTEx_TISSUE_NAME.v8.egenes.txt.gz (from GTEx_Analysis_v8_eQTL.tar), GTEx_TISSUE_NAME.v8.sqtl_signifpairs.txt.gz, GTEx_TISSUE_NAME.v8.sgenes.txt.gz (from GTEx_Analysis_v8_sQTL.tar)). 

Row example (tab-delimited): 
Adipose_Subcutaneous	Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt.gz	Adipose_Subcutaneous.v8.egenes.txt.gz	Adipose_Subcutaneous.v8.sqtl_signifpairs.txt.gz	Adipose_Subcutaneous.v8.sgenes.txt.gz

`-ldkb`: Kilobase window around each GWAS variant in which to search for LD (Default: 5000) (=5Mb)

`-development_skip`: If True, assumes this is a development run where the proxy/trait table have already been generated under output name. Loads the files and skips the LD calculations (Default: FALSE)

`-plink_memory_allocation`: Amount of memory to allocate for PLINK running step only. For annotating ~1000 variants we use: 16384; for annotating 200k variants we recommend: 100152.

### Output files:

`*trait_table.tsv`: Table with list of traits and GWAS variants associated with the traits, and additional information if provided in the GWAS variant input file (`--GWAS_variants' above), such as GWAS variant p-value and effect size and sample ancestry.

`*eQTL_table.tsv`: Contains list of GWAS variants annotated with eVariants that are in LD with them, the eQTL target gene, the distance of the eVariant to the target gene's TSSm and an LD proxy variant for the GWAS variant if was needed to look up the GWAS variant in the eQTL study.

`*eQTL_tissue_table.tsv`: Contains the eQTL p-value and effect size of all eVariants that are in LD with any of the GWAS variants tested, and all variant-gene-tissue combinations for which the eVariant is significant (FDR<5%).

`*sQTL_table.tsv`: Contains list of GWAS variants annotated with sVariants that are in LD with them, the sQTL target gene and intron excision cluster, the distance of the sVariant to the target gene's TSS and an LD proxy variant for the GWAS variant if was needed to look up the GWAS variant in the sQTL study.

`*sQTL_tissue_table.tsv`: Contains the sQTL p-value and effect size of all sVariants that are in LD with any of the GWAS variants tested, and all gene-intron excision cluster-tissue combinations for which the sVariant is significant (FDR<5%).

`*proxy_table.tsv`: Contains LD proxy variants for the inputted GWAS variants and their r2 relative to the GWAS variants taken from a reference panel (e.g., 1000 Genomes Project) if GWAS variant is not found in the QTL study (e.g., GTEx).

## License
The GWASVar2Gene software package is distributed under the terms of the BSD 3-Clause License. See LICENSE.txt file for more details.

## Citation
Rouahana*, Wang* et al., BioRxiv 2021.

Last updated: November 3, 2021
