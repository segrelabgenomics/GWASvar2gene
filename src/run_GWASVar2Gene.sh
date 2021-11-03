module purge
#conda activate py39
#module load miniconda/4.7.12
#conda deactivate
#conda activate pandas
module load PLINK/1.9b-6.10

python -u src/GWASVar2Gene.py --GWAS_variants output/out_GWAS_variant_IDs_added.tsv \
 	                   --primary_plink_directory data/PATH_TO_QTL_STUDY_PLINK_FILES \
        	           -primary_plink_prefix GTEx \
                           -primary_ldr 0.8 \
                           -primary_subset data/SAMPLE_ID_QTL_STUDY.txt \
                           -secondary_plink_directory data/PATH_TO_1000_GENOMES_PROJECT_PLINK_FILES \
                           -secondary_plink_prefix 1000_genomes \
                           -secondary_subset data/1KG_EUR_sample_IDs.txt \
                           -GENCODE data/gencode.v26.annotation.gtf \
                           -output_dir output/GWAS_var_GTEx_QTL_Gene_mapping_DATE \
                           -output_prefix Sample_GWAS_variants_r208 \
                           -GENCODE /example/gencode.v26.annotation.gtf \
                           -tissue_file data/QTL_FILENAMES_PER_TISSUE.txt \
			   # For 1000 GWAS variants:
                           -plink_memory_allocation 16384
			   # For GWAS catalog (~100k variants):
                           # -plink_memory_allocation 100152

