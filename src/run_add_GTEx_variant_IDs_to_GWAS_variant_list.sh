#!/bin/bash

module purge
module load pandas/2.5.0

./add_GTEx_variant_ID_to_GWAS_variant_list.py
