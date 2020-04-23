#!/bin/bash
# Find row and print columns in tab-delimited file

#tax_id	GeneID	PubMed_ID
# wildcard=[0-9.]+
# tax_id=3702
# PubMed_ID=31054974
# GI1=831748
# GI2=831889

# Intersection
# patternGI1="^$tax_id\t$GI1\t$wildcard$"
# patternGI2="^$tax_id\t$GI2\t$wildcard$"
# awk -F "\t" -v pat="$patternGI1" '$0 ~ pat {print $1 "\t" $2 "\t" $3}' data/gene2pubmed | \
#   awk -F "\t" -v pat="$patternGI2" '$0 ~ pat {print $1 "\t" $2 "\t" $3}'

# awk -F "\t" -v pat="^$tax_id\t$GI1\t$wildcard$" '$0 ~ pat {print $3}' data/gene2pubmed | \
#   while read pmid; do
#     echo "$pmid"
#     #awk -F "\t" -v pat="^$tax_id\t$GI2\t$pmid$" '$0 ~ pat {print $3; exit}' data/gene2pubmed
#   done


# for GI in $GI1 $GI2; do
#   patternHits="^$tax_id\t$GI\t$PubMed_ID$"
#   patternIDCounts="^$tax_id\t$GI\t$wildcard$"
#   awk -F "\t" -v pat="$patternHits" '$0 ~ pat {print $1 "\t" $2 "\t" $3; exit}' data/gene2pubmed
#   awk -F "\t" -v pat="$patternIDCounts" '$0 ~ pat {print $1 "\t" $2 "\t" $3}' data/gene2pubmed | wc -l
# done

# patternAllGene="^$tax_id\t$wildcard\t$PubMed_ID$"
# awk -F "\t" -v pat="$patternAllGene" '$0 ~ pat {print $1 "\t" $2 "\t" $3;}' data/gene2pubmed


## gene_info
awk -F "\t" -v pat="^9606\t9518\t.+" '$0 ~ pat {print; exit}' data/gene_info