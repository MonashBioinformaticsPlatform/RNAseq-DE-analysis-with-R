#!/bin/bash
echo "======================"
echo "Downloading hg19 reference (8.4gb)."
echo "This may take hours."
echo "======================"
echo ""
echo ""
echo ""

wget http://biotraining.erc.monash.edu/ref_data/hg19.00.b.array
wget http://biotraining.erc.monash.edu/ref_data/hg19.00.b.tab
wget http://biotraining.erc.monash.edu/ref_data/hg19.fa
wget http://biotraining.erc.monash.edu/ref_data/hg19.files
wget http://biotraining.erc.monash.edu/ref_data/hg19.genes.gtf
wget http://biotraining.erc.monash.edu/ref_data/hg19.log
wget http://biotraining.erc.monash.edu/ref_data/hg19.reads
