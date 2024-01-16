#!/bin/bash


#Directory paths
rawdata=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/rawdata
contig=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/contings
workdir=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/work
outdir=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/out

mkdir ${workdir}
mkdir ${outdir}

data=${rawdata}/var.diplo-0.1.vcf
wolcontig=${contig}/wolbachia-contigs.bed
allcontig=${contig}/pilon_round4.bed

# extract list of samples belonging to Brno group and Vienna group that will be compared
bcftools query -l ${data} | grep 'Brno' > ${outdir}/brnolist.txt
bcftools query -l ${data} | grep 'Vienna' > ${outdir}/viennalist.txt

# gzip the vcf file to be able to work with bcftools & index it
bgzip -c ${data} > ${data}.gz
tabix -p vcf ${data}.gz

# filter out all Wolbachia contigs from the bed file containing all contings
bedtools intersect -v -a ${allcontig} -b ${wolcontig} > ${allcontig}.nowol.bed

# keep the filtered contings in the vcf file
bcftools view -R ${allcontig}.nowol.bed ${data}.gz > ${data}.NoWol.vcf
