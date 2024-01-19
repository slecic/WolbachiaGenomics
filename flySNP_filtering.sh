#!/bin/bash


#Directory paths
rawdata=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/rawdata
contig=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/contings
work=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/work
outdir=/mnt/raid10/cerasi/slecic/cerasi/cline/fst/out

mkdir ${outdir}

wolcontig=${contig}/wolbachia-contigs.bed
allcontig=${contig}/pilon_round4.bed

RUN=$1

# extract list of samples belonging to Brno group and Vienna group that will be compared
bcftools query -l ${rawdata}/${RUN} | grep 'Brno' > ${outdir}/brnolist.txt
bcftools query -l ${rawdata}/${RUN} | grep 'Vienna' > ${outdir}/viennalist.txt

# gzip the vcf file to be able to work with bcftools & index it
bgzip -c ${rawdata}/${RUN} > ${rawdata}/${RUN}.gz
tabix -p vcf ${rawdata}/${RUN}.gz

# filter out all Wolbachia contigs from the bed file containing all contings
bedtools intersect -v -a ${allcontig} -b ${wolcontig} > ${allcontig}.nowol.bed

# keep only contings in the filteretd bed file
bcftools view -R ${allcontig}.nowol.bed ${rawdata}/${RUN}.gz > ${work}/${RUN}.NoWol.vcf

# keep only variants that are of type biallelic SNP (if SNPs were called with freebayes this will also include MNVs)
bcftools view --types snps -m 2 -M 2 ${work}/${RUN}.NoWol.vcf > ${work}/${RUN}.NoWol.snp.vcf
