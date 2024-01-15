import os
import pandas as pd
import numpy as np
import subprocess
import fileinput
import re
import argparse

# Synopsis
#python3 extractgeno.py --infile file.vcf  

##################################### make genotype file from a vcf file #################################################

#input file arguments - a vcf file
parser = argparse.ArgumentParser(description='Extract SNP genotypes from a VCF and encode them as 0,1,2...')
parser.add_argument("--infile")
args = parser.parse_args()

def vcf_extract(vcf_file):
    ### extract chrom, pos, ref, genotype from the vcf file keeping the header
    cmd=f"vcf-to-tab <{vcf_file}"
    ### write this into a tab delimited file
    with open('flygeno.txt', 'w') as outfile:
        outfile.write(os.popen(cmd).read()+"\n")

vcf_extract(args.infile)

### names of the first four columns
cols = ['#CHROM','POS', 'REF']

dt = pd.read_table('flygeno.txt')
dtm = pd.melt(dt, id_vars=cols, var_name='POP', value_name='GEN')
print(dtm)
dtm[['BASE1','BASE2']] = dtm['GEN'].str.split('[:/]', expand=True)
#dtm = dtm.drop(columns = ['COUNTS'])
print(dtm)

#g=re.search("^.*G$", row['REF'] and row['BASE1' and row['BASE2']])

### Add a column for each base 
def genencode(row):

    ### base on the base in REF and BASE1 and BASE2 add 0/0 for reference homozygote, 0/1 encodes for the heterozygote, 1/1 encodes for alternative homozygote; otherwise add zero.
    ### The function considers bi-allelic SNPs, but MNVs are also considered if these haven't been filtered out.
    if re.match("^.*A$", row['REF']) and re.match("^.*A$", row['BASE1']) and re.match("^.*A$", row['BASE2']):
        return ('0/0')
    elif re.match("^.*A$", row['REF']) and re.match("^.*A$", row['BASE1']) and re.match("^.*(?!A)$", row['BASE2']):
        return ('0/1')
    elif re.match("^.*A$", row['REF']) and re.match("^.*(?!A)$", row['BASE1']) and re.match("^.*A$", row['BASE2']):
        return ('1/0')
    elif re.match("^.*A$", row['REF']) and re.match("^.*(?!A)$", row['BASE1']) and re.match("^.*(?!A)$", row['BASE2']):
        return ('1/1')
    elif re.match("^.*T$", row['REF']) and re.match("^.*T$", row['BASE1']) and re.match("^.*T$", row['BASE2']):
        return ('0/0')
    elif re.match("^.*T$", row['REF']) and re.match("^.*T$", row['BASE1']) and re.match("^.*(?!T)$", row['BASE2']):
        return ('0/1')
    elif re.match("^.*T$", row['REF']) and re.match("^.*(?!T)$", row['BASE1']) and re.match("^.*T$", row['BASE2']):
        return ('1/0')
    elif re.match("^.*T$", row['REF']) and re.match("^.*(?!T)$", row['BASE1']) and re.match("^.*(?!T)$", row['BASE2']):
        return ('1/1')
    elif re.match("^.*C$", row['REF']) and re.match("^.*C$", row['BASE1']) and re.match("^.*C$", row['BASE2']):
        return ('0/0')
    elif re.match("^.*C$", row['REF']) and re.match("^.*C$", row['BASE1']) and re.match("^.*(?!C)$", row['BASE2']):
        return ('0/1')
    elif re.match("^.*C$", row['REF']) and re.match("^.*(?!C)$", row['BASE1']) and re.match("^.*C$", row['BASE2']):
        return ('1/0')
    elif re.match("^.*C$", row['REF']) and re.match("^.*(?!C)$", row['BASE1']) and re.match("^.*(?!C)$", row['BASE2']):
        return ('1/1')
    elif re.match("^.*G$", row['REF']) and re.match("^.*G$", row['BASE1']) and re.match("^.*G$", row['BASE2']):
        return ('0/0')
    elif re.match("^.*G$", row['REF']) and re.match("^.*G$", row['BASE1']) and re.match("^.*(?!G)$", row['BASE2']):
        return ('0/1')
    elif re.match("^.*G$", row['REF']) and re.match("^.*(?!G)$", row['BASE1']) and re.match("^.*G$", row['BASE2']):
        return ('1/0')
    elif re.match("^.*G$", row['REF']) and re.match("^.*(?!G)$", row['BASE1']) and re.match("^.*(?!G)$", row['BASE2']):
        return ('1/1')
    else:
        return (0)

def addcode():

    dtm['GENcode'] = dtm.apply (lambda row: genencode(row), axis=1)

    return (dtm)

dtm_gencode = addcode()
print(dtm_gencode)

### save to a text file
dtm_gencode.to_csv(r'gencode.txt', header=True, index=None, sep=' ', mode='w')

### delete the two BASE columns 
dtm_gencode.drop(['BASE1', 'BASE2'], axis = 1, inplace=True)
print(dtm_gencode)

### Add a column for each base 
def SNPencode(row):

    ### encode SNP genotypes as 0,1,2; here if genotype 0/0 give 0; if genotype is either 0/1 or 1/0 give 1 (for 1 SNP); if genotype is 1/1 add 2 (for two SNPs)
	if row['GENcode']=='0/0':
		return ('0')
	elif row['GENcode']=='0/1':
		return ('1')
	elif row['GENcode']=='1/0':
		return ('1')
	elif row['GENcode']=='1/1':
		return ('2')
	else:
		return ('err')

def addSNP():

    dtm_gencode['SNPencode'] = dtm_gencode.apply (lambda row: SNPencode(row), axis=1)

    return (dtm_gencode)

dtm_gencode_SNPcode = addSNP()
print(dtm_gencode_SNPcode)

### save to a text file
dtm_gencode_SNPcode.to_csv(r'SNPcode.csv', header=False, index=None, sep=',', mode='w')
