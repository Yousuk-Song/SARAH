#!/home/song7602/miniconda3/envs/py311/bin/python

import sys
import os

vcf = sys.argv[1]
open_vcf = open(vcf, 'r')

def remove_chr(chrom):
	return chrom.split('chr')[-1]

edited_vcf = vcf.replace('.vcf', '.edited.vcf')
fo = open(edited_vcf, 'w')
for line in open_vcf:
	if line[0] == '#':
		fo.write(line)
		continue
	col = line.rstrip().split('\t')
	chrom = col[0]
	fo.write('\t'.join([remove_chr(chrom)] + col[1:]) + '\n')
	del line
	del col
fo.close()

os.system(f'bgzip -c {edited_vcf} > {edited_vcf}.gz')

#os.system(f'tabix -p vcf ${edited_vcf}.gz')



