#!/home/song7602/miniconda3/envs/py311/bin/python

import sys
import os

vcf = sys.argv[1]

os.system(f'gunzip {vcf}')

unzip_vcf = vcf.replace('.vcf.gz', '.vcf')

open_vcf = open(unzip_vcf, 'r')

fo = open(unzip_vcf.replace('.dose.vcf', '.imputed.vcf'), 'w')
for line in open_vcf:
	if line[0] == '#':
		fo.write(line)
		continue
	col = line.rstrip().split('\t')
	chrom = col[0]
	fo.write('\t'.join([f'chr{chrom}'] + col[1:]) + '\n')
	del line
	del col
fo.close()



