#!/data2/home/song7602/.conda/envs/py311/bin/python3.11

import os
import sys
import threading 

ScriptDir = sys.path[0]
python = '/home/song7602/miniconda3/envs/py311/bin/python'
common_dir = '/home/song7602/KNCC/Data/KRG1722/KRG1722_common_SNV+indels'

blood_wgs = sys.argv[1]
tumor_wgs = sys.argv[2]
tumor_hic = sys.argv[3]

chrom = sys.argv[4]

order = []

# common variant overlap with KRG1722 common variant for blood_wgs, tumor_wgs and tumor_hic + LOH filtering (VAF > 0.1)
import CommonVariantOverlap

Threads = []
for bam in [blood_wgs, tumor_wgs, tumor_hic]:
	t1 = threading.Thread(target=CommonVariantOverlap.overlap, args=(bam, chrom, common_dir))
	Threads.append(t1)

for thread in Threads:
	thread.start()

for thread in Threads:
	thread.join()

# output vcf's name
tumor_hic_vcf = tumor_hic.split('/')[-1].replace('.bam', '.KRG1722.vcf')
tumor_wgs_vcf = tumor_wgs.split('/')[-1].replace('.bam', '.KRG1722.vcf')
blood_wgs_vcf = blood_wgs.split('/')[-1].replace('.bam', '.KRG1722.vcf')

# combine vcfs from hi-c and wgs

def combine_vcf(tumor_hic_vcf, tumor_wgs_vcf, name1, name2):
	open_tumor_hic_vcf = open(tumor_hic_vcf, 'r')
	open_tumor_wgs_vcf = open(tumor_wgs_vcf, 'r')

	combined_vcf = tumor_hic_vcf.replace(name1, name2)
	fo = open(combined_vcf, 'w')

	tumor_hic_vcf_D = {}
	for line in open_tumor_hic_vcf:
		if line[0] == '#':
			fo.write(line)
			continue
		col = line.rstrip().split('\t')
		chrom = col[0]
		pos = int(col[1])
		ref = col[3]
		alt = col[4]
		info = col[7]
		tumor_hic_vcf_D[(chrom, pos)] = [ref, alt, info]
		del line
	open_tumor_hic_vcf.close()

	for line in open_tumor_wgs_vcf:
		if line[0] != '#':
			col = line.rstrip().split()
			chrom = col[0]
			pos = int(col[1])
			ref = col[3]
			alt = col[4]
			info = col[7]
			if tumor_hic_vcf_D.get((chrom,pos), -1) == -1: #  exists only in wgs
				continue
			[ref, alt, info] = tumor_hic_vcf_D[(chrom, pos)]
			fo.write('\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', info, 'GT', '0/1']) + '\n')
		del line
	open_tumor_wgs_vcf.close()
	fo.close()

# Combine: Tumor Hi-C + Tumor WGS
combine_vcf(tumor_hic_vcf, tumor_wgs_vcf, '.Tumor.HiC', '.Tumor.HiC-WGS')
combined_vcf = tumor_hic_vcf.replace('.Tumor.HiC', '.Tumor.HiC-WGS')

# Combine: Combined Tumor + Blood WGS
combine_vcf(combined_vcf, blood_wgs_vcf , '.Tumor.HiC-WGS', '.Tumor.HiC-WGS.Blood.WGS')
final_germline_vcf = combined_vcf.replace('.Tumor.HiC-WGS', '.Tumor.HiC-WGS.Blood.WGS')

# HAPCUT
def hapcut(vcf, bam):
	hapcutfragment = vcf.replace('.vcf', '.fragment')
	hapcutfile = vcf.replace('.vcf', '.haplotype')

	order1 = ['extractHAIRS', '--hic 1', '--bam', bam, '--VCF', vcf, '--out', hapcutfragment, '--indels 1']
	order2 = ['HAPCUT2', '--hic 1', '--fragments', hapcutfragment, '--VCF', vcf, '--output', hapcutfile]

	os.system(' '.join(order1))
	os.system(' '.join(order2))

#hapcut(final_germline_vcf, tumor_hic)







