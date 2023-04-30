#!/data2/home/song7602/.conda/envs/py311/bin/python3.11


import sys

name = sys.argv[1]
chrom = sys.argv[2]

def count_line(vcf):
	n = 0
	open_vcf = open(vcf, 'r')
	for line in open_vcf:
		if line[0] == '#':
			continue
		n += 1
	open_vcf.close()
	return n

Tumor_HiC_and_Blood_WGS = count_line(f'{chrom}.{name}.Hi-C.Tumor_∩_KRG1722_∩_Blood.vcf')

Tumor_HiC_or_WGS_Tumor_and_Blood_WGS = count_line(f'{chrom}.{name}.Hi-C_∪_WGS.Tumor_∩_KRG1722_∩_Blood.vcf')

Blood_WGS = count_line(f'{chrom}.{name}.WGS.Blood_∩_KRG1722.vcf')

Tumor_HiC = count_line(f'{chrom}.{name}.Hi-C.Tumor_∩_KRG1722.vcf')

Tumor_HiC_and_Tumor_WGS = count_line(f'{chrom}.{name}.Hi-C_∩_WGS.Tumor_∩_KRG1722.vcf')

Tumor_WGS_and_Blood_WGS = count_line(f'{chrom}.{name}.WGS.Tumor_∩_KRG1722_∩_Blood.vcf')

Tumor_HiC_and_Tumor_WGS_and_Blood_WGS = count_line(f'{chrom}.{name}.Hi-C_∩_WGS.Tumor_∩_KRG1722_∩_Blood.vcf')

Tumor_HiC_or_Tumor_WGS = count_line(f'{chrom}.{name}.Hi-C_∪_WGS.Tumor_∩_KRG1722.vcf')

Tumor_WGS = count_line(f'{chrom}.{name}.WGS.Tumor_∩_KRG1722.vcf')


n1 = Tumor_HiC_and_Tumor_WGS_and_Blood_WGS
n2 = Tumor_HiC_and_Tumor_WGS - Tumor_HiC_and_Tumor_WGS_and_Blood_WGS
n3 = Tumor_HiC_and_Blood_WGS - Tumor_HiC_and_Tumor_WGS_and_Blood_WGS
n4 = Tumor_WGS_and_Blood_WGS - Tumor_HiC_and_Tumor_WGS_and_Blood_WGS
n5 = Tumor_WGS - (n1 + n2 + n4)
n6 = Blood_WGS - (n1 + n3 + n4)
n7 = Tumor_HiC - (n1 + n2 + n3)

print(f'(1) : {n1}')
print(f'(2) : {n2}')
print(f'(3) : {n3}')
print(f'(4) : {n4}')
print(f'(5) : {n5}')
print(f'(6) : {n6}')
print(f'(7) : {n7}')



