#!/data2/home/song7602/.conda/envs/py311/bin/python3.11

import sys

name = sys.argv[1]
chrom = sys.argv[2]

fo = open(f'{name}.{chrom}.longest_hap.barplot.csv', 'w')
fo.write(',Longest block length,Number of phased variants')

def get_info(info, hap):
	import subprocess
	p1 = subprocess.Popen(['python', '/data2/home/song7602/Tools/GetGermline/GetLongestBlockInfo.py',hap], stdout = subprocess.PIPE)
	[longest_length, phased_numb] = p1.communicate()[0].decode("utf-8").split('\n')[1].split(' / ')
	[longest_length, phased_numb] = (longest_length.split(' : ')[-1], phased_numb.split(' : ')[-1])
	p1.stdout.close()
	fo.write(f'{info},{longest_length},{phased_numb}\n')
	

# Overlap in 3 sequencing platform:
#1. (Tumor Hi-C ∩ Tumor WGS ∩ Blood WGS) ∩ KRG1722
get_info('1. (Tumor Hi-C ∩ Tumor WGS ∩ Blood WGS) ∩ KRG1722', f'{chrom}.{name}.Hi-C_∩_WGS.Tumor_∩_KRG1722_∩_Blood.haplotype')

# Overlap in 2 sequencing platform:
#2.(Tumor Hi-C ∩ Blood WGS) ∩ KRG1722
get_info('2.(Tumor Hi-C ∩ Blood WGS) ∩ KRG1722', f'{chrom}.{name}.Hi-C.Tumor_∩_KRG1722_∩_Blood.haplotype')

#3. (Tumor WGS ∩ Blood WGS) ∩ KRG1722
get_info('3. (Tumor WGS ∩ Blood WGS) ∩ KRG1722', f'{chrom}.{name}.WGS.Tumor_∩_KRG1722_∩_Blood.haplotype')

#4. (Tumor Hi-C ∪ Tumor WGS) ∩ Blood WGS ∩ KRG1722
get_info('4. (Tumor Hi-C ∪ Tumor WGS) ∩ Blood WGS ∩ KRG1722', f'{chrom}.{name}.Hi-C_∪_WGS.Tumor_∩_KRG1722_∩_Blood.haplotype')

# Overlap in 1 sequencing platform:
#5. (Tumor Hi-C ∪ Tumor WGS) ∩ KRG1722
get_info('5. (Tumor Hi-C ∪ Tumor WGS) ∩ KRG1722', f'{chrom}.{name}.Hi-C_∪_WGS.Tumor_∩_KRG1722.haplotype')

#6. (Tumor Hi-C ∩ Tumor WGS) ∩ KRG1722
get_info('6. (Tumor Hi-C ∩ Tumor WGS) ∩ KRG1722', f'{chrom}.{name}.Hi-C_∩_WGS.Tumor_∩_KRG1722.haplotype')

#7. Blood WGS ∩ KRG1722
get_info('7. Blood WGS ∩ KRG1722', f'{chrom}.{name}.WGS.Blood_∩_KRG1722.haplotype')

#8. Tumor WGS ∩ KRG1722
get_info('8. Tumor WGS ∩ KRG1722', f'{chrom}.{name}.WGS.Tumor_∩_KRG1722.haplotype')

#9. Tumor Hi-C ∩ KRG1722
get_info('9. Tumor Hi-C ∩ KRG1722', f'{chrom}.{name}.Hi-C.Tumor_∩_KRG1722.haplotype')



