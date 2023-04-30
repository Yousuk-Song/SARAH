#!/data2/home/song7602/.conda/envs/py311/bin/python3.11

import sys

hap = sys.argv[1]

def longest_block(hap):
	open_hap = open(hap, 'r')
	D = {}
	for line in open_hap:
		line = line.rstrip()
		if line[0] == 'B':
			cols = line.split()
			phased_numb = int(cols[6])
			block_len = int(cols[8])
			D[line] = (block_len, phased_numb)
		del line
	longest_len = 0
	block_info = ''
	phased_numb = 0
	for key in D:
		if longest_len < D[key][0]:
			longest_len = D[key][0]
			phased_numb = D[key][1]
			block_info = key
	print(block_info)
	print(f'Longest block length : {longest_len} / Phased variants numb : {phased_numb}')

def Percent(subject, total):
	return round((subject/total)*100, 2)

def hap_to_vcf(hap):
	return hap.split('.haplotype')[0] + '.vcf'

def VcfNumb(vcf):
	open_vcf = open(vcf, 'r')
	total_variant_numb = 0
	for line in open_vcf:
		line = line.rstrip()
		if line[0] != '#':
			total_variant_numb += 1
		del line
	open_vcf.close()
	return total_variant_numb

def N50_block_len(hap):
	open_hap = open(hap, 'r')
	block_len_list = []
	block_numb = 0
	total_phased_numb = 0
	for line in open_hap:
		line = line.rstrip()
		if line[0] == 'B':
			cols = line.split()
			phased_numb = int(cols[6])
			block_numb += 1
			block_len = int(cols[8])
			block_len_list.append((block_len, phased_numb))
			total_phased_numb += phased_numb
			continue
		del line
	open_hap.close()
	block_len_list.sort(reverse=True)
	total_len = sum([i[0] for i in block_len_list])
	sum_len = 0
	phased_numb = 0
	for i in block_len_list:
		len_s = i[0]
		phased_numb = i[1]
		sum_len += len_s
		if sum_len > total_len/2:
			break
	return [len_s, total_phased_numb, phased_numb]

def print_results(n50_len, n50_phased_numb, phased_numb, total_numb):
	print(f'N50 block length(bp) : {n50_len} / Phased variants numb : {n50_phased_numb}')
	print('Number of Phased variants : ', phased_numb)
	print('Number of Unphased variants : ', total_numb-phased_numb)
	print('Number of Total variants : ', total_numb)
	print('Phased/Total variants rate(%) : ', Percent(phased_numb,total_numb))
	print()

res_D = {}

longest_block(hap)
[same_n50_len, same_phased_numb, same_n50_phased_numb] = N50_block_len(hap)
same_total_numb = VcfNumb(hap_to_vcf(hap))
print_results(same_n50_len, same_n50_phased_numb, same_phased_numb, same_total_numb)






