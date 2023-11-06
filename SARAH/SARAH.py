#!/data2/home/song7602/miniconda3/bin/python3.11

"""
The purpose of this python3 script is to separate WGS and Hi-C bam files into alleles based on haplotypes. (for archon http://devmh.iptime.org:25881/main/main.do)
Author: Yousuk Song
Last updated date: 2023.09.13
"""

import argparse
import datetime
import os
import sys
import pysam
import time
import multiprocessing as mp

# Dependency PATH
gatk = '/data2/home/song7602/1.Scripts/gatk-4.2.6.1/gatk'
bamtools = '/home/eaststar0/miniconda3/bin/bamtools'
samtools = '/home/eaststar0/miniconda3/bin/samtools'
cooler = '/data2/home/song7602/miniconda3/bin/cooler'
ToolDir = sys.path[0]

def parse_arguments():
	parser = argparse.ArgumentParser(description='Type input for SARAH: (1)haplotype, (2)bam file (3)chromosome name, (4)threads') #change
	required = parser.add_argument_group('required arguments')
	required.add_argument('-n', '--NAME',
                              type=str,
                              required=True,
                              help='Sample_name and Tumor/Normal Informations (ex) HNT_112.Tumor')

	required.add_argument('-hap', '--HAPLOTYPE',
			      required=True,
			      type=str,
                              help='Chromosome-speparated haplotype')

	required.add_argument('-b', '--BAM',
                              type=str,
			      required=True,
                              help='Input bam file')

	required.add_argument('-c', '--CHROM',
                              type=str,
			      required=True,
                              help='Target chromomsome name for phasing')

	required.add_argument('-s', '--SCORE',
			      default=100, 
                              type=int,
                              help='Mimimum hapcut quality score for phasing')

	required.add_argument('-r', '--REFERENCE',
  			        default='hg38',
				type=str,
				help='Reference Genome (hg19, hg38)')
	
	required.add_argument('-d', '--DATA',
  			      default='HiC',
                              help='Data type of input bam file (WGS, HiC)')

	optional = parser.add_argument_group('optional arguments')
	optional.add_argument('-t', '--THREADS',
                              default=10,
                              help='Threads for multi-proecessing')
	optional.add_argument('-N', '--NUMBER',
			      type=int,
                              default= 1,
                              help='Number of blocks for phasing (1, 2)')
	args = parser.parse_args()
	return args


def get_current_datetime():
	return str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))


def print_log(msg):
	print("[" + get_current_datetime() + "] " + str(msg))


def check_if_finished(output, message):
	import time
	finish = 1
	while finish:
		if not os.path.isfile(output):
			if message:
				print(message)
			time.sleep(1)
		else:
			finish = 0

def longest_block(chrom, hap, numb):
	global major_hap 
	### 2-1. make haplotype directory ###
	print_log('Thank You For Using SARAH (Seperation of Aligned Reads into Alleles based-on Haplotypes)')
	print_log(f"Getting longest block from {hap} starting...")
	open_hap = open(hap, 'r')

	block_dict = {}
	block_name = ''
	n = 0
	for line in open_hap:
		line = line.rstrip()
		if line[0] == 'B':
			if n == 0:
				n = 1
				block_name = line
		elif line[0] != '*':
			col = line.split()
			if n == 1:
				block_name += ' ' + 'chromosome' + ' ' + col[3]
				n = 2
				block_dict[block_name] = [line]
			else:
				block_dict[block_name] += [line]
		elif line[0] == '*':
			if n == 2:
				block_dict[block_name] += [line]
				block_name = ''
				n = 0
	for block_name in list(block_dict.keys()):
		if block_name.split()[-1] != chrom:
			del block_dict[block_name]
	top2_block_list = sorted(list(block_dict.keys()), key = lambda x:int(x.split()[8]), reverse = True)[:2]
	top1_block = top2_block_list[0]
	top2_block = top2_block_list[1]

	### 2-2. get longest block ###
	# longest haplotype blocks for phasing
	
	top1_hap = f'{hap}_{chrom}_top1_haplotype'
	top2_hap = f'{hap}_{chrom}_top2_haplotype'

	out_top1_hap = open(top1_hap, 'w')
	out_top1_hap.write(top1_block + '\n')
	for line in block_dict[top1_block]:
		out_top1_hap.write(line + '\n')
	out_top1_hap.close()
	if numb == 2:
		out_top2_hap = open(top2_hap, 'w')
		out_top2_hap.write(top2_block + '\n')
		for line in block_dict[top2_block]:
			out_top2_hap.write(line + '\n')
		out_top2_hap.close()
	del block_dict
	print_log(f"Getting longest block from {hap} finished")


def align_seq(read):
	init = 0
	next_init = 0
	seq = ''
	align = ''
	valid_cigar = []

	pos = read.pos - 1 # to build ins, del dict
	ins_seq = ''
	del_seq = ''
	del_seq2 = ''
	for cigar_info in read.cigar:
		if cigar_info[0] in [0,1,2]:
			valid_cigar.append(cigar_info)
	for cigar_info in valid_cigar:
		cigar = cigar_info[0]
		cigar_len = int(cigar_info[1])
		if cigar == 0: # match or mismatch
			align += 'M'*cigar_len
			for base_pos in range(init, init + cigar_len):
				base = read.query_alignment_sequence[base_pos]
				seq += base
				next_init += 1
			init = next_init
		elif cigar == 1: # insertion
			align += 'I'*cigar_len
			for base_pos in range(init, init + cigar_len):
				base = read.query_alignment_sequence[base_pos]
				seq += base
				next_init += 1
			init = next_init
		elif cigar == 2: # deletion
			seq += '-'*cigar_len
			align += 'D'*cigar_len
	return [seq, align]


def avoid_soft_clipping(read, single_pos):
	is_valid = 0
	if read.cigarstring != None:
		start_cigar = read.cigar[0]
		end_cigar = read.cigar[-1]
		if start_cigar[0] == 4: # front
			if read.pos <= single_pos: # front valid
				if end_cigar[0] == 4: # end
					if single_pos <= read.pos + sum([i[1] for i in read.cigar[:-1] if i[0] != 1]): # back valid
						is_valid = 1
				else: # only front
					is_valid = 1

		elif end_cigar[0] == 4: # not front, only back
			if single_pos <= read.pos + sum([i[1] for i in read.cigar[:-1]]): # back valid
				is_valid = 1
		else: # no soft clip
			is_valid = 1
	return is_valid

def SplitBam(bam, hap, process_dir, chrom, pos, ref, alt, order, score, name):
	hap_name = hap.split('/')[-1]
	a_bam = f'{process_dir}/A.{hap_name}.{name}.{chrom}.{pos}.bam'
	b_bam = f'{process_dir}/B.{hap_name}.{name}.{chrom}.{pos}.bam'
	c_bam = f'{process_dir}/C.{hap_name}.{name}.{chrom}.{pos}.bam'
	_pos = pos - 1
	del_dict = {}
	ins_dict = {}
	pos_dict = {}
	if len(ref) > len(alt): # deletion -> '-' in alt
		alt = alt + '-'*(len(ref)-1)
	elif len(ref) < len(alt): # insertion -> 'I' in ref
		ref = ref + 'I'*(len(alt)-1)
	if order == ['1', '0']:
		[A,B] = [alt, ref]
	elif order == ['0', '1']:
		[A,B] = [ref, alt]
	if '-' in alt: # del_dict
		del_dict[_pos] = {'A':A, 'B':B}
	if 'I' in ref: # ins_dict
		ins_dict[_pos] = {'A':A, 'B':B}
	if [A,B] != ['', '']:
		pos_dict[_pos] = {'A':A, 'B':B}
	samfile = pysam.AlignmentFile(bam, 'rb')
	A_bam = pysam.AlignmentFile(a_bam, 'wb', template = samfile)
	B_bam = pysam.AlignmentFile(b_bam, 'wb', template = samfile)
	C_bam = pysam.AlignmentFile(c_bam, 'wb', template = samfile)

	for single_pos in list(pos_dict.keys()):
		A_B_base = [pos_dict[single_pos][x] for x in ['A', 'B']]
		a_base = A_B_base[0]
		b_base = A_B_base[1]
		for read in samfile.fetch(chrom, single_pos, single_pos+1):
			if read.flag & 1 and read.flag & 4 == read.flag & 8 == read.flag & 256 == read.flag & 512 == read.flag & 1024 == 0:
#			if read.flag & 8 == read.flag & 256 == read.flag & 512 == read.flag & 1024 == 0:
				if avoid_soft_clipping(read, single_pos):
					query_pos = single_pos - read.pos
					align_res = align_seq(read)
					query_seq = align_res[0]
					cigar_seq = align_res[1]
					query_pos += cigar_seq[:query_pos].count('I')
					if query_pos >= 0:
						base_allele = "c" # initial value
						indel = ''
						if single_pos in ins_dict: # ins pos
							base = query_seq[query_pos:query_pos + len(ins_dict[single_pos]['A'])]
							base_cigar = cigar_seq[query_pos:query_pos + len(ins_dict[single_pos]['A'])]
							if base not in [ins_dict[single_pos]['A'], ins_dict[single_pos]['B']]:
								if base_cigar != 'M'*len(base_cigar): # it's ins. but not alternative ins
									pass
								else: # no ins  
									if 'I' in pos_dict[single_pos]['A']:	
										base_allele = 'a'
									elif 'I' in pos_dict[single_pos]['B']:
										base_allele = 'b'
							else: # found matched query indel, but confusing
								if 'I' in base_cigar: # real ins
									if base == ins_dict[single_pos]['A']:
										base_allele = 'a'
									elif base == pos_dict[single_pos]['B']:
										base_allele = 'b'

								else: # no ins 
									if 'I' in pos_dict[single_pos]['A']:	
										base_allele = 'a'
									elif 'I' in pos_dict[single_pos]['B']:
										base_allele = 'b'

						elif single_pos in del_dict: # del
							indel = 'del'
							base = query_seq[query_pos:query_pos + len(del_dict[single_pos]['A'])]
						
							if base == del_dict[single_pos]['A']:
								base_allele = "a"
							elif base == del_dict[single_pos]['B']:
								base_allele = "b"

						else: # no indel
							base = query_seq[query_pos]
							if base == pos_dict[single_pos]['A']:
								base_allele = "a"
							elif base == pos_dict[single_pos]['B']:
								base_allele = "b"
						if base_allele != "c": # single base match with haps
							real_pos_base_dic = {} 
							all_allele_match = 1
							L = []
							for all_pos in list(pos_dict.keys()): # every hap variants in single read
								all_query_pos = all_pos - read.pos
								if 0 <= all_query_pos < len(query_seq):
									indel = ''
									all_base_allele = 'c'
									if all_pos in ins_dict: # all_ins
										all_base = query_seq[all_query_pos:all_query_pos + len(ins_dict[all_pos]['A'])]
										all_base_cigar = cigar_seq[all_query_pos:all_query_pos + len(ins_dict[all_pos]['A'])]
										if all_base not in [ins_dict[all_pos]['A'], ins_dict[all_pos]['B']]:
											if all_base_cigar != 'M'*len(all_base_cigar): # it's ins. but not alternative ins
												pass
											else: # no ins
												if 'I' in pos_dict[all_pos]['A']:
													all_base_allele = 'a'
												elif 'I' in pos_dict[all_pos]['B']:
													all_base_allele = 'b'			
										else: # found matched query indel, but confusing
											if 'I' in all_base_cigar: # real ins
												if all_base == ins_dict[all_pos]['A']:
													all_base_allele = 'a'
												elif all_base == pos_dict[all_pos]['B']:
													all_base_allele = 'b'
											else: # no ins 
												if 'I' in pos_dict[all_pos]['A']:
													all_base_allele = 'a'
												elif 'I' in pos_dict[all_pos]['B']:
													all_base_allele = 'b'
													
									elif all_pos in del_dict: # all del
										indel = 'del'
										all_base = query_seq[all_query_pos:all_query_pos + len(del_dict[all_pos]['A'])]
										if all_base == del_dict[all_pos]['A']:
											all_base_allele = "a"
										elif all_base == del_dict[all_pos]['B']:
											all_base_allele = "b"
									else: # no indel
										all_base = query_seq[all_query_pos]
										if all_base == pos_dict[all_pos]['A']:
											all_base_allele = "a"
										elif all_base == pos_dict[all_pos]['B']:
											all_base_allele = "b"

									if all_base_allele != base_allele:
										all_allele_match = 0

							if all_allele_match:
								if base_allele == "a": # allele match with hap_A
									A_bam.write(read)
								elif base_allele == "b": # allele match with hap_B
									B_bam.write(read)
								else:
									C_bam.write(read)
							else:
								C_bam.write(read)
						else:
							C_bam.write(read)
					del align_res
			del read
		del A_B_base
		del single_pos

	del pos_dict #del read_list
	del del_dict
	del ins_dict

	samfile.close()
	A_bam.close()
	B_bam.close()
	C_bam.close()

	os.system(f'{samtools} sort -o {a_bam.replace(".bam", ".sorted.bam")} {a_bam}')
	os.system(f'{samtools} index {a_bam.replace(".bam", ".sorted.bam")}')
	os.system(f'rm {a_bam}')

	os.system(f'{samtools} sort -o {b_bam.replace(".bam", ".sorted.bam")} {b_bam}')
	os.system(f'{samtools} index {b_bam.replace(".bam", ".sorted.bam")}')
	os.system(f'rm {b_bam}')

	os.system(f'{samtools} sort -o {c_bam.replace(".bam", ".sorted.bam")} {c_bam}')
	os.system(f'{samtools} index {c_bam.replace(".bam", ".sorted.bam")}')
	os.system(f'rm {c_bam}')
	os.system("ps aux | grep defunct | awk '{print $2}' | xargs kill -9 2> /dev/null")
	print_log(f'{chrom} {pos}')

def split_phased_bam(phased_bam_for_merge, phased_split_bam, phased_bam_list, split_n):
	n1 = 0
	n2 = 0
	for line in phased_bam_list:
		n1 += 1
		phased_split_bam.append(line)
		if n1 == split_n:
			n2 += 1
			phased_bam_for_merge[n2] = phased_split_bam[:]
			phased_split_bam.clear()
			n1 = 0
	if phased_split_bam != []:
		n2 += 1
		phased_bam_for_merge[n2] = phased_split_bam[:]
		phased_split_bam.clear()

# merging bam files 
def merge_bam(allele, hap_name, input_phased_split_bam, numb, process_dir,bamtools, name): # for limiting the number of files per bamtools process
	print(f'merging {len(input_phased_split_bam)} bamfiles into {numb}st {allele} merged bam...')

	input_for_bamtools = f'{process_dir}/{allele}.{hap_name}.{numb}.input.txt'

	open_input_for_bamtools = open(input_for_bamtools, 'w')
	for line in input_phased_split_bam:
		input_bam = f'{process_dir}/{allele}.{line}'
		open_input_for_bamtools.write(input_bam + '\n')
	open_input_for_bamtools.close()

	bam_0 = f'{process_dir}/{allele}.{hap_name}.{name}.bam_0.bam.numb:{numb}.bam'
	bamtools_order = [bamtools, 'merge', '-list', input_for_bamtools, '-out', bam_0]
	os.system(' '.join(bamtools_order))
	check_if_finished(bam_0, 'merging bamfiles not finished yet...')

	os.system(f'{samtools} index {bam_0}')
	os.system(f'rm {input_for_bamtools}')

def remove_dup(input_bam, rmdup_bam):
	order = [samtools, 'rmdup', '-S', input_bam, rmdup_bam]
	os.system(' '.join(order))
	os.system(f'rm {input_bam}')
	os.system(f'rm {input_bam}.bai')
	os.system(f'{samtools} index {rmdup_bam}')

# seperate chromosome by hap snp
def Sarah(chr_hap, bam, chrom, name, score, data_type, bamtools, multi_thread):
	n = 1
	print_log('Building SNP-Overlap Bam Files')
	multi_thread = int(multi_thread)
	open_hap = open(chr_hap, 'r')
	hap_name = chr_hap.split('/')[-1]
	if data_type == 'HiC':
		process_dir = f'{chrom}.{name}.phasing.hic.dir'
	else:
		process_dir = f'{chrom}.{name}.phaisng.wgs.dir'
	os.system(f'mkdir {process_dir}')

	phased_bam_list = []
	hap_line = []
	allele_numb = 0
	processes = []
	for line in open_hap:
		line = line.rstrip()
		if line[0] not in ["*", "B"]:
			if '-' not in [line.split()[1], line.split()[2]]:
				hap_score =  float(line.split()[-2]) 
				if hap_score >= score:
					hap_line.append(line.split())
					allele_numb += 1
					col = line.split()
					chrom = col[3]
					pos = int(col[4])
					ref = col[5]
					alt = col[6]
					order = col[1:3]

					phased_bam = f'{hap_name}.{name}.{chrom}.{pos}.sorted.bam'
					phased_bam_list.append(phased_bam)
					SplitBam(bam, chr_hap, process_dir, chrom, pos, ref, alt, order, score, name)
					print_log(f"{n} variants phasing done!")
					n += 1
#					process = mp.Process(target=SplitBam, args=(bam, hap, process_dir, chrom, pos, ref, alt, order, score, name))
#					processes.append(process)
#	split_process = [processes[i:i+multi_thread] for i in range(0, len(processes), multi_thread)]
#	split_numb = 0
#	for process_list in split_process:
#		for process in process_list:
#			process.start()
#		for process in mp.active_children():
#			process.join()
#		split_numb += multi_thread
#		print_log(f"{split_numb} variants phasing done!")
	

	print_log('done')
	#### Merging Bam Files ####
	print_log('Separation of Allele-Specific Bam Files Finished')
	print_log('Merging Bam Files ...')


	split_n = 1000
	phased_bam_for_merge = {}
	phased_split_bam = []
	split_phased_bam(phased_bam_for_merge, phased_split_bam, phased_bam_list, split_n)

	for numb in phased_bam_for_merge:
		merge_bam('A', hap_name, phased_bam_for_merge[numb], numb, process_dir, bamtools, name)
		check_if_finished(f'{process_dir}/A.{hap_name}.{name}.bam_0.bam.numb:{numb}.bam', 'slow down...')
		merge_bam('B', hap_name, phased_bam_for_merge[numb], numb, process_dir, bamtools, name)
		check_if_finished(f'{process_dir}/B.{hap_name}.{name}.bam_0.bam.numb:{numb}.bam', 'slow down...')
		merge_bam('C', hap_name, phased_bam_for_merge[numb], numb, process_dir, bamtools, name)
		check_if_finished(f'{process_dir}/C.{hap_name}.{name}.bam_0.bam.numb:{numb}.bam', 'slow down...')

	A_merged_bam0 = f'{process_dir}/A.{hap_name}.{name}.bam_0.bam'
	B_merged_bam0 = f'{process_dir}/B.{hap_name}.{name}.bam_0.bam'
	C_merged_bam0 = f'{process_dir}/C.{hap_name}.{name}.bam_0.bam'

	if data_type == 'HiC':
		rmdup_bam = f'{hap_name}.{name}.hic.bam'
	else:
		rmdup_bam = f'{hap_name}.{name}.wgs.bam'

	A_rmdup_bam = f'{process_dir}/A.{rmdup_bam}'
	B_rmdup_bam = f'{process_dir}/B.{rmdup_bam}'
	C_rmdup_bam = f'{process_dir}/C.{rmdup_bam}'

	A_merge_list = f'{process_dir}/A.{hap_name}.merged.list.txt'
	B_merge_list = f'{process_dir}/B.{hap_name}.merged.list.txt'
	C_merge_list = f'{process_dir}/C.{hap_name}.merged.list.txt'

	A_merge_list_file = open(A_merge_list, 'w')
	B_merge_list_file = open(B_merge_list, 'w')
	C_merge_list_file = open(C_merge_list, 'w')


	if len(phased_bam_for_merge) == 1: # total variant number under 1000
		os.system(f'mv {A_merged_bam0}.numb:{numb}.bam {A_merged_bam0}\n')
		os.system(f'mv {B_merged_bam0}.numb:{numb}.bam {B_merged_bam0}\n')
		os.system(f'mv {C_merged_bam0}.numb:{numb}.bam {C_merged_bam0}\n')

		os.system(f'mv {A_merged_bam0}.numb:{numb}.bam.bai {A_merged_bam0}.bai\n')
		os.system(f'mv {B_merged_bam0}.numb:{numb}.bam.bai {B_merged_bam0}.bai\n')
		os.system(f'mv {C_merged_bam0}.numb:{numb}.bam.bai {C_merged_bam0}.bai\n')

		os.system(f'{samtools} index {A_merged_bam0}')
		os.system(f'{samtools} index {B_merged_bam0}')
		os.system(f'{samtools} index {C_merged_bam0}')

		remove_dup(A_merged_bam0, A_rmdup_bam)
		remove_dup(B_merged_bam0, B_rmdup_bam)
		remove_dup(C_merged_bam0, C_rmdup_bam)
		
	else: # total variant number over 1000
		rm_list = []
		for numb in phased_bam_for_merge:
			A_merge_list_file.write(f'{A_merged_bam0}.numb:{numb}.bam\n')
			B_merge_list_file.write(f'{B_merged_bam0}.numb:{numb}.bam\n')
			C_merge_list_file.write(f'{C_merged_bam0}.numb:{numb}.bam\n')
			rm_list += [f'{A_merged_bam0}.numb:{numb}.bam', f'{A_merged_bam0}.numb:{numb}.bam.bai']
			rm_list += [f'{B_merged_bam0}.numb:{numb}.bam', f'{B_merged_bam0}.numb:{numb}.bam.bai']
			rm_list += [f'{C_merged_bam0}.numb:{numb}.bam', f'{C_merged_bam0}.numb:{numb}.bam.bai']
		A_merge_list_file.close()
		B_merge_list_file.close()
		C_merge_list_file.close()

		bamtools_order_A = [bamtools, 'merge', '-list', A_merge_list, '-out', A_merged_bam0]
		bamtools_order_B = [bamtools, 'merge', '-list', B_merge_list, '-out', B_merged_bam0]
		bamtools_order_C = [bamtools, 'merge', '-list', C_merge_list, '-out', C_merged_bam0]

		os.system(' '.join(bamtools_order_A))
		check_if_finished(A_merged_bam0, 'merging not finished yet')
		os.system(f'{samtools} index {A_merged_bam0}')
		check_if_finished(A_merged_bam0 + '.bai', 'merging not finished yet')

		os.system(' '.join(bamtools_order_B))
		check_if_finished(B_merged_bam0, 'merging not finished yet')
		os.system(f'{samtools} index {B_merged_bam0}')
		check_if_finished(B_merged_bam0 + '.bai', 'merging not finished yet')

		os.system(' '.join(bamtools_order_C))
		check_if_finished(C_merged_bam0, 'merging not finished yet')
		os.system(f'{samtools} index {C_merged_bam0}')
		check_if_finished(C_merged_bam0 + '.bai', 'merging not finished yet')

		remove_dup(A_merged_bam0, A_rmdup_bam)
		check_if_finished(A_rmdup_bam + '.bai', 'merging not finished yet')

		remove_dup(B_merged_bam0, B_rmdup_bam)
		check_if_finished(B_rmdup_bam + '.bai', 'merging not finished yet')

		remove_dup(C_merged_bam0, C_rmdup_bam)
		check_if_finished(C_rmdup_bam + '.bai', 'merging not finished yet')

		os.system(f"rm {' '.join(rm_list)}")

	for bam_f in phased_bam_list:
		os.system(f'rm {process_dir}/A.{bam_f} {process_dir}/A.{bam_f}.bai')
		os.system(f'rm {process_dir}/B.{bam_f} {process_dir}/B.{bam_f}.bai')
		os.system(f'rm {process_dir}/C.{bam_f} {process_dir}/C.{bam_f}.bai')

	os.system(f'rm {process_dir}/A.{hap_name}.merged.list.txt {process_dir}/B.{hap_name}.merged.list.txt {process_dir}/C.{hap_name}.merged.list.txt ')

	print_log('\ndone\n')
	print_log('Merging Bam Files Finished\n')
	if data_type != 'HiC':
		os.system(f'mv {process_dir}/A*wgs.bam ./')
		os.system(f'mv {process_dir}/A*wgs.bam.bai ./')
		os.system(f'mv {process_dir}/B*wgs.bam ./')
		os.system(f'mv {process_dir}/B*wgs.bam.bai ./')
		os.system(f'mv {process_dir}/C*wgs.bam ./')
		os.system(f'mv {process_dir}/C*wgs.bam.bai ./')
		os.system(f'rm -rf {process_dir}')

	print_log('Separation of Allele-Specific Bam Files Processes are Done\n')

def GetMate(bam, phased_bam, read_chrom, mate_chrom, process_dir):
	get_mate_phased_bam = f'{process_dir}/{phased_bam.replace(".bam", f".mate.{mate_chrom}.bam")}'
	sorted_get_mate_bam = get_mate_phased_bam.replace(".bam", ".sort.bam")

	samfile = pysam.AlignmentFile(f'{process_dir}/{phased_bam}', "rb") # bam
	template_file = pysam.AlignmentFile(bam, "rb") 
	outfile = pysam.AlignmentFile(get_mate_phased_bam, "wb", template = samfile)

	bam_name = bam.split("/")[-1].split(f"{read_chrom}.")[-1]

	for read in samfile:
		if -1 not in [read.tid, read.rnext] and read.is_read1:
			if read.flag & 1 and read.flag & 4 == read.flag & 8 == read.flag & 256 == read.flag & 1024 == 0:
				next_ref = read.next_reference_name
				if next_ref == mate_chrom:
					try:
						mate_read = template_file.mate(read)
						if mate_read.flag & 4 == mate_read.flag & 8 == mate_read.flag & 256 == mate_read.flag & 512 == mate_read.flag & 1024 == 0:
							print_log(f'{read_chrom}\t{read.pos}\t{mate_chrom}\t{mate_read.pos}')
							outfile.write(mate_read)
					except ValueError:
						print_log(f'{read_chrom}\t{read.pos}\tNo mate reads found')
						continue
		del read
	template_file.close()
	samfile.close() 
	outfile.close() # warning: please make sure to have the pysam alignmentfile close before using gatk merge!!
	os.system(f'{samtools} sort -o {sorted_get_mate_bam} {get_mate_phased_bam}')
	os.system(f'{samtools} index {sorted_get_mate_bam}')
	os.system(f'rm {get_mate_phased_bam}')


def mate_merge(gatk, phased_bam, out_list, process_dir):

	final_out0 = phased_bam.replace('.bam', '.get.mate.merged_0.bam')
	final_out = final_out0.replace('.get.mate.merged_0.bam', '.get.mate.merged.bam')

	order1 = [gatk, 'MergeSamFiles']
	for input_bam in out_list:
		order1.append(f'I={input_bam}')
	order1.append(f'I={process_dir}/{phased_bam}')

	order1.append(f'O={final_out0}')
	os.system(' '.join(order1)) # picard merge bam files
	check_if_finished(final_out0, f'merging {final_out0} not finished yet')

	order2 = [samtools, 'rmdup', '-S', final_out0, final_out]
	os.system(' '.join(order2)) # samtools rmdup
	check_if_finished(final_out, f'rmdup of {final_out} not finished yet')

	order3 = ['rm'] + out_list + [i+'.bai' for i in out_list] + [final_out0]
	os.system(' '.join(order3)) # remove input files

mate_list = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9',
	'chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17',	
	'chr18','chr19','chr20','chr21','chr22', 'chrX', 'chrY', 'HPV16', 'HPV33', 'HPV35']

def get_chrom_mate(gatk, bam, chr_hap, read_chrom, name):
	process_dir = f'{read_chrom}.{name}.phasing.hic.dir'
	print_log('SARAH Get_Mate_Reads_By_Chromosomes Initializes')
	print_log('Getting Mate Reads Starting...')

	hap_bam = f'{chr_hap}.{name}.hic.bam'
	A_bam = 'A.' + hap_bam
	B_bam = 'B.' + hap_bam
	C_bam = 'C.' + hap_bam
	out_list_a = []
	out_list_b = []
	out_list_c = []

	processes = []
	for mate_chrom in mate_list:
#		GetMate(bam, A_bam, read_chrom, mate_chrom, process_dir)
#		GetMate(bam, B_bam, read_chrom, mate_chrom, process_dir)
#		GetMate(bam, C_bam, read_chrom, mate_chrom, process_dir)
		processes.append(mp.Process(target=GetMate, args=(bam, A_bam, read_chrom, mate_chrom, process_dir)))
		processes.append(mp.Process(target=GetMate, args=(bam, B_bam, read_chrom, mate_chrom, process_dir)))
		processes.append(mp.Process(target=GetMate, args=(bam, C_bam, read_chrom, mate_chrom, process_dir)))
			
		output_a = f'{process_dir}/{A_bam.replace(".bam", ".mate")}.{mate_chrom}.sort.bam'
		output_b = f'{process_dir}/{B_bam.replace(".bam", ".mate")}.{mate_chrom}.sort.bam'
		output_c = f'{process_dir}/{C_bam.replace(".bam", ".mate")}.{mate_chrom}.sort.bam'

		out_list_a.append(output_a)
		out_list_b.append(output_b)
		out_list_c.append(output_c)
	for process in processes:
		process.start()
	for process in mp.active_children():
		process.join()

	print_log('Getting Mate Reads Finished')
	print_log('Merging Bam Files Starting...')

	a_final = A_bam.replace(".bam", ".get.mate.merged.bam")
	b_final = B_bam.replace(".bam", ".get.mate.merged.bam")
	c_final = C_bam.replace(".bam", ".get.mate.merged.bam")

	mate_merge(gatk, A_bam, out_list_a, process_dir)
	check_if_finished(a_final, f'{a_final} not yet merged')
	mate_merge(gatk, B_bam, out_list_b, process_dir)
	check_if_finished(b_final, f'{b_final} not yet merged')
	mate_merge(gatk, C_bam, out_list_c, process_dir)
	check_if_finished(c_final, f'{c_final} not yet merged')

	# sorting
	a_sort_final = a_final.replace(".bam", ".sort.bam")
	b_sort_final = b_final.replace(".bam", ".sort.bam")
	c_sort_final = c_final.replace(".bam", ".sort.bam")

	os.system(f'{samtools} sort -o {a_sort_final} {a_final}')
	os.system(f'{samtools} sort -o {b_sort_final} {b_final}')
	os.system(f'{samtools} sort -o {c_sort_final} {c_final}')

	# indexing 
	os.system(f'{samtools} index {a_sort_final}')
	os.system(f'{samtools} index {b_sort_final}')
	os.system(f'{samtools} index {c_sort_final}')

	# remove prev
	os.system(f'rm {a_final}'); os.system(f'rm {b_final}'); os.system(f'rm {c_final}');
	os.system(f'rm -rf {process_dir}')
	
	print_log(f"Merging Block {chr_hap}'s Bam Files Finished")
	print_log(f"Getting Mate Reads Of Block {chr_hap} is Done")

def bam2cool(bam, ref):
	ToolDir = sys.path[0]
	out_prefix=bam.split('/')[-1].split('.bam')[0]
	bin_size=1000
	res = '1kb'
	ncores=8
	max_split=2
	chrom_size = f'{ToolDir}/bin/hg19.chr.size.txt'
	os.system(f'{sys.path[0]}/bam2pairs.pl -c {chrom_size} {bam} {out_prefix}')

	pairs_file= f'{out_prefix}.bsorted.pairs.gz'
	os.system(f'cooler cload pairix -p {ncores} -s {max_split} {chrom_size}:{bin_size} {pairs_file} {out_prefix}.cool')
	os.system(f'cooler zoomify -o {out_prefix}.mcool -r "1000,5000,10000,25000,50000,100000,200000,500000,1000000,2000000,5000000" {out_prefix}.cool') # Raw mcool file

def final_bam(allele, chr_hap, name):
	return f'{allele}.{chr_hap}.{name}.hic.get.mate.merged.sort.bam'


if __name__ == '__main__':
	args = parse_arguments()
	longest_block(args.CHROM, args.HAPLOTYPE, args.NUMBER)
	top1_hap = f"{args.HAPLOTYPE.split('/')[-1]}_{args.CHROM}_top1_haplotype"
	top2_hap = f"{args.HAPLOTYPE.split('/')[-1]}_{args.CHROM}_top2_haplotype"
	if args.DATA in ['WGS', 'HiC']:
		if args.DATA == 'HiC':
			# top1
			Sarah(top1_hap, args.BAM, args.CHROM, args.NAME, args.SCORE, args.DATA, bamtools, args.THREADS)
			get_chrom_mate(gatk, args.BAM, top1_hap, args.CHROM, args.NAME)

			if args.NUMBER == 2: # top2
				Sarah(top2_hap, args.BAM, args.CHROM, args.NAME, args.SCORE, args.DATA, bamtools, args.THREADS)
				get_chrom_mate(gatk, args.BAM, top2_hap, args.CHROM, args.NAME)
#			A_bam = final_bam('A', top1_hap, args.NAME)
#			B_bam = final_bam('B', top1_hap, args.NAME)
#			process1 = mp.Process(target=bam2cool, args=(A_bam, args.REFERENCE))
#			process2 = mp.Process(target=bam2cool, args=(B_bam, args.REFERENCE))
#			process1.start()
#			process2.start()
#			for process in mp.active_children():
#				process.join()
#			A_mcool = A_bam.split('/')[-1].split('.bam')[0] + '.mcool'
#			B_mcool = B_bam.split('/')[-1].split('.bam')[0] + '.mcool'
#			cutoff = 0.5
#			sv_call(A_mcool, cutoff)
#			sv_call(B_mcool, cutoff)
	




