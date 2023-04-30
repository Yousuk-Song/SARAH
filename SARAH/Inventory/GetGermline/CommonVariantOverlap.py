
import sys
import pysam
import os

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


def GetVAF(samfile, chrom, pos, ref, alt, common_variant):
	is_del = 0
	is_ins = 0
	VAF_D = {}
	read_numb = 0
	if len(ref) > len(alt): # deletion -> '-' in alt
		alt = alt + '-'*(len(ref)-1)
		is_del = 1
	elif len(ref) < len(alt): # insertion -> 'I' in ref
		ref = ref + 'I'*(len(alt)-1)
		is_ins = 1
	for read in samfile.fetch(chrom, pos, pos+1):
		if avoid_soft_clipping(read, pos):
			query_pos = pos - read.pos 
			align_res = align_seq(read)
			query_seq = align_res[0]
			cigar_seq = align_res[1]
			query_pos += cigar_seq[:query_pos].count('I')
			if query_pos >= 0:
				base = ''
				read_numb += 1
				if is_ins or is_del:
					if query_pos + len(alt) <= len(query_seq):
						base = query_seq[query_pos:query_pos + len(alt)]
						base_cigar = cigar_seq[query_pos:query_pos + len(alt)]
						if is_ins:
							if base not in [alt, ref]:
								if 'I' not in base_cigar: #  no ins -> just base
									base = base[0]

							else: # found matched query indel, but confusing
								if 'I' not in base_cigar: #  no ins -> just base
									base = base[0]
					else:
						base = query_seq[query_pos:]
						base_cigar = cigar_seq[query_pos:]
						if is_ins:
							if base not in [alt[:len(base)], ref[:len(base)]]:
								if 'I' not in base_cigar: #  no ins -> just base
									base = base[0]
							else:
								if 'I' not in base_cigar: #  no ins -> just base
									base = base[0]
								else:
									if base == alt[:len(base)]:
										base = alt
									elif base == ref[:len(base)]:
										base = ref
						elif is_del:
							if base == alt[:len(base)]:
								base = alt
							elif base == ref[:len(base)]:
								base = ref
				else: # no indel
					base = query_seq[query_pos]
				if base != '':
					if base in VAF_D:
						VAF_D[base] += 1
					else:
						VAF_D[base] = 1
		del read
	return [VAF_D, read_numb]

def ratio_cal(count_list):
	count_alt1 = count_list[0]
	count_alt2 = count_list[1]
	total = count_alt1+count_alt2
	if count_alt1 > count_alt2:
		return count_alt1/total
	else: 
		return count_alt2/total


def VAF_cal(read_numb, ALT):
	return round(ALT/read_numb, 5)


def SomaticRate(D, alt, ref):
	if len(alt) < len(ref):
		alt = alt + '-'*(len(ref)-len(alt))
	elif len(ref) > len(alt):
		ref = ref + '-'*(len(alt)-len(ref))
	somatic_n = 0
	for base in D:
		if base not in [alt, ref]:
			somatic_n += D[base]
	return round(somatic_n/(D[alt]+somatic_n), 4)

def DataToHist(data, bin_numb):
	import matplotlib.pyplot as plt
	import scipy.stats as ss
	fig = plt.hist(data, bins=bin_numb)
	plt.title(f'{vcf.replace(".vcf", "")} VAF Distribution')
	plt.xlabel("Variant Allele Frequency (%)")
	plt.ylabel("Variant Counts")
	plt.savefig(vcf.replace(".vcf", f".hist.{bin_numb}.png"))

def overlap(bam, chrom, common_dir):
	ScriptDir = sys.path[0]
	vaf_threshold = 0.1

	common_variant = f'{common_dir}/variants1722_cmm_snv+indels.{chrom}.sort.txt'
	samfile = pysam.AlignmentFile(bam, "rb")
	open_variant = open(common_variant, 'r')

	fo = open(bam.split('/')[-1].replace('.bam', '.KRG1722.vcf'), 'w')
	for line in open(f'{ScriptDir}/bin/VcfTemplate.vcf', 'r'):
		fo.write(line)

	vcf_dict = {}
	n = 0
	for line in open_variant:
		line = line.rstrip()
		if not n:
			n = 1
			continue
		col = line.split('\t')
		(chrom, pos, ref) = (col[0], int(col[1]) - 1, col[4])
		alt_dict = dict([[i.strip().split(':')[0], float(i.strip().split(':')[1])] for i in col[6].split(',')[:-1]])
		if len(alt_dict) == 1:
			alt = list(alt_dict.keys())[0]
			[VAF_D, read_count] = GetVAF(samfile, chrom, pos, ref, alt, common_variant)
			# remove 'I' for ins's ref
			if ref[1:] == 'I'*(len(ref)-1):
				ref = ref[0]
			if alt[1:] == 'I'*(len(alt)-1):
				alt = alt[0]

			# add '-' for del's alt
			if len(ref) > len(alt):
				if ref[0] == alt:
					alt += '-'*(len(ref)-1)

			(ref_count, alt_count) = (VAF_D.get(ref, -1), VAF_D.get(alt, -1))
			# remove '-' for del's alt
			if alt[1:] == '-'*(len(alt)-1):
				alt = alt[0]
			if len(VAF_D) == 1: # Homo
				pass
			else: # Hetero
				if -1 not in [ref_count, alt_count]:
					vcf_pos = pos+1
					vaf = VAF_cal(read_count, alt_count)
					if vaf_threshold < vaf < 1 - vaf_threshold: # filter out Homo or LOH
#						print('Hetero')
#						print(chrom, vcf_pos, ref, alt)
#						print(VAF_D)
#						print(f'ref : {ref_count}, alt : {alt_count}')
#						print(f'Read counts : {read_count}')
#						print(f'VAF : {vaf}')
#						print()
						vcf_line_list = [chrom, str(vcf_pos), '.', ref, alt, '.', '.', \
						f'AC={alt_count};AF={vaf};AN=2;DP={read_count}', 'GT', '0/1']
						vcf_line = '\t'.join(vcf_line_list) + '\n'
						if (chrom, pos) not in vcf_dict:
							vcf_dict[(chrom, pos)] = [vcf_line]
						else:
							vcf_dict[(chrom, pos)] += [vcf_line]
						del vcf_line_list
						del vcf_line
		del line

	for key in vcf_dict:
		if len(vcf_dict[key]) == 1:
			fo.write(vcf_dict[key][0])
	print(bam.split('/')[-1].replace('.bam', '.KRG1722.vcf'), 'finished')
	fo.close()








