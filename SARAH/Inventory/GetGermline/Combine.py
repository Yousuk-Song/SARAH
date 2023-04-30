
def combine_vcf(tumor_hic_vcf, tumor_wgs_vcf, name1, name2): 
	open_tumor_hic_vcf = open(tumor_hic_vcf, 'r')
	open_tumor_wgs_vcf = open(tumor_wgs_vcf, 'r')
	
	combined_vcf = tumor_hic_vcf.replace('.Tumor.Hi-C', '.Tumor_Hi-C_and_Tumor_WGS')
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
			[ref, alt, info] = hic_D[(chrom, pos)]
			fo.write('\t'.join([chrom, str(pos), '.', ref, alt, '.', '.', info, 'GT', '0/1']) + '\n')
		del line
	open_tumor_wgs_vcf.close()
	fo.close()



			
			
	



