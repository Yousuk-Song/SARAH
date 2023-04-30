#!/data2/home/song7602/.conda/envs/py311/bin/python3.11

import sys
import os

vcf = sys.argv[1]
bam = sys.argv[2]

# 1. HAPCUT

hapcutfragment = vcf.replace('.vcf', '.fragment')
hapcutfile = vcf.replace('.vcf', '.haplotype')

order1 = ['extractHAIRS', '--hic 1', '--bam', bam, '--VCF', vcf, '--out', hapcutfragment, '--indels 1']
order2 = ['HAPCUT2', '--hic 1', '--fragments', hapcutfragment, '--VCF', vcf, '--output', hapcutfile]

os.system(' '.join(order1))
os.system(' '.join(order2))


