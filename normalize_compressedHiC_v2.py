'''
Takes a reduced Hi-C matrix file output by HiC_count_bin_linkages_singleChr_diagonal and 
normalizes based on the bin counts, vanilla coverage normalization.

'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	(options, args) = parser.parse_args()
	return options
            
			

options = parse_options()
bin_counts = {}
arb_constant = 10000000
infile = options.filename

file1 = open(infile, 'r')
file_stem = re.sub('.txt', '', infile)
outfile = open(file_stem + '_VCnorm.txt', 'w')


for line in file1:
	line = line.rstrip()
	if (line[0] == '#'):
		line = line[1:]
		(chr, bin, count) = line.split('\t')
		if (chr not in bin_counts):
			bin_counts[chr] = {}
		bin_counts[chr][bin] = int(count)
	
	else:
		(chr, bin1, bin2, count) = line.split('\t')
		outfile.write(chr + '\t' + bin1 + '\t' + bin2 + '\t')
		if (bin1 not in bin_counts[chr]):
			print (chr + '\t' + bin1)
		if (bin2 not in bin_counts[chr]):
			print (chr + '\t' + bin2)
		if (bin_counts[chr][bin1] != 0 and bin_counts[chr][bin2] != 0):
			count_norm = float(count) / bin_counts[chr][bin1] / bin_counts[chr][bin2] * arb_constant
			outfile.write(str(count_norm) + '\n')
		else:
			outfile.write('NA\n')

				
file1.close()
outfile.close()