'''
Takes "diagonal bin coverage" from HiC_count_bin_linkages_singleChr_diagonal_v1a.py. Generates a series of 
Hi-C matrices of supplied width about diagonal, makes files for consecutive regions covering the entire
genomic region represented by the count file. This script is used to generate the panels for manual
boundary calling.

'''
from optparse import OptionParser
import sys
import re
import gzip

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-s", "--file_stem", dest="file_stem",
					  help="stem for outfiles", metavar="STEM")
	parser.add_option("-w", "--width", dest="width",
					  help="num bins on either side of diagonal to store, in bins", metavar="WIDTH")

	(options, args) = parser.parse_args()
	return options

def Generate_map(chr, bin1, bin2, bin_size):
	pos1 = bin1 * bin_size
	pos2 = bin2 * bin_size + bin_size - 1
	outfile_local = options.file_stem + '/' + chr + '_' + str(pos1) + '_' + str(pos2) + '.txt'
	outfile = open(outfile_local, 'w')
	#print colnames
	for k in range (bin1, bin2 + 1):
		if (k != bin1):
			outfile.write('\t')
		outfile.write(chr + '_' + str(k))
	outfile.write('\n')
		
	for i in range (bin1, bin2 + 1):
		outfile.write(chr + '_' + str(i)) #row name
		for j in range(bin1, bin2 + 1):
			count = '0.0'
			if (i in bin_counts[chr]):
				if(j in bin_counts[chr][i]):
					count = bin_counts[chr][i][j]
			outfile.write('\t' + count)
		outfile.write('\n')
			
	outfile.close()			

def add_count(bin_counts, chr, bin1, bin2, count):
	if (chr in bin_counts):
		if (bin1 in bin_counts[chr]):
			bin_counts[chr][bin1][bin2] = count
		else:
			bin_counts[chr][bin1] = {}
			bin_counts[chr][bin1][bin2] = count
	else:
		bin_counts[chr] = {}
		bin_counts[chr][bin1] = {}
		bin_counts[chr][bin1][bin2] = count

options = parse_options()
#bin_counts = MakeBins()
bin_counts = {}
width = int(options.width)
bin_size = int(options.bin_size)

f = options.filename
if (f[-2:] == 'gz'):
	infile = gzip.open(f, 'rt')
else:
	infile = open(options.filename,'r')

for line in infile:
	line = line.rstrip()
	if (line[0] != '#'):
		(chr, bin1, bin2, count) = line.split('\t')
		bin1 = int(bin1)
		bin2 = int(bin2)
		if (abs(bin1 - bin2) < (width + 2)): #just storing an extra couple bins to avoid thinking about end problems
			add_count(bin_counts, chr, bin1, bin2, count)
		#bin_counts[chr][bin1][bin2] = count

infile.close()


for chr in bin_counts:
	chr = re.sub('chr', '', chr)
	for i in range(0, 10000):
		start_bin = i * width
		end_bin = start_bin + width - 1
		if (end_bin in bin_counts[chr]):
			Generate_map(chr, start_bin, end_bin, bin_size)
		else:
			break

