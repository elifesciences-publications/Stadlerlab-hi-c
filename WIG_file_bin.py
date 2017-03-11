"""

"""

from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="WIG file, with 4 columns", metavar="FILE")
	
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
	
	parser.add_option("-o", "--outfile_stem",
						dest="outfile_stem", default='none',
					  help="output file stem")

	(options, args) = parser.parse_args()
	return options

def add_line(chr, bin, value, genome):
	if (chr in genome):
		if (bin in genome[chr]):
			genome[chr][bin] += value
		else:
			genome[chr][bin] = value
	else:
		genome[chr] = {}
		genome[chr][bin] = value	


def update_max_bin(chr, bin, max_bin):
	if (chr in max_bin):
		if(bin > max_bin[chr]):
			max_bin[chr] = bin
	else:
		max_bin[chr] = bin

options = parse_options()
bin_size = int(options.bin_size)
genome = {}
max_bin = {}

infile = open(options.filename, 'r')

for line in infile:
	line = line.rstrip()
	if (not re.match('track', line)):
		(chr, posL, posR, value) = line.split()
		bin = int(int(posL) / bin_size)
		add_line(chr, bin, float(value), genome)
		update_max_bin(chr, bin, max_bin)

infile.close()

outfile = open(options.outfile_stem + '_binned' + str(int(bin_size / 1000)) + 'kb.txt','w')

for chr in genome.keys():
	for i in range(0, max_bin[chr] + 1):
		outfile.write(chr + '_' + str(i) + '\t')
		if (i in genome[chr]):
			outfile.write(str(genome[chr][i]))
		else:
			outfile.write('0')
		outfile.write('\n')

outfile.close()

