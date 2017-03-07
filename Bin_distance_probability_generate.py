'''

'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="file",
					  help="paired map file", metavar="FILE")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100000,
					  help="bin size in bp", metavar="BINSIZE")
	parser.add_option("-o", "--outfile", dest="outfile", 
					  help="outfile", metavar="OUTFILE")

	(options, args) = parser.parse_args()
	return options


options = parse_options()
bin_size = int(options.bin_size)
bin_distance_counts = {}
counts_total = 0
bin_distance_max = 0

infile = open(options.file, 'r')

for line in infile:
	line = line.rstrip()
	items = line.split()
	chr1 = items[2]
	pos1 = int(items[3])
	chr2 = items[5]
	pos2 = int(items[6])
	bin1 = int(pos1 / bin_size)
	bin2 = int(pos2 / bin_size)
	if (chr1 == chr2 and bin1 != bin2):
		bin_diff = abs(bin1 - bin2)
		if (bin_diff > bin_distance_max):
			bin_distance_max = bin_diff
		counts_total += 1
		if (bin_diff in bin_distance_counts):
			bin_distance_counts[bin_diff] += 1
		else:
			bin_distance_counts[bin_diff] = 1

infile.close()

outfile = open(options.outfile, 'w')

for i in range(1,bin_distance_max + 1):
	outfile.write(str(i) + '\t')
	if (i in bin_distance_counts):
		outfile.write (str(bin_distance_counts[i] / counts_total))
	else:
		outfile.write('0')
	outfile.write ('\n')

outfile.close()






