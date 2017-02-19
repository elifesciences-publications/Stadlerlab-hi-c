'''
This is a naive approach to call boundaries in Hi-C data based on directionality--the relative numbers
contacts with loci to the right vs. left. It takes as input a WIG file with directionality assigned to 
bins. User supplies bin size (depends on input wig file), width (number of bins to left and right to 
include in window), and threshold.

The algorithm is to find the average directional bias (log10 of R/L) within the window of supplied width
to the left and right, with the bin of interest included in teh right-facing window. If both left and 
right are properly biased over the threshold, the bin is assigned a value of 1, otherwise it gets a 0. 
Output is a WIG file.

v2b: writes a separate file containing coordinates of boundaries in list form.
'''
from optparse import OptionParser
from math import log
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Directional bias file", metavar="FILE")
	parser.add_option("-b", "--binsize", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-w", "--width", dest="width",
					  help="width about diagonal, in bins", metavar="WIDTH")
	parser.add_option("-t", "--threshold", dest="threshold",
					  help="threshold for calling boundary", metavar="THRESH")
	parser.add_option("-n", "--track_name", dest="track_name",
					  help="name for browser track", metavar="TRACKNAME")

	(options, args) = parser.parse_args()
	return options

# Adds up counts within window specified
def sum_counts(chr, start, end, counts, bin_size):
	sum = 0
	for i in range(start, end + 1):
		bin = i * bin_size
		if (bin in counts[chr]):
			sum += counts[chr][bin]
	return(sum)

def update_max(max_bins, chr, bin_number):
	if (chr in max_bins):
		if(bin_number > max_bins[chr]):
			max_bins[chr] = bin_number
	else:
		max_bins[chr] = bin_number

def add_count(chr, bin, count, counts):
	if (chr in counts):
		counts[chr][bin] = count
	else:
		counts[chr] = {}
		counts[chr][bin] = count


options = parse_options()
bin_size = int(options.bin_size)
counts = {}
max_bins = {}
infile = open(options.filename,'r')
for line in infile:
	line = line.rstrip()
	if (line[0:5] != 'track'):
		items = line.split('\t')
		(chr, bin, end, count) = items
		bin = int(bin)
		count = float(count)
		bin_number = bin / bin_size # bit of confusion here betwee "bin number", which is 0, 1, 2, and the positions of the bins, 0, 500, 1000...
		update_max(max_bins, chr, bin_number)
		add_count(chr, bin, count, counts)
		
outfile_stem = re.sub('.WIG', '', options.filename)
outfile = open(outfile_stem + '_boundaries_w' + options.width + '_t' + options.threshold + '.WIG','w')
outfile2 = open(outfile_stem + '_boundaries_w' + options.width + '_t' + options.threshold + '.txt','w')
outfile.write('track type=wiggle_0 name="boundaries_' + options.track_name + '" description="boundaries_' + options.track_name + '"\n')
chromosomes = ('2L','X','3L','4','2R','3R')

for chr in chromosomes:
	chr = 'chr' + chr
	max_bin = int(max_bins[chr])
	width = int(options.width)
	inblock = False
	block_start = 0
	block_end = 0
	for bin1 in range(0, max_bin + 1):
		left_avg = sum_counts(chr, bin1 - width, bin1 - 1, counts, bin_size) / width
		right_avg = sum_counts(chr, bin1, bin1 + width - 1, counts, bin_size) / width
		right_thresh = float(options.threshold)
		left_thresh = -1 * right_thresh
		signal='0'
		if (left_avg <= left_thresh and right_avg >= right_thresh):
			signal = '1'
			if (inblock):
				block_end = bin1
			else:
				inblock = True
				block_start = bin1
				block_end = bin1
		else:
			if(inblock):
				start = str(block_start * bin_size)
				end = str((block_end * bin_size) + bin_size - 1)
				outfile2.write(chr + ':' + start + '-' + end + '\n')
			inblock = False
		bin_pos = bin1 * bin_size
		outfile.write(chr + '\t' + str(bin_pos) + '\t' + str(bin_pos + bin_size - 1) + '\t' + signal + '\n')
		

outfile.close()
outfile2.close()
infile.close()

