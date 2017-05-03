'''

'''
from optparse import OptionParser
from math import log
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="Reduced bin file", metavar="FILE")
	parser.add_option("-b", "--bin_size", dest="bin_size",
					  help="Bin size in bp", metavar="BINSIZE")
	parser.add_option("-w", "--window_size", dest="window_size",
					  help="window size (on either side of home bin), in bins", metavar="WINDOW")
	parser.add_option("-s", "--step_size", dest="step_size",
					  help="step size, in bins", metavar="STEP")

	(options, args) = parser.parse_args()
	return options

# Initializes container based on the fly chromosome sizes (Dm3)            
def make_bins(chromosomes, sizes, bin_size):
	bin_counts = {}
	for i in range(0,6):
		num_bins = int(sizes[i] / bin_size)
		chr = chromosomes[i]
		bin_counts[chr] = {}
		for j in range(0, num_bins + 1):
			bin_counts[chr][j] = 0
	return bin_counts

# Reads reduced format into bin container
def read_boundary_scores(filename, chromosomes, sizes, bin_size):
	bin_counts = make_bins(chromosomes, sizes, bin_size)
	infile = open(options.filename,'r')
	for line in infile:
		line = line.rstrip()
		if (line[0:5] != 'track'): 
			(chr, pos1, pos2, count) = line.split('\t')
			bin_num = int(int(pos1) / bin_size)
			bin_counts[chr][bin_num] = float(count)
	infile.close()
	return bin_counts

# Runs through and scores each genomic position by it's left/right deal
def slide_window(bin_counts, sizes, chromosomes, window_size, step_size, bin_size):	
	window_scores = []
	for i in range(0,6):
		chr = chromosomes[i]
		max_bin = int(sizes[i] / bin_size)
		for bin_home in range(window_size, max_bin - window_size, step_size): #walking along chromosome
		#########################################
		# print some bin numbers to make sure we're not skipping anything
		##########################################
			max_value = -10000
			bin_of_max = 0
			for bin in range(bin_home - window_size, bin_home + window_size):
				if (bin_counts[chr][bin] > max_value):
					max_value = bin_counts[chr][bin]
					bin_of_max = bin
			if (bin_of_max != 0):
				pos_start = bin_of_max * bin_size
				pos_end = pos_start + bin_size - 1
				line = chr + '\t' + str(pos_start) + '\t' + str(pos_end) + '\t' + str(max_value)
				window_scores.append(line)
				#print(window_scores)
			else:
				print ('bad bin at ' + chr + ' ' + str(bin_home))
	return window_scores



def merge_print_window_scores(window_scores, merge_dist, outfile):
	skip_me = False
	for i in range (0, len(window_scores) - 1):
		if (not skip_me):
			line1 = window_scores[i].split()
			line2 = window_scores[i+1].split()
			chr1 = line1[0]
			pos1 = int(line1[1])
			chr2 = line2[0]
			pos2 = int(line2[1])
			skip_next = False
			#print(pos2)
			# logic: if they're not mergeable, print 1 and go to next iteration where position 2 will be 1. If they are mergable and 2 has the higher score,
			# do not print 1, go to next iteration where 2 is 1. If they are mergable and 1 has higher score, print 1 and delete 2, next iteration starts 
			# with +2 position
			if (chr1 == chr2):
				if (abs(pos1 - pos2) <= merge_dist):
					score1 = line1[3]
					score2 = line2[3]
					if (score1 > score2):
						outfile.write(window_scores[i] + '\n')
						skip_me = True
					#if (score2 > score1): do nothing
				else:
					outfile.write(window_scores[i] + '\n')
		else:
			skip_me = False


# Main section calling everything
options = parse_options()

sizes = (23011544, 22422827, 24543557, 1351857, 21146708, 27905053)
chromosomes = ('chr2L','chrX','chr3L','chr4','chr2R','chr3R')
window_size = int(options.window_size)
step_size = int(options.step_size)
bin_size = int(options.bin_size)
merge_dist = 2000

bin_counts = read_boundary_scores(options.filename, chromosomes, sizes, bin_size)
window_scores = slide_window(bin_counts, sizes, chromosomes, window_size, step_size, bin_size)
outfilename = re.sub('.WIG', '', options.filename) + 'windowMaximums_w' + str(window_size) + 's' + str(step_size) + '.WIG'
outfile = open(outfilename,'w')
name = 'windowMaximums_w' + str(window_size) + '_s' + str(step_size)
outfile.write ('track type=wiggle_0 name="' + name + '" description="' + name + 'test"' + '\n')
merge_print_window_scores(window_scores, merge_dist, outfile)
outfile.close()

