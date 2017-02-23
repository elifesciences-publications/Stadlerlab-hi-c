'''

'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-p", "--positions_file", dest="pos_file",
					  help="file with positions", metavar="POSFILE")
	parser.add_option("-f", "--features_file", dest="features_file",
					  help="files with features, comma-delimited", metavar="FEATFILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile", metavar="OUTFILE")
	parser.add_option("-w", "--window", dest="window", default=20,
					  help="window in bins", metavar="WINDOW")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100,
					  help="bin size bp", metavar="BINSIZE")


	(options, args) = parser.parse_args()
	return options

def add_position(line, genome, window, bin_size, counts):
	line = line.rstrip()
	(chr, positions) = line.split(':')
	(pos1, pos2) = positions.split('-')
	bin_Lmost = int(int(pos1) / bin_size)
	
	for i in range(-1 * window, window):
		bin = bin_Lmost + i
		if (bin in genome[chr]):
			counts[i] += genome[chr][bin]
	

def add_feature(line, genome, bin_size):
	(chr, pos1, pos2, value) = line.split()
	pos1 = int(pos1)
	pos2 = int(pos2)
	value = float(value)
	bin = int(pos1 / bin_size)

	if (chr in genome):
		if (bin in genome[chr]):
			genome[chr][bin] += value
		else:
			genome[chr][bin] = value
	else:
		genome[chr] = {}
		genome[chr][bin] = value
			




options = parse_options()
genome = {}
counts = {}
window = int(options.window)
bin_size = int(options.bin_size)
for i in range(-1 * window, window):
	counts[i] = 0

# first load feature file WIG
featurefile = open(options.features_file, 'r')

for line in featurefile:
	line = line.rstrip()
	if (re.search('track', line)):
		pass
	else:
		add_feature(line, genome, bin_size)

featurefile.close()

# Next deal with positions file
pos_file = open(options.pos_file, 'r')

for line in pos_file:
	add_position(line, genome, window, bin_size, counts)

pos_file.close()



outfile = open(options.outfile,'w')

for i in range(-1 * window, window):
	outfile.write(str(i) + '\t' + str(counts[i]) + '\n')

outfile.close()




