'''
Takes a file of positions (boundaries etc) in UCSC chr:pos1-pos2 format, calculates the distance to the nearest boundary for each 
individual element
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()

	parser.add_option("-f", "--file", dest="filename",
					  help="files with positions in chr:pos1-pos2 format", metavar="FILE")

	

	(options, args) = parser.parse_args()
	return options

options = parse_options()
elements = {}

infile = open(options.filename, 'r')
out_stem = re.sub('.txt', '', options.filename)
outfile = open(out_stem + '_distToNearest.txt','w')

for line in infile:
	line = line.rstrip()
	elements[line] = 1

for line1 in elements:
	(chr1, positions) = line1.split(':')
	(pos1L, pos1R) = positions.split('-')
	min_dist = 100000000
	for line2 in elements:
		if (line1 != line2):
			(chr2, positions) = line2.split(':')
			(pos2L, pos2R) = positions.split('-')
			dist = abs(int(pos2L) - int(pos1L))
			if (chr2 == chr1 and dist < min_dist):
				min_dist = dist

	if (min_dist < 100000000):
		outfile.write(line1 + '\t' + str(min_dist) + '\n')



infile.close()
outfile.close()

