'''
Stupid little script that takes an existing WIG file (presumably DNase) and just makes it uniform
by setting all values to 1. I wrote it with a mind to doing my own observed/expected matrices but
it didn't work well for some reason. I'll have to look into that, but in the meantime, I'll just
use juicebox for that.
'''
from optparse import OptionParser
import sys
import re
import random
from random import randint
from random import shuffle
from numpy import percentile

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="file",
					  help="JuiceBox file matrix dump", metavar="FILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile", metavar="FILE")
	(options, args) = parser.parse_args()
	return options


options = parse_options()

infile = open(options.file, 'r')

for line in infile:
	line = line.rstrip()
	if (not re.match('track', line)):
		items = line.split()
		items[3] = '1'
		s = '\t'
		line = s.join(items)
	print(line)	

	

infile.close()