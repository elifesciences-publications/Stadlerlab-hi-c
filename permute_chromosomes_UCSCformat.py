'''
Permutes the chromosomes for a file of positions in UCSC format (chr:pos1-pos2). Ensures it doesn't pick the same chromosome as
original. It's a quick and dirty way to make matched control regions. Works generally because the various chromosome arms are 
similarly sized so random positions probably exist on other chromosomes. A more sophisticated approach is necessary for anything 
more than this. Just quick and dirty.
'''
from optparse import OptionParser
from math import log
import sys
import re
import random

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="positions file in UCSC chr:pos1-pos2 format", metavar="FILE")

	(options, args) = parser.parse_args()
	return options

def get_rand_chr(chr):
	chromosomes = ['chrX','chr2L','chr2R','chr3L','chr3R','chr4']
	for i in range(1,100):
		new_chr = chromosomes[random.randint(0,4)]
		if (new_chr != chr):
			return(new_chr)



options = parse_options()

infile = open(options.filename, 'r')

for line in infile:
	line = line.rstrip()
	(chr, rest) = line.split(':')
	new_chr = get_rand_chr(chr)
	print(new_chr + ':' + rest)

infile.close()