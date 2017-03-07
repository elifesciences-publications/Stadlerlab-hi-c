'''

'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-d", "--dnase_file", dest="dnase_file",
					  help="WIG file of DNase accessibility", metavar="DNASEFILE")
	parser.add_option("-c", "--chromosome", dest="chromosome",
					  help="chromosome", metavar="CHROMOSOME")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100000,
					  help="bin size bp", metavar="BINSIZE")
	parser.add_option("-i", "--individual", dest="individual", default=False,
					  help="Print individual lines? (any character)", metavar="INDIV")
	parser.add_option("-g", "--gff", dest="gff", default=False,
					  help="Is position file in gff format?", metavar="GFF")

	(options, args) = parser.parse_args()
	return options


options = parse_options()




