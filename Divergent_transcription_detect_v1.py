'''

'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--boundary_file", dest="boundary_filename",
					  help="boundary_file", metavar="BOUNDFILE")
	parser.add_option("-t", "--transcript_file", dest="transcript_filename",
					  help="BED fiel of transcripts", metavar="TFILE")
	parser.add_option("-r", "--range", dest="range",
					  help="bp allowed between boundary and TSS", metavar="RANGE")
	parser.add_option("-c", "--range_close", dest="range_close",
					  help="bp required between boundary and TSS", metavar="RANGECLOSE")

	(options, args) = parser.parse_args()
	return options
            
options = parse_options()
chromosomes = ('chr2L', 'chr2R', 'chr3L', 'chr3R', 'chr4', 'chrX', 'chrM', 'chrU')

boundaries = {}
for i in chromosomes:
	boundaries[i] = {}

boundary_file = open(options.boundary_filename, 'r')
for line in boundary_file:
	line = line.rstrip()
	items = line.split()
	chr = items[0]
	posL = float(items[1])
	posR = float(items[2])
	middle = (posL + posR) / 2
	boundaries[chr][middle] = 1
boundary_file.close()

Lm = 0
Lp = 0
Rm = 0
Rp = 0
bp_range = int(options.range)
range_req = int(options.range_close)

transcript_file = open(options.transcript_filename, 'r')
for line in transcript_file:
	line = line.rstrip()
	items = line.split()
	chr = items[0]
	posL = float(items[1])
	posR = float(items[2])
	strand = items[5]
	
	if (strand == '-'): TSS = posR
	if (strand == '+'): TSS = posL

	for bound_pos in boundaries[chr]:
		if (abs(bound_pos - TSS) < bp_range and abs(bound_pos - TSS) > range_req):
			if (TSS < bound_pos):
				if (strand == '+'): Lp += 1
				else: Lm += 1
			else:
				if (strand == '+'): Rp += 1
				else: Rm += 1
transcript_file.close()

print ('Lm: ' + str(Lm) + '\n' + 'Lp: ' + str(Lp) + '\n' + 'Rm: ' + str(Rm) + '\n' + 'Rp: ' + str(Rp))



