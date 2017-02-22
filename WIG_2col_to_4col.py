'''
Converts a two-column (fixed span) WIG file to a 4-column WIG format
'''
from optparse import OptionParser
import sys
import re
from Bio import SeqIO

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="wig_file",
					  help="WIG file", metavar="FILE")

	(options, args) = parser.parse_args()
	return options

options = parse_options()

infile = open(options.wig_file,'r')
outfilename = ''
#print(options.wig_file)
if re.search('WIG',options.wig_file):
	outfilename = re.sub('.WIG','_4col.WIG', options.wig_file)
elif re.search('wig',options.wig_file):
	outfilename = re.sub('.wig','_4col.WIG', options.wig_file)
else:
	sys.exit('Does not appear to be a WIG file!')
outfile = open(outfilename, 'w')

curr_chr = ''
curr_span = 0
track_line_written = False

#outfile.write('track type=wiggle_0 name="test" description="test"\n')

for line in infile:
	line = line.rstrip()
	if (line[0:5] == 'track'):
		if(not track_line_written):
			outfile.write(line + '\n')
			track_line_written = True
	elif(line[0:4] == 'vari'):
		(blank, chr, span) = line.split()
		curr_chr = chr[6:]
		if (not re.match('chr', curr_chr)):
			curr_chr = 'chr' + curr_chr
		curr_span = int(span[5:])
		#print(curr_chr + ' ' + curr_span)
	else:
		#print(line)
		(pos, value) = line.split()
		pos = int(pos)
		outfile.write(curr_chr + '\t' + str(pos) + '\t' + str(pos + curr_span - 1) + '\t' + value + '\n')
infile.close()