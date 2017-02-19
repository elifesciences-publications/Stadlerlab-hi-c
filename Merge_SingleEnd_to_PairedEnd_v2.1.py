# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 11:52:44 2014
Takes two bowtie outputs from single end alignments of paired reads, matches
based on read name and prints paired end alignments in format:
read name \t R1_strand \t R1_chromosome \t R1_Lmost \t R2_strand \t R2_chromosome \t R2_Lmost

Also filters away reads that map to the same chromosome within a set distance, currently
0 bp, set with distance_minimum variable
@author: MStadler
"""
import sys
from optparse import OptionParser
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file1", dest="filename1",
					  help="First mapping file", metavar="FILE1")
	parser.add_option("-s", "--file2", dest="filename2",
					  help="Second mapping file", metavar="FILE2")

	(options, args) = parser.parse_args()
	return options
		
options = parse_options()

file1_data = {}
goodPair_count = 0
total = 0
distance_minimum = 0

outfilename = ''
if (re.search('R1', options.filename1)):
	outfilename = options.filename1
elif (re.search('R1', options.filename2)):
	outfilename = options.filename2
outfilename = re.sub('_R1_','_', outfilename)
outfilename = re.sub('.bowtie', '', outfilename) + '_pairMerged_' + str(distance_minimum) + 'bp.txt'

outfile = open(outfilename, "w")
file1 = open(options.filename1, 'r')

for line in file1:
    pair_name = line.split()[0]
    splitline= line.split('\t')
    file1_data[pair_name] = "\t".join(splitline[1:4])
file1.close()

file2 = open(options.filename2, 'r')

for line in file2:
	total += 1
	pair_name = line.split()[0]
	splitline1= line.split('\t')
	if pair_name in file1_data:
		splitline2 = file1_data[pair_name].split('\t')
		chr1 = splitline1[2]
		pos1 = int(splitline1[3])
		chr2 = splitline2[1]
		pos2 = int(splitline2[2])
		if (chr1 == chr2 and abs(pos1 - pos2) < distance_minimum):
			pass
		else:
			outfile.write(pair_name + "\t" + file1_data[pair_name] + "\t" + "\t".join(splitline1[1:4]) + '\n')
			goodPair_count += 1

file2.close()


pct = float(goodPair_count) / float(total) * 100
# Optional stats reporting in python2.6 language.
#print >> sys.stderr, str(goodPair_count) + ' of ' + str(total) + ' reads have aligned pairs, or ' + str(pct) + ' percent'
