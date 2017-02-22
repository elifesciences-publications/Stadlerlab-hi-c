'''
Takes two types of files. First is a list of boundary elements (or, hell, anything else), the second is a series of WIG files
with scores for ChIP (or, again, anything else). Given these regions (and user-defined flanks), simply sums up the values for 
each ChIP track within this file, prints a matrix of positions x ChIP enrichment

boundary file: UCSC format of chr:pos1-pos2
WIG files must be 4-position: chr pos1 pos2 value
'''
from optparse import OptionParser
import sys
import re

def parse_options():
	parser = OptionParser()
	parser.add_option("-p", "--positions_file", dest="pos_file",
					  help="file with positions", metavar="POSFILE")
	parser.add_option("-f", "--features_file", dest="features_files",
					  help="files with features, comma-delimited", metavar="FEATFILE")
	parser.add_option("-o", "--outfile", dest="outfile",
					  help="outfile", metavar="OUTFILE")
	parser.add_option("-3", "--three_extension", dest="three_extend", default=0,
					  help="3' extension in bp", metavar="THREE")
	parser.add_option("-5", "--five_extension", dest="five_extend", default=0,
					  help="5' extension in bp", metavar="FIVE")

	(options, args) = parser.parse_args()
	return options

def add_position(line, genome, counts, extension5, extension3):
	line = line.rstrip()
	counts[line] = {}
	(chr, positions) = line.split(':')
	(pos1, pos2) = positions.split('-')
	start = int(pos1) - extension5
	end = int(pos2) + extension3
	if (chr not in genome):
		genome[chr] = {}
	for i in range(start, end):
		genome[chr][i] = line

def add_feature_count(line, genome, counts, track_name):
	hits = {}
	(chr, pos1, pos2, value) = line.split()
	pos1 = int(pos1)
	pos2 = int(pos2)
	value = float(value)

	for i in range(pos1, pos2):
		if (chr in genome):
			if(i in genome[chr]):
				hits[genome[chr][i]] = 1
	for hit in hits:
		#add_value(counts, hit, value)
		counts[hit][track_name] += value
			

def add_value(counts, hit, value):
	pass

def add_track(track_name, tracks, counts):
	for item in counts:
		counts[item][track_name] = 0
	tracks.append(track_name)

options = parse_options()
genome = {}
counts = {}
tracks = []
extension3 = int(options.three_extend)
extension5 = int(options.five_extend)

pos_file = open(options.pos_file, 'r')

for line in pos_file:
	add_position(line, genome, counts, extension5, extension3)

pos_file.close()

feature_file_list = options.features_files
feature_files = feature_file_list.split(',')

for feature_file in feature_files:
	file2 = open(feature_file, 'r')
	track_name = ''
	for line in file2:
		line = line.rstrip()
		if (re.search('track', line)):
			track_name = re.search('name="(.*?)"', line).group(1)
			add_track(track_name, tracks, counts)
		else:
			add_feature_count(line, genome, counts, track_name)

	file2.close()

outfile = open(options.outfile,'w')

outfile.write(tracks[0]) #write first item of tracks
for track in tracks[1:]:
	outfile.write('\t' + track)
outfile.write('\n')

for item in counts:
	outfile.write(item)
	for track in tracks:
		outfile.write('\t' + str(counts[item][track]))
	outfile.write('\n')

outfile.close()




