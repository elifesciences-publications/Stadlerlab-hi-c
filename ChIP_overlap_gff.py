'''
Takes two types of files. First is a list of boundary elements (or, hell, anything else), the second is a series of gff3 files
with regions of ChIP enrichment (or, again, anything else). Given these regions (and user-defined flanks), simply marks a region
as 0 (does not have a gff "peak" within supplied boundaries) or 1 (overlaps with at least 1 gff peak). Prints a matrix

boundary file: UCSC format of chr:pos1-pos2
Written to work with modencode files, not sure if there are variable format issues for other sources.
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
	parser.add_option("-n", "--names", dest="track_names",
					  help="names of tracks, comma delim", metavar="NAMES")

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
	items = line.split()
	if (len(items) == 9):
		chr = items[0]
		if (not re.match('chr',chr)):
			chr = 'chr' + chr
		pos1 = int(items[3])
		pos2 = int(items[4])

		for i in range(pos1, pos2):
			if (chr in genome):
				if(i in genome[chr]):
					hits[genome[chr][i]] = 1
		for hit in hits:
			counts[hit][track_name] = 1
			

def add_value(counts, hit, value):
	pass

def add_track(track_name, counts):
	for item in counts:
		counts[item][track_name] = 0

options = parse_options()
genome = {}
counts = {}
extension3 = int(options.three_extend)
extension5 = int(options.five_extend)

pos_file = open(options.pos_file, 'r')

for line in pos_file:
	add_position(line, genome, counts, extension5, extension3)

pos_file.close()

feature_file_list = options.features_files
feature_files = feature_file_list.split(',')
track_names = options.track_names.split(',')

for i in range(0,len(feature_files)):
	feature_file = feature_files[i]
	track_name = track_names[i]
	add_track(track_name, counts)
	file2 = open(feature_file, 'r')
	for line in file2:
		line = line.rstrip()
		add_feature_count(line, genome, counts, track_name)

	file2.close()

outfile = open(options.outfile,'w')

outfile.write(track_names[0]) #write first item of tracks
for track in track_names[1:]:
	outfile.write('\t' + track)
outfile.write('\n')

for item in counts:
	outfile.write(item)
	for track in track_names:
		outfile.write('\t' + str(counts[item][track]))
	outfile.write('\n')

outfile.close()




