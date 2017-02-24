'''
Takes a file of positions (boundaries etc) and a gff of features and calculates the distance from the 
leftmost location of the position to the nearest feature, reporting it if it's within some window.
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
	parser.add_option("-w", "--window", dest="window", default=0,
					  help="window to search in bp", metavar="WINDOW")

	

	(options, args) = parser.parse_args()
	return options

def calc_position(line, genome, window, outfile):
	#outfile = open(outfilename, 'w')
	(chr, positions) = line.split(':')
	(pos1, pos2) = positions.split('-')
	boundary_Lmost = int(pos1)
	dist_to_nearest_L = search_left(genome, chr, boundary_Lmost, window)
	dist_to_nearest_R = search_right(genome, chr, boundary_Lmost, window)

	if (dist_to_nearest_L < window + 1):
		if (dist_to_nearest_L <= dist_to_nearest_R):
			outfile.write(line + '\t' + str(dist_to_nearest_L) + '\n')
			#return(dist_to_nearest_L)
	if (dist_to_nearest_R < window + 1):
		if (dist_to_nearest_R < dist_to_nearest_L):
			outfile.write(line + '\t' + str(dist_to_nearest_R) + '\n')
			#return(dist_to_nearest_R)
	#return('NA')
		
	#outfile.close()

def search_left(genome, chr, boundary_Lmost, window):
	for i in range(0,window):
		test_pos = boundary_Lmost - i
		if(test_pos in genome[chr]):
			return(i)
	return(window + 1)

def search_right(genome, chr, boundary_Lmost, window):
	for i in range(0,window):
		test_pos = boundary_Lmost + i
		if(test_pos in genome[chr]):
			return(i)
	return(window + 1)

def add_feature(line, genome):
	items = line.split()
	if (len(items) == 9):
		chr = items[0]
		#print(chr)
		if (not re.match('chr',chr)):
			chr = 'chr' + chr
		pos1 = int(items[3])
		pos2 = int(items[4])

		for i in range(pos1, pos2):
			if (chr in genome):
				genome[chr][i] = 1
			else:
				#print(chr)
				genome[chr] = {}
				genome[chr][i] = 1

# Main			

options = parse_options()
genome = {}
window = int(options.window)

# Read feature file, assign to genomic position
featurefile = open(options.features_file, 'r')

for line in featurefile:
	line = line.rstrip()
	add_feature(line, genome)

featurefile.close()

# Go through position (boundaries) file and calc distances
pos_file = open(options.pos_file, 'r')
outfile = open(options.outfile, 'w')
for line in pos_file:
	line = line.rstrip()
	calc_position(line, genome, window, outfile)
	####outfile.write('testing')
	#if (dist != 'NA'):
		#print('booya')
		#outfile.write(line + '\t' + str())

pos_file.close()
outfile.close()

'''
outfile = open(options.outfile,'w')

outfile.write(track_names[0]) #write first item of tracks
for track in track_names[1:]:
	outfile.write('\t' + track)
outfile.write('\n')

for item in counts:
	#print(item + '\t' + str(counts[item]['CTCF']))
	outfile.write(item)
	for track in track_names:
		#print(track)
		outfile.write('\t' + str(counts[item][track]))
	outfile.write('\n')

outfile.close()
'''



