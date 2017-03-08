'''

'''
from optparse import OptionParser
import sys
import re
import random
from random import randint
from numpy import percentile

def parse_options():
	parser = OptionParser()
	parser.add_option("-d", "--dnase_file", dest="dnase_file",
					  help="WIG file of DNase accessibility", metavar="DNASEFILE")
	parser.add_option("-p", "--distance_file", dest="distance_file",
					  help="file with distance probabilities", metavar="DISTFILE")
	parser.add_option("-c", "--chromosome", dest="chromosome",
					  help="chromosome", metavar="CHROMOSOME")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100000,
					  help="bin size bp", metavar="BINSIZE")
	parser.add_option("-a", "--a", dest="a", default=50,
					  help="a, free ends = (ad)^x", metavar="A")
	parser.add_option("-x", "--x", dest="x", default=0.8,
					  help="x, free ends = (ad)^x", metavar="X")
	parser.add_option("-i", "--iterations", dest="iterations", default=10000000,
					  help="number of times to iterate", metavar="ITERATIONS")
	parser.add_option("-o", "--outfile", dest="outfile", default='HiC_simulation.txt',
					  help="outfile", metavar="OUTFILE")

	(options, args) = parser.parse_args()
	return options

# Read a file of probabilities of seeing a read at various bin distances. Returns an array with number of 
# entries proportional to probability, so randomly selecting member of this array returns 
def read_distances(file):
	distances = []
	expansion_factor = 10000 #multiply the probability by this, take int value. Means anything with p < 1/expansion factor will be 0, also array will be ~ expansion_factor long (aprox bc of rounding)
	dist_file = open(file, 'r')
	for line in dist_file:
		line = line.rstrip()
		(dist, prob) = line.split()
		prob = float(prob)
		dist = int(dist)
		iterations = int(prob * expansion_factor)
		distances = distances + ([dist] * iterations) #this ([n] * i) notation makes an array of n repeated i times, the + is a concatenation of lists
	dist_file.close()
	return(distances)

def read_dnase(file, chr, bin_size):
	dnase = {}
	max_bin = 0
	dnase_file = open(file, 'r')
	for line in dnase_file:
		if (not re.match('track', line)):
			line = line.rstrip()
			(chr2, posL, posR, value) = line.split()
			posL = int(posL)
			value = float(value)
			if (chr2 == chr):
				bin = int(int(posL) / bin_size)
				if (bin in dnase):
					dnase[bin] += value
				else:
					dnase[bin] = value
				if (bin > max_bin):
					max_bin = bin
	dnase_file.close()
	return(dnase, max_bin)

def select_bin_from_dist(bin1, distances, max_bin):
	for i in range(0,1000): #this is really a while true, I jsut don't want to get stuck in an infinite loop during testing
		dist = random.choice(distances)
		direction = randint(0,1)
		if (direction == 0): #left
			bin2 = bin1 - dist
			if (bin2 > 0):
				return(bin2)
		if (direction == 1): #right
			bin2 = bin1 + dist
			if (bin2 <= max_bin):
				return(bin2)

def generate_linkages(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, a, x):
	# get free ends
	free_ends_bin1 = get_free_ends(bin1, dnase_values, dnase_norm_factor, a, x)
	free_ends_bin2 = get_free_ends(bin2, dnase_values, dnase_norm_factor, a, x)
	# randomly pair ends, reducing list each time
	pair_ends(bin1, bin2, free_ends_bin1, free_ends_bin2, count_matrix)
	return int((free_ends_bin1 + free_ends_bin2) / 2)

def get_free_ends(bin, dnase_values, dnase_norm_factor, a, x):
	norm_dnase = dnase_values[bin] / dnase_norm_factor
	num_free_ends = int(a * (norm_dnase ** x))
	if (num_free_ends >= 1): #if a region never has a free end the program crashes with a key error, plus this just doesn't make sense
		return (num_free_ends)
	else:
		return (1)

def pair_ends (bin1, bin2, free_ends_bin1, free_ends_bin2, count_matrix):
	ends = ([bin1] * free_ends_bin1) + ([bin2] * free_ends_bin2)
	for i in range(0, int(len(ends) / 2)): #correct number of iterations to use all ends, leave one dangling if necessary ::shrug::
		end1 = ends.pop(randint(0,len(ends) - 1)) 
		end2 = ends.pop(randint(0,len(ends) - 1))
		add_linkage(end1, end2, count_matrix)
		if (end1 != end2):
			add_linkage(end2, end1, count_matrix)
	
def add_linkage(end1, end2, count_matrix):
	if (end1 in count_matrix):
		if (end2 in count_matrix[end1]):
			count_matrix[end1][end2] += 1
		else:
			count_matrix[end1][end2] = 1
	else:
		count_matrix[end1] = {}
		count_matrix[end1][end2] = 1

def print_matrix(count_matrix, max_bin, outfilename):
	outfile = open(outfilename, 'w')
	for i in range (0, max_bin):
		for j in range (0, max_bin):
			if (j != 0):
				outfile.write('\t')
			if (j in count_matrix[i]):
				outfile.write(str(count_matrix[i][j]))
			else:
				outfile.write('0')
		outfile.write('\n')
	outfile.close()


options = parse_options()
a = int(options.a)
x = float(options.x)
bin_size = int(options.bin_size)
distance_array = read_distances(options.distance_file)
(dnase_values, max_bin) = read_dnase(options.dnase_file, options.chromosome, bin_size) 
dnase_norm_factor = percentile(list(dnase_values.values()), 95)
count_total = 0
iterations = int(options.iterations)
count_matrix = {}

for i in range(0, iterations):
	bin1 = randint(0, max_bin)
	bin2 = select_bin_from_dist(bin1, distance_array, max_bin)
	links_formed = generate_linkages(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, a, x)

print_matrix(count_matrix, max_bin, options.outfile)




