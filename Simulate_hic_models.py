'''
This script simulates data for two different models of how Hi-C works. The model assumes that all loci are uniformly
distributed throughout the nucleus, and the probability of two loci (regions defined by bin_size b) contacting each
other is governed only by the linear chromosomal distance between them. 

Both models begin by repeatedly drawing two regions from the genome. The first region is drawn at random, the second is drawn
based on the empirically-derived (from Hi-C data...I get that this is circular but it really shouldn't matter too much)
probability of observing Hi-C linkages between two regions at various distances. This drawing simulates finding two regions 
adjacent in the nucleus, and the two models then determine the number and composition of Hi-C products that form as a
result of this juxtaposition.

In the standard model, a probability of "visibility" is next calculated for each region. This captures the likelihood
that the region will be accessible to nuclease digestion and ligation. Here, we are modeling this as a function of DNase
accessibility. The probability that a region is visible is determined as the normalized (0 to 1) DNase value raised to 
the power Z (model parameter supplied by user). This power serves to either compress (exponents less than 1) or expand
(exponents greater than 1) the dependence of visibility on DNase signal. For example, taking two regions with DNase values
of region A=0.8 and region B=0.4, for z=1 A is twice as likely to be visible as B, for z=2 it is four times more likely,
for z=0.5 it is ~1.4 times more likely. Once the individual probabilities of visibility are calculated, the probability
of forming a product is determined by multiplying these two probabilities together (assumption of independence). A single
linkage between the two regions is then either formed or not formed based on a biased coin flip according to this total
probability.

In the Eisen model, DNase signal is used to generate a discrete number of "free ends" for each region. The conversion of 
normalized (0 to 1) DNase signal to free ends is free_ends = a * (DNase)^x [sorry for the variable names...b was taken up]
a and x are user supplied parameters. x has the same effect as Z in the standard model, to either compress or expand the 
sensitivity of free end count to DNase signal. Once the number of free ends of each region are determined, pairs are then
formed randomly. All pairings are allowed: A-A, B-B, and A-B. The counts for each are incorporated into the count matrix.

There are a bunch more notes I need to add, but that's a decent description.
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
	parser.add_option("-d", "--dnase_file", dest="dnase_file",
					  help="WIG file of DNase accessibility", metavar="DNASEFILE")
	parser.add_option("-p", "--distance_file", dest="distance_file",
					  help="file with distance probabilities", metavar="DISTFILE")
	parser.add_option("-c", "--chromosome", dest="chromosome",
					  help="chromosome to simulate", metavar="CHROMOSOME")
	parser.add_option("-b", "--bin_size", dest="bin_size", default=100000,
					  help="bin size in bp", metavar="BINSIZE")
	parser.add_option("-a", "--a", dest="a", default=50,
					  help="a, free ends = a * (DNase)^x", metavar="A")
	parser.add_option("-x", "--x", dest="x", default=0.8,
					  help="x, free ends = a * (DNase)^x", metavar="X")
	parser.add_option("-i", "--iterations", dest="iterations", default=10000000,
					  help="number of times to iterate", metavar="ITERATIONS")
	parser.add_option("-o", "--outfile", dest="outfile", default='HiC_simulation.txt',
					  help="outfile", metavar="OUTFILE")
	parser.add_option("-s", "--standard_model", dest="standard_model", default=False,
					  help="use standard model? (any nonzero value)", metavar="STDMOD")
	parser.add_option("-z", "--z", dest="z", default=0.5,
					  help="Exponent for converting normalized DNase score to prob. of visibility prob(visibility) = (DNase)^z [standard model only]", metavar="Z")

	(options, args) = parser.parse_args()
	return options

# Read a file of probabilities of seeing a read at various bin distances. Returns an array with number of 
# entries proportional to probability, so randomly selecting member of the array selects a bin according to 
# distance distribution
def read_distances(file):
	distances = []
	expansion_factor = 100000 #multiply the probability by this, take int value. Means anything with p < 1/expansion factor will be 0, also array will be ~ expansion_factor long (aprox bc of rounding)
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

# Reads WIG file of DNase accessibility (could be any WIG file). Sums values in each bin defined by -b. Also finds the max bin encountered.
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

# For a given bin, draws a distance from the distance distribution. Then randomly picks left or right direction. Checks to
# make sure the bin at that distance and direction is on the chromosome, returns the new bin if so, iterates until it finds
# one if not. No reason for this to fail, but I added failure mode with reporting
def select_bin_from_dist(bin1, distances, max_bin):
	for i in range(0,1000): #this is really a while true, I just don't want to get stuck in an infinite loop during testing
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
	print('Failed to find a suitable partner for ' + str(bin1))
	return(False)

# Mostly calls other things. Gets num free ends for both, randomly pairs them, returns the number of links formed
def generate_linkages(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, a, x):
	# get free ends
	free_ends_bin1 = get_free_ends(bin1, dnase_values, dnase_norm_factor, a, x)
	free_ends_bin2 = get_free_ends(bin2, dnase_values, dnase_norm_factor, a, x)
	# randomly pair ends, reducing list each time
	pair_ends(bin1, bin2, free_ends_bin1, free_ends_bin2, count_matrix)
	return int((free_ends_bin1 + free_ends_bin2) / 2)

# Calculates the number of free ends for a given bin based on scaled DNase data. Currently this is deterministic, but it could be
# made random by modeling free end count as poisson process lambda = num_free_ends below.
def get_free_ends(bin, dnase_values, dnase_norm_factor, a, x):
	norm_dnase = dnase_values[bin] / dnase_norm_factor
	num_free_ends = int(a * (norm_dnase ** x))
	if (num_free_ends >= 1): #if a region never has a free end the program crashes with a key error, plus this just doesn't make sense
		return (num_free_ends)
	else:
		return (1)

# Builds an array with entries of bin1 and bin2, with the number of each matching the calculated number of free ends. Randomly shuffles the array,
# then repeatedly extracts pairs of array entries, records them as a linkage in count matrix. Iterates until array is empty (or singlton left)
# This is assumption of ligation driven to completion, but could probably also tweak this by modeling the number of iterations in some
# random way if we are concerned about that assumption.
def pair_ends (bin1, bin2, free_ends_bin1, free_ends_bin2, count_matrix):
	ends = ([bin1] * free_ends_bin1) + ([bin2] * free_ends_bin2)
	shuffle(ends) #surprising syntax. Doesn't return the shuffled list, rather shuffles it in place
	for i in range(0, int(len(ends) / 2)): #correct number of iterations to use all ends, leave one dangling if necessary ::shrug::
		end1 = ends.pop() 
		end2 = ends.pop()
		add_linkage(end1, end2, count_matrix)
		if (end1 != end2):
			add_linkage(end2, end1, count_matrix)

# Just handles hash checks to add a linkage to count matrix	
def add_linkage(end1, end2, count_matrix):
	if (end1 in count_matrix):
		if (end2 in count_matrix[end1]):
			count_matrix[end1][end2] += 1
		else:
			count_matrix[end1][end2] = 1
	else:
		count_matrix[end1] = {}
		count_matrix[end1][end2] = 1

# Gets probability of visibility of two bins, calculates the overall prob of linkage, flips the biased coin
# and then adds (or doesn't) the linkage to the count matrix
def generate_linkage_standardmodel(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, z):
	prob_bin1 = get_visible_prob_standardmodel(bin1, dnase_values, dnase_norm_factor, z)
	prob_bin2 = get_visible_prob_standardmodel(bin2, dnase_values, dnase_norm_factor, z)
	prob_combined = prob_bin1 * prob_bin2
	success = flip_coin(prob_combined)
	#print(str(prob_combined) + ' ' + str(success))
	if(success):
		add_linkage(bin1, bin2, count_matrix)
		if (bin1 != bin2):
			add_linkage(bin2, bin1, count_matrix)

# Flipping a biased coin according to prob of success p.
def flip_coin(p):
	return True if random.random() < p else False

# Calc visibility of a bin by scaling DNase signal according to z
def get_visible_prob_standardmodel(bin, dnase_values, dnase_norm_factor, z):
	norm_dnase = dnase_values[bin] / dnase_norm_factor
	prob = norm_dnase ** z
	if(prob > 1):
		return(1)
	return(prob)

# Prints the count matrix
def print_matrix(count_matrix, max_bin, outfilename, chr):
	outfile = open(outfilename, 'w')

	for i in range(0, max_bin):
		if (i != 0):
			outfile.write('\t')
		outfile.write(chr + '_' + str(i))
	outfile.write('\n')
	for i in range (0, max_bin):
		outfile.write(chr + '_' + str(i))
		for j in range (0, max_bin):
		#	if (j != 0):
			outfile.write('\t')
			if (j in count_matrix[i]):
				outfile.write(str(count_matrix[i][j]))
			else:
				outfile.write('0')
		outfile.write('\n')
	outfile.close()

# MAIN 
options = parse_options()
a = int(options.a)
x = float(options.x)
chromosome = options.chromosome
if (not re.match('chr', chromosome)):
	chromosome = 'chr' + chromosome
bin_size = int(options.bin_size)
distance_array = read_distances(options.distance_file)
(dnase_values, max_bin) = read_dnase(options.dnase_file, chromosome, bin_size) 
dnase_norm_factor = percentile(list(dnase_values.values()), 95)
count_total = 0
iterations = int(options.iterations)
count_matrix = {}
standard_model = options.standard_model
z = float(options.z)

for i in range(0, iterations):
	if (i % 10000000 == 0):
		print(i)
	bin1 = randint(0, max_bin)
	bin2 = select_bin_from_dist(bin1, distance_array, max_bin)
	if(bin2):
		if (standard_model):
			generate_linkage_standardmodel(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, z)
		else:	
			links_formed = generate_linkages(bin1, bin2, dnase_values, dnase_norm_factor, count_matrix, a, x)
		
print_matrix(count_matrix, max_bin, options.outfile, chromosome)




