'''
Takes a HI-C paired map file and a bin size, calculate the probability of observing links between two bins
separated by n bins. Prints as distance [tab] probability
'''
from optparse import OptionParser
import sys
import re
from numpy import power
from numpy import random

def parse_options():
	parser = OptionParser()
	parser.add_option("-a", "--a", dest="a", default=2,
					  help="a of P ~ aX^a-1", metavar="A")


	(options, args) = parser.parse_args()
	return options


options = parse_options()
a = int(options.a)

for i in range(0,1000):
	value = float(random.power(a,1))
	inverse = 1 / value
	print(inverse)






