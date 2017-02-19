"""

"""

from optparse import OptionParser
import sys

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--file", dest="filename",
					  help="paired alignment FILE", metavar="FILE")
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
	parser.add_option("-w", "--window_size",
					  dest="window",
					  help="window size")		
					    
	(options, args) = parser.parse_args()
	return options

def Add_read(bin, container, chr):
	if (container.has_key(chr)):
		if (container[chr].has_key(bin)):
			container[chr][bin] = container[chr][bin] + 1
		else:
			container[chr][bin] = 1
	else:
		container[chr] = {}
		container[chr][bin] = 1
	
	if(max_bin.has_key(chr)):
		if (bin > max_bin[chr]):
			max_bin[chr] = bin
	else:
		max_bin[chr] = bin

def PrintArray(container, filename):
	fileout = open(filename, 'w')
	fileout.write('track type=wiggle_0 name="test"\n')
	for chr in ('X', '2L', '2R', '3L', '3R'):
		fileout.write('variableStep chrom=chr' + chr + '\n')
		for i in range (0, max_bin[chr]):
			coordinate = str((i * bin_size) + (bin_size / 2))
			if (container[chr].has_key(i)):
				 fileout.write(coordinate + '\t' + str(container[chr][i]))
			else:
				fileout.write(coordinate + '\t0')
			fileout.write('\n')

options = parse_options()
bin_size = int(options.bin_size)
window = int(options.window)
left_counts = {}
right_counts = {}
max_bin = {}
left_file = 'test_left.txt'
right_file = 'test_right.txt'

file1 = open(options.filename, 'r')

for line in file1:
	line = line.rstrip()
	items = line.split()
	(chr1, chr2) = (items[2],items[5],)
	Lmost1 = int(items[3])
	Lmost2 = int(items[6])
	if (chr1 == chr2 and abs(Lmost1 - Lmost2) <= window):
		bin1 = int(Lmost1 / bin_size)
		bin2 = int(Lmost2 / bin_size)
		if (bin1 < bin2):
			Add_read(bin1, right_counts, chr1)
			Add_read(bin2, left_counts, chr1)
		if (bin1 > bin2):	
			Add_read(bin1, left_counts, chr1)	
			Add_read(bin2, right_counts, chr1)

PrintArray(left_counts,left_file)
PrintArray(right_counts,right_file)