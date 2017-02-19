"""
Script builds a binned contact matrix for only a subsection of a single chromosome. Allows
high-res (small bin) without blowing your computer up.

v1a: supports multiple comma-delimited input files
"""

from optparse import OptionParser
import sys

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--files", dest="filenames",
					  help="paired alignment FILES, sep. with comma", metavar="FILE")
	parser.add_option("-b", "--bin_size",
					   dest="bin_size", default=1000000,
					  help="bin size")
	parser.add_option("-c", "--chr",
					  dest="chromosome",
					  help="chromosome")
	parser.add_option("-s", "--start",
					  dest="start",
					  help="start")
	parser.add_option("-e", "--end",
					  dest="end",
					  help="end")				  

	(options, args) = parser.parse_args()
	return options

def Add_read(bin1, bin2):
	if (bin_bin_counts.has_key(bin1)):
		if (bin_bin_counts[bin1].has_key(bin2)):
			bin_bin_counts[bin1][bin2] += 1
		else:
			bin_bin_counts[bin1][bin2] = 1
	else:
		temp = {}
		temp[bin2] = 1
		bin_bin_counts[bin1] = temp

def write_matrix(bin1, bin2):
	#determine if bin-bin counts exist, print
	if (bin_bin_counts.has_key(bin1)):
		#test = bin_bin_counts[bin1]
		if (bin_bin_counts[bin1].has_key(bin2)):
			sys.stdout.write('\t' + str(bin_bin_counts[bin1][bin2]))
		else:
			sys.stdout.write("\t0")
	else:
		sys.stdout.write("\t0")

options = parse_options()
bin_size = int(options.bin_size)
chosen_chromosome = options.chromosome
start = int(options.start)
end = int(options.end)
bin_bin_counts = {}
max_bin = 0

files = options.filenames
files = files.split(',')

for file in files:
	file1 = open(file, 'r')
	for line in file1:
		line = line.rstrip()
		items = line.split()
		(chr1, chr2) = (items[2], items[5])
		if (chosen_chromosome == chr1 == chr2):
			(Lmost1, Lmost2) = (int(items[3]), int(items[6]))
			if (Lmost1 >= start and Lmost1 <= end and Lmost2 >= start and Lmost2 <= end):
				Lmost1 = Lmost1 - start
				Lmost2 = Lmost2 - start
				bin1 = int(Lmost1 / bin_size)
				bin2 = int(Lmost2 / bin_size)
				if (max(bin1, bin2) > max_bin):
					max_bin = max(bin1, bin2)
					#sys.stderr.write(str(max_bin) + '\n')
				if (bin1 <= bin2): 
					Add_read(bin1, bin2)
				else:
					Add_read(bin2, bin1)

#print options.filename
for x in range(0, max_bin):
	for y in range(0, x):
		write_matrix(y,x)
		#sys.stdout.write()
	if(x+1 <= max_bin):
		for y in range(x+1, max_bin):
			write_matrix(x,y)
	sys.stdout.write("\n")	