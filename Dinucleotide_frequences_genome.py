'''

'''
from optparse import OptionParser
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

def parse_options():
	parser = OptionParser()
	parser.add_option("-f", "--fasta_file", dest="fastafile",
					  help="fasta file", metavar="FILE")
	parser.add_option("-o", "--out_stem", dest="out_stem",
					  help="outfile stem", metavar="STEM")
	parser.add_option("-b", "--bin_size", dest="bin_size",
					  help="bin size", metavar="BINSIZE")

	(options, args) = parser.parse_args()
	return options

def calc_dn_freqs(dinuc, dinuc_rc, outfilestem, genome, bin_size):
	outfile = open(outfilestem + '_dinucFreqs_' + str(bin_size) + 'bp_' + dinuc + '.WIG', 'w')
	outfile.write('track type=wiggle_0 name="' + dinuc + '_freq" description="' + dinuc + '"\n')
	for chr in genome:
		print(dinuc + ' ' + chr)
		for pos in range(0, int(len(genome[chr]) / bin_size)):
			start = pos * bin_size
			end = (pos * bin_size) + bin_size
			seq = str(genome[chr][start:end])
			dinuc_pattern = '(?=' + dinuc + ')'
			dinuc_rc_pattern = '(?=' + str(dinuc_rc) + ')' 
			count = len(re.findall(dinuc_pattern, seq)) + len(re.findall(dinuc_rc_pattern, seq))
			outfile.write('chr' + chr + '\t' + str(start) + '\t' + str(end) + '\t' + str(count) + '\n')


	outfile.close()

options = parse_options()
bin_size = int(options.bin_size)
genome = {}
chromosomes = ['X', '2L', '2R', '3L', '3R', '4']
for record in SeqIO.parse(options.fastafile, "fasta"):
	if (str(record.name) in chromosomes):
		genome[record.id] = record.seq

nts = ['A', 'C', 'T', 'G']

for nt1 in nts:
	for nt2 in nts:
		dinuc = nt1 + nt2
		dinuc_rc = Seq(dinuc,generic_dna)
		dinuc_rc = dinuc_rc.reverse_complement()
		#print (dinuc + dinuc_rc)
		calc_dn_freqs(dinuc, dinuc_rc, options.out_stem, genome, bin_size)



