use warnings;
use strict;

while(<>){
	chomp;
	my @line = split/\t/;
	if ($line[5] eq '+'){
		$line[2] = $line[1] + 1;
	}
	if ($line[5] eq '-'){
		$line[1] = $line[2];
		$line[2]++
	}
	#print join("\t",@line);
		#print "\n";
	print "$line[0]\t.\t.\t$line[1]\t$line[2]\t.\t.\t.\t.\n";
}