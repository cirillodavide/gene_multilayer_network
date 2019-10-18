use warnings;
use strict;
use Algorithm::Combinatorics qw(combinations);
use Data::Dumper;
use POSIX qw(strftime);

my $date = strftime "%d-%m-%Y", localtime;

# download Reactome lowest levels
my $url = 'https://reactome.org/download/current/NCBI2Reactome.txt';
system("wget $url -P src");

my $file = 'src/NCBI2Reactome.txt';
my %reactome;
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	my($gene,$pathway,$ec,$taxo) = (split/\t/,$_)[0,1,4,5];
	next unless $taxo eq "Homo sapiens";
	next unless $gene=~/^\d+$/;
	#next if $ec !~ m/EXP|IDA|IPI|IMP|IGI|IEP/; # they are all IEA and TAS...
	push @{$reactome{$pathway}}, $gene unless grep{$gene eq $_}@{$reactome{$pathway}};
}
close FILE;

open OUT, '>', 'networks/Reactome_pathways.'.$date.'.gr' or die $!;
my %seen;
foreach my $pathway(keys %reactome){
	next if scalar @{$reactome{$pathway}} < 2;
	my $genes = [@{$reactome{$pathway}}];
	my $iter = combinations($genes, 2);
	while (my $c = $iter->next) {
		my $g1 = @$c[0];
		my $g2 = @$c[1];
		my @arr = ($g1,$g2);
		my @s_arr = sort{$a<=>$b}@arr;
		push @s_arr, $pathway;
		$seen{"@s_arr"} = 1;
	}
}
foreach my $tag(keys %seen){
	print OUT $tag,"\n";
}
close OUT;
