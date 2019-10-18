use warnings;
use strict;
use Data::Dumper;
use POSIX qw(strftime);

my $date = strftime "%d-%m-%Y", localtime;

# check current release here: https://downloads.thebiogrid.org/BioGRID

my $release = '3.5.177';
my $url = 'https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-'.$release.'/BIOGRID-ALL-'.$release.'.tab2.zip';
my $zip = 'src/BIOGRID-ALL-'.$release.'.tab2.zip';
my $file = $zip;
$file =~ s/zip$/txt/g;
system("curl -o $zip $url");
system("unzip $zip -d src");
system("rm $zip");

# create the graph

my %biogrid;
open FILE, $file or die $!;
while(<FILE>){
	chomp;
	next if $. < 2;
	my($int1,$int2,$expType,$taxo1,$taxo2) = (split/\t/,$_)[1,2,12,15,16];
	next unless $taxo1 eq "9606" and $taxo2 eq "9606";
	push @{$biogrid{$int1}{$expType}}, $int2 unless grep{$int1 eq $_}@{$biogrid{$int2}{$expType}} or grep{$int2 eq $_}@{$biogrid{$int1}{$expType}};
}
close FILE;

open OUT, '>', 'networks/BioGRID_interactions.'.$date.'.gr' or die $!;
foreach my $k(keys %biogrid){
	foreach my $e(keys %{$biogrid{$k}}){
		foreach my $v(@{$biogrid{$k}{$e}}){
			print OUT join(" ",$k,$v,$e),"\n";
		}
	}
}
close OUT;