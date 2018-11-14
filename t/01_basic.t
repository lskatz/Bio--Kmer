use strict;
use warnings;
use File::Basename qw/dirname/;
use File::Temp qw/tempdir/;
use FindBin qw/$RealBin/;
use IO::Uncompress::Gunzip qw/gunzip $GunzipError/;
use Data::Dumper qw/Dumper/;

use Test::More tests => 1;

use lib "$RealBin/../lib";
use_ok 'Bio::Kmer';

diag "Made it to line ".__LINE__;

# expected histogram
my @correctCounts=(
  0,
  16087,
  17621,
  12868,
  6857,
  3070,
  1096,
  380,
  105,
  17,
# 6,
);
diag "Made it to line ".__LINE__;

# expected query results
my %query=(
  TTGGAGCA => 3,
  TTGGAGCT => 6,
  TTGGAGCTA=> -1, # invalid
  AAAAAAAA => 0,  # not found
);
diag "Made it to line ".__LINE__;

# Test pure perl
my $kmer=Bio::Kmer->new(dirname($0)."/../data/rand.fastq.gz",{kmerlength=>8});
diag "Made it to line ".__LINE__;
my $hist=$kmer->histogram();
diag "Made it to line ".__LINE__;
for(my $i=0;$i<@correctCounts;$i++){
  diag "d Frequency: ".$$hist[$i]." <=> $correctCounts[$i]";
  note "n Frequency: ".$$hist[$i]." <=> $correctCounts[$i]";
  print"p Frequency: ".$$hist[$i]." <=> $correctCounts[$i]\n";
  print STDERR "p Frequency: ".$$hist[$i]." <=> $correctCounts[$i]\n";
  is $$hist[$i], $correctCounts[$i], "Freq of $i checks out";
}
for my $query(keys(%query)){
  is $query{$query}, $kmer->query($query), "Queried for $query{$query}";
}
$kmer->close();

# Test subsampling: a subsample should have fewer kmers than
# the full set but more than 0.
my $subsampleKmer=Bio::Kmer->new(dirname($0)."/../data/rand.fastq.gz",{kmerlength=>8,sample=>0.1});
my $subsampleHist=$kmer->histogram();
my $subsampleKmerHash=$subsampleKmer->kmers();
my $numSubsampledKmers = scalar(keys(%$subsampleKmerHash));
my $numKmers = scalar(keys(%{ $kmer->kmers() }));

ok(($numSubsampledKmers > 0), "Subsample kmers, and there are a nonzero count of results.");

ok(($numSubsampledKmers < $numKmers), "Subsample kmers fewer than full count of kmers");


