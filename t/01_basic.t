use strict;
use warnings;
use File::Basename qw/dirname/;

use Test::More tests => 11;

use_ok 'Bio::Kmer';

my $kmer=Bio::Kmer->new(dirname($0)."/../data/rand.fastq.gz",{kmerlength=>8});
$kmer->count();
my $hist=$kmer->histogram();

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

for(my $i=0;$i<@correctCounts;$i++){
  is $$hist[$i], $correctCounts[$i], "Freq of $i checks out";
}

