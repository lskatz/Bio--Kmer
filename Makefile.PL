use strict;
use warnings;
 
use 5.008;
 
use ExtUtils::MakeMaker;
WriteMakefile
(
  NAME         => 'Bio::Kmer',
  VERSION_FROM => 'lib/perl5/Bio/Kmer.pm',
  PREREQ_PM    => {
    'File::Basename' => 0,
    'File::Temp'     => 0,
    'Data::Dumper'   => 0,
    'IO::Uncompress::Gunzip' => 0,
    'threads'        => 0,
    'threads::shared'=> 0,
    'Thread::Queue'  => 0,
    'Exporter'       => 0,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
);

