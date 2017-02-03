use strict;
use warnings;
 
use 5.008;
 
use ExtUtils::MakeMaker;
WriteMakefile
(
  NAME         => 'Bio::Kmer',
  VERSION_FROM => 'lib/Bio/Kmer.pm',
  PREREQ_PM    => {
    'File::Basename' => 0.1,
    'File::Temp'     => 0.1,
    'Data::Dumper'   => 0.1,
    'IO::Uncompress::Gunzip' => 0.1,
    'threads'        => 0.1,
    'threads::shared'=> 0.1,
    'Thread::Queue'  => 0.1,
    'Exporter'       => 0.1,
  },
  BUILD_REQUIRES => {
    'Test::More'     => 0.47,
  },
);

