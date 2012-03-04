use ExtUtils::MakeMaker;

WriteMakefile(
    'AUTHOR'   => 'Rutger Vos <rutgeraldo@gmail.com>',
    'NAME'     => 'Bio::Phylo::Beagle',
    'ABSTRACT' => 'wrapper around Beagle-lib',
    'VERSION'  => '0.01',
    #'C'        => [ 'beagle_wrap.c' ],
    #'CCFLAGS'  => `pkg-config --cflags hmsbeagle-1`,
    #'LDFLAGS'  => `pkg-config --libs hmsbeagle-1`,
);