use ExtUtils::MakeMaker;

WriteMakefile(
    'AUTHOR'   => 'Rutger Vos <rutgeraldo@gmail.com>',
    'NAME'     => 'Bio::Phylo::Beagle',
    'ABSTRACT' => 'wrapper around Beagle-lib',
    'VERSION'  => '0.01',
    'OBJECT'   => 'beagle_wrap.o',
    # Add extra include paths here, if needed
	'INC'    => `pkg-config --cflags hmsbeagle-1`,
	# Add shared and static libraries here, if needed
	'LDFROM' => 'beagle_wrap.o ' . `pkg-config --libs hmsbeagle-1`,
);