swig -perl beagle.i
# `perl -MExtUtils::Embed -e ccopts` doesn't work on OSX, adds -arch ppc
gcc -fPIC -c beagle_wrap.c -I. -I/System/Library/Perl/5.10.0/darwin-thread-multi-2level/CORE/ `pkg-config --cflags hmsbeagle-1`
# `perl -MExtUtils::Embed -e ldopts` doesn't work on OSX
gcc -bundle -o beagle.bundle beagle_wrap.o -arch x86_64 -L/usr/local/lib  -L/System/Library/Perl/5.10.0/darwin-thread-multi-2level/CORE -lperl -ldl -lm -lutil -lc `pkg-config --libs hmsbeagle-1`
