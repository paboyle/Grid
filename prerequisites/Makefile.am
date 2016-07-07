FFTFLAGS=$(filter-out -std=c++11, $(CXXFLAGS) )

EIGENVER=3.2.8
EIGEN=eigen$(EIGENVER)
EIGENTAR=$(EIGEN).tar.bz2
EIGENURL=https://bitbucket.org/eigen/eigen/get/$(EIGENVER).tar.bz2

FFTWVER=3.3.4
FFTW=fftw-$(FFTWVER)
FFTWTAR=fftw-$(FFTWVER).tar.gz
FFTWURL=http://www.fftw.org/$(FFTWTAR)

all: Eigen FFTW headerlist

$(top_srcdir)/prerequisites/$(EIGENTAR): 
	curl -v $(EIGENURL) -o $(top_srcdir)/prerequisites/$(EIGENTAR)

$(top_srcdir)/prerequisites/$(FFTWTAR): 
	curl -v $(FFTWURL) -o $(top_srcdir)/prerequisites/$(FFTWTAR)

Eigen:  $(top_srcdir)/prerequisites/$(EIGENTAR)
	tar xvf $(top_srcdir)/prerequisites/$(EIGENTAR)
	- rm -rf  ../include/Grid/Eigen
	mv eigen-eigen*/Eigen .
	echo EFILES=`find Eigen -type f -name '*.h' ` > $(top_srcdir)/lib/Eigen.inc
	mv Eigen ../include/Grid/
	touch Eigen

FFTW: $(top_srcdir)/prerequisites/$(FFTWTAR)
	tar xvf $(top_srcdir)/prerequisites/$(FFTWTAR)
	cd $(FFTW) &&	./configure --prefix=@abs_top_builddir@/prerequisites/fftwinstall CFLAGS="$(FFTFLAGS)" CC=$(CC) LDFLAGS="$(LDFLAGS)" && make all install
	cp -pr fftwinstall/include/fftw3.h ../include/Grid/
	cp -pr fftwinstall/lib/libfftw3.a  ../lib/
	touch FFTW

headerlist:
	cd $(top_srcdir) && ./scripts/filelist
	touch headerlist