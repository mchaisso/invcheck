HTSLIB=htslib
SAMTOOLS=samtools
BLASR=blasr/common

ST_OPS=-I $(SAMTOOLS)  -L $(SAMTOOLS)  -lhts -lm -lz
all:  screenInversions

htslib/libhts.a:
	cd htslib && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 && make -j 8 

samtools/samtools:
	cd samtools && autoheader && autoconf && ./configure --disable-lzma --disable-bz2 && make -j 8 

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread
