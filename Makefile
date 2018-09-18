SAMTOOLS=htslib
BLASR=blasr/common

ST_OPS=-I $(SAMTOOLS)  -L $(SAMTOOLS)  -lhts -lm -lz
all: sdplite screenInversions 

$(SAMTOOLS)/htslib/libhts.a:
	cd $(SAMTOOLS) && autoheader && autoconf && ./configure --prefix=htslib && make

screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -std=c++98 -O3 $< -o $@ -I$(SAMTOOLS)/htslib -I $(BLASR) -lpthread

sdplite: TestSDPAlignLite.cpp InversionAlign.h
	g++ -g $< -o $@ -I $(BLASR)
