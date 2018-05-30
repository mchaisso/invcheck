SAMTOOLS=htslib
SAMTOOLS=$(HOME)/software/samtools-1.2
BLASR=blasr/common

ST_OPS=-I $(SAMTOOLS)  -L $(SAMTOOLS)  -lhts -lm -lz
all: testInversionAlign sdplite screenInversions screenSplitReads

htslib/libhts.a:
	cd htslib && make

realignReverse: RealignReverse.cpp htslib/libhts.a
	g++ RealignReverse.cpp -I $(BLASR)  $(ST_OPS) -o $@ 


testInversionAlign: TestInversionAlign.cpp InversionAlign.h
	g++ -O3 TestInversionAlign.cpp -o $@ -I $(BLASR)


screenInversions: ScreenInversions.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread


screenSplitReads: ScreenSplitReads.cpp InversionAlign.h
	g++ -O3 $< -o $@ -I $(BLASR) -lpthread -I $(SAMTOOLS) -I $(SAMTOOLS)/htslib-1.2.1  -L $(SAMTOOLS) -lbam  -lz -L$(SAMTOOLS)/htslib-1.2.1 -l hts

sdplite: TestSDPAlignLite.cpp InversionAlign.h
	g++ -g $< -o $@ -I $(BLASR)
