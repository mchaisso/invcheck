//#include "htslib/bam.h"
#include "htslib/sam.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>

//
// Blasr source includes.
//
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment.h"

int main(int argc, char* argv[]) {
	
	if (argc < 4) {
		cout << "usage: realignReverse out.bed reference.fasta [-in <filename> -in  <filename> ... ]  ... [options ]" << endl;
		cout << " Prints coverage by base." << endl;
		cout << " Options: " << endl;
		cout << "   -q q   min mapping quality (30)." << endl;
		exit(1);
	}
	string outFileName;
	outFileName = argv[1];

	int argi = 2;
	int minQuality = 30;
  int bin = 50;
	vector<string> bamFileNames;
	bool skipZero = false;
	while ( argi < argc ) {
		if (strcmp(argv[argi], "-q") == 0) {
			++argi;
			minQuality = atoi(argv[argi]);
			++argi;
		}
		else if (strcmp(argv[argi], "-in") == 0) {
			++argi;
			bamFileNames.push_back(argv[argi]);
			++argi;
		}
	}
	

	ofstream outFile(outFileName.c_str());
	int i;
	int bamI = 0;
	vector<vector<int> > coverage;
	bam_header_t *header;
	int readIndex =0 ;
	long totalNumBases = 0;

	for (bamI = 0; bamI < bamFileNames.size(); bamI++) {
		cerr << bamFileNames[bamI] << endl;
		BGZF *in;
		in = bam_open(bamFileNames[bamI].c_str(), "rb");
		

		header = sam_header_read(in);
		bam1_t *b =  bam_init1();
	
		while (sam_read1(in, b) >= 0) {
		
			int tStart = b->core.pos;
			int tEnd   = b->core.pos + GetTLen(b);
			if (b->core.qual >= minQuality) {
				vector<int>* v = &coverage[b->core.tid];
				for (i = tStart; i < tEnd; i++) {
					(*v)[i/bin]+=1;
				}
			}
			++readIndex;
			totalNumBases += tEnd - tStart;
			if (readIndex % 1000000 == 0) {
				cerr << "processed " << readIndex << " reads " << totalNumBases / 1000000000 << "Gb" << endl;
			}

			bam_destroy1(b);
			b = bam_init1();
		}
		bam_close(in);
	}

