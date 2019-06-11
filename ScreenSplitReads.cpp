#include "samtools/sam.h"
#include "stdlib.h"
#include <vector>
#include <string>
#include <set>
#include <iostream>
#include <fstream>
#include <map>
#include <sstream>
#include <numeric>

#include "FASTAReader.h"
#include "FASTASequence.h"
#include "InversionAlign.h"

using namespace std;
typedef struct {  
    int beg, end;  
    samfile_t *in;  
} tmpstruct_t;  
  

class Output {
public:
	ofstream *outFilePtr;
	int minClipping;
	int minq;
	string chrom;
};


int GetTLen(const bam1_t *b) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int tlen = 0;
	int i;
	for (i = 0; i < len; i++) {
		int op = cigar[i] & 0xf;
		int oplen = cigar[i] >> 4;
		if (op == BAM_CMATCH or op == BAM_CDEL) {
			tlen += oplen;
		}
	}
	return tlen;
}

int GetStrand(const bam1_t *b) {
	if (b->core.flag & 16) {
		return 1;
	}
	else {
		return 0;
	}
}


int GetClipping(const bam1_t *b, int &leftClip, int &rightClip) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int qlen = 0;
	int i;
	leftClip = rightClip = 0;
	
	if (len > 0) {
		int first = 0;
		int op = cigar[first] & 0xf;
		int oplen = cigar[first] >> 4;
		if (op == BAM_CHARD_CLIP) {
			first++;
			op = cigar[first] & 0xf;
			oplen = cigar[first] >> 4;
		}
		if (op == BAM_CSOFT_CLIP) {
			leftClip = oplen;
		}
		int last = len - 1;
		if (last > first) {
			op = cigar[last] & 0xf;
			oplen = cigar[last] >> 4;
			if ( op == BAM_CHARD_CLIP) {
				last--;
			}
			if (last > first) {
				op = cigar[last] & 0xf;
				oplen = cigar[last] >> 4;
				if (op == BAM_CSOFT_CLIP) {
					rightClip = oplen;
				}
			}
		}
	}
	return 0;
}

int numProcessed = 0;

const int MIN_ALIGNED_LENGTH = 500;
int numFound = 0;
int WriteHardStop(const bam1_t *b, void *data) {
  Output *output = (Output*) data;
	// 
	// determine overlap with exon
	// 
	if (b->core.qual < output->minq) {
		return 0;
	}
	int leftClipping, rightClipping;
	GetClipping(b, leftClipping, rightClipping);
	int tStart = b->core.pos;
	int tEnd   = b->core.pos + GetTLen(b);
	if (b->core.l_qseq - leftClipping - rightClipping  < MIN_ALIGNED_LENGTH ) {
		return 0;
	}
	if (leftClipping > output->minClipping or rightClipping > output->minClipping) {
		(*output->outFilePtr) << output->chrom << "\t" << tStart << "\t" << tEnd << "\t" << bam1_qname(b) << "\t" << leftClipping << "\t" << rightClipping << "\t" << GetStrand(b) << endl;
	}		
	++numProcessed;
}


int GetQLen(const bam1_t *b) {
	uint32_t *cigar = bam1_cigar(b);
	int len = b->core.n_cigar;
	int qlen = 0;
	int i;
	for (i = 0; i < len; i++) {
		int op = cigar[i] & 0xf;
		int oplen = cigar[i] >> 4;
		if (op == BAM_CMATCH or op == BAM_CINS) {
			qlen += oplen;
		}
	}
	return qlen;
}

class Entry {
public:
	string seq;
	int refStart, refEnd, qStart, qEnd;
	int mapqv;
	int strand;
	string readName;
	int startClip, endClip;
	Entry(string s, int rs, int re, int qs, int qe, int m, int st): 
		seq(s), refStart(rs), refEnd(re), qStart(qs), qEnd(qe), mapqv(m), strand(st) {};
	Entry() {
		seq = "";
		refStart= refEnd = qStart = qEnd = mapqv = 0;
	}
	Entry(bam1_t *b, int startClipP, int endClipP) {
		int qlen = bam_cigar2qlen(&b->core, bam1_cigar(b));
		seq.resize(qlen);
		int i;
		uint8_t *bitSeq = bam1_seq(b);
		for (i = 0; i < seq.size(); i++) { seq[i] = bam_nt16_rev_table[bam1_seqi(bitSeq, i)];}
		refStart = b->core.pos;
		refEnd   = refStart + GetTLen(b);
		qStart = startClipP;
		qEnd   = seq.size() - endClipP;
		mapqv = b->core.qual;
		strand = GetStrand(b);
		readName = bam1_qname(b);
		startClip = startClipP;
		endClip = endClipP;
	}
		
};

InvParameters params;

int main(int argc, char* argv[]) {
	

	bam_index_t *idx;
	bam_plbuf_t *buf;

	bamFile bamFile;

	if (argc < 4) {
		cout << "usage:  input.bam genome minClipping maxSpan out.bed " << endl;
		exit(1);
	}

	samfile_t *in;  
	int argi = 1;
	in = samopen(argv[argi++], "rb", 0);

	string genomeFileName = argv[argi++];
	int minClipping = atoi(argv[argi++]);
	int maxSpan    = atoi(argv[argi++]);
	string outFileName = argv[argi++];
	int minQV = 20;
	while (argi < argc) {
		if (strcmp(argv[argi], "-q") == 0) {
			minQV = atoi(argv[++argi]);
		}
		++argi;
	}

	int printed = 0;

	vector<FASTASequence> genome;
	map<string, FASTASequence*> chromMap;
	FASTAReader genomeReader;

	genomeReader.Initialize(genomeFileName);
	genomeReader.ReadAllSequences(genome);
	cerr << "done reading genome." << endl;
	int i;
	for (i = 0; i < genome.size(); i++) {
		chromMap[genome[i].title] = &genome[i];
	}

  idx = bam_index_load(argv[1]);

	map<string, Entry> alignments;

	ofstream outFile(outFileName.c_str());

	bam1_t *entry = new bam1_t;

	int ref;
	string prevChrom = "";

	bam1_t *b = bam_init1();	
	int readIndex = 0;
	while (bam_read1(in->x.bam, b) > 0) {
		++readIndex;
		if (readIndex % 1000 == 0) {
			cerr << "processed " << readIndex << endl;
		}
		if (b->core.tid < 0) {
			continue;
		}
		
		string chrom = in->header->target_name[b->core.tid];

		if (chrom != prevChrom) {
			cerr << "Resetting alignments!" << endl;
			alignments.clear();
		}
		prevChrom = chrom;
		string qName = bam1_qname(b);
		
		int mapq = b->core.qual;
		map<string, Entry>::iterator queryIt;
		queryIt = alignments.find(qName);

		

		//
		// If this hasn't been seen yet, add the alignmetn.
		//
		if (queryIt == alignments.end()) {
			int leftClipping,rightClipping;
			GetClipping(b, leftClipping, rightClipping);
			if (leftClipping > minClipping or rightClipping > minClipping) {
				alignments[qName] = Entry(b,leftClipping, rightClipping);
			}
		}
		else {

			if (mapq < minQV) {
				alignments.erase(queryIt);
				continue;
			}
			//
			// Otherwise, actually do something with this alignment.
			//
	 
			int leftClipping, rightClipping;
			GetClipping(b, leftClipping, rightClipping);
			if (leftClipping > minClipping or rightClipping > minClipping) {
				//
				// Check the gap between the two sequences
				//
				Entry curEntry(b, leftClipping,rightClipping);
				Entry prevEntry = queryIt->second;
				if (prevEntry.refEnd < curEntry.refStart or 
						prevEntry.refStart > curEntry.refEnd) {
					//
					// Found a potential breakpoint. for now just print.
					//
					if (prevEntry.strand != curEntry.strand) {
						// 
						// Found potential inversion
						//
						
						// Establish the boundaries of the inversion.
						
						int lStart, lEnd, gStart, gEnd;
						int lStrand, gStrand;

						Entry left, right;
						if (prevEntry.refEnd < curEntry.refStart) {
							left = prevEntry;
							right = curEntry;
						}
						else {
							left = curEntry;
							right = prevEntry;
						}
						int bpStart, bpEnd;
						int alnStart, alnEnd;
						int refRegionStart, refRegionEnd;
						int invRefStart, invRefEnd;
						int order = 0; // 0 implies ref then inversion. 1 inversion then ref.

						FASTASequence read;
						read.CopyTitle(prevEntry.readName);
						read.Resize(prevEntry.seq.size());
						memcpy(read.seq, (Nucleotide*)prevEntry.seq.c_str(), prevEntry.seq.size());

						if (left.startClip > minClipping and right.startClip > minClipping and 
								left.endClip <= minClipping and right.endClip <= minClipping) {
							//   ....|---->......|<--#.... 
							// Should swap to 
							bpStart = left.refStart;
							bpEnd = right.refStart;

							invRefStart = right.refStart;
							invRefEnd   = right.refEnd;
							refRegionStart = bpEnd - (left.refEnd - left.refStart);
							refRegionEnd   = right.refEnd;
							order = 1;

							if (left.strand == 0) {
								read.ReverseComplementSelf();
							}
						}
						else if (left.startClip <= minClipping and right.startClip <= minClipping and 
										 left.endClip > minClipping and right.endClip > minClipping) {
							bpStart = left.refEnd;
							bpEnd = right.refEnd;

							invRefStart = left.refStart;
							invRefEnd   = left.refEnd;
							refRegionStart = left.refStart;
							refRegionEnd   = bpStart + (right.refEnd - right.refStart);
						}
						else {
							continue;
						}
									


						//
						// Create the inverted sequence.
						//
						FASTASequence donor;
						FASTASequence ref;
						int inversionLength = bpEnd - bpStart;
						if (inversionLength > 100000) {
							continue;
						}
						ref.ReferenceSubstring(*chromMap[chrom], bpStart, inversionLength);
						ref.MakeRC(donor);
						FASTASequence invertedLocus;
						int invRefLen = invRefEnd - invRefStart;
						invertedLocus.Resize(invRefLen + inversionLength);
						
						if (order == 0) {
							memcpy(invertedLocus.seq, &chromMap[chrom]->seq[invRefStart], invRefLen);
							memcpy(&invertedLocus.seq[invRefLen], donor.seq, donor.length);
						}
						else {
							memcpy(&invertedLocus.seq[0], donor.seq, donor.length);
							memcpy(&invertedLocus.seq[donor.length], &chromMap[chrom]->seq[invRefStart], invRefLen);
						}							
						
						invertedLocus.ToUpper();
							
						FASTASequence refLocus;
						refLocus.ReferenceSubstring(*chromMap[chrom], refRegionStart, refRegionEnd - refRegionStart);

							
						int invScore, refScore;
						invScore = SDPAlignLite(read, invertedLocus, params);
						refScore = SDPAlignLite(read, refLocus, params);


/*
						ref.CopyTitle("inverted_ref");
						donor.CopyTitle("donor");
						ofstream readOut("read.fasta");
						
						read.PrintSeq(readOut);
						refLocus.CopyTitle("ref");
						refLocus.PrintSeq(readOut);
						invertedLocus.CopyTitle("inversion");
						read.PrintSeq(readOut);
						invertedLocus.PrintSeq(readOut);

						read.PrintSeq(readOut);
						donor.PrintSeq(readOut);
						read.PrintSeq(readOut);
						ref.PrintSeq(readOut);

						FASTASequence expRef;
						expRef.ReferenceSubstring(*chromMap[chrom], refRegionStart - 40000, refRegionEnd - refRegionStart + 80000);
						expRef.CopyTitle("expref");
						read.PrintSeq(readOut);
						expRef.PrintSeq(readOut);
						readOut.close();
						delete[] ref.title;
						exit(0);
*/
						cout << chrom << "\t" << bpStart << "\t" << bpEnd << "\t" << refScore << "\t" << invScore << "\t" << invScore - refScore << "\t" << left.readName << endl;
						invertedLocus.Free();
						read.Free();
						donor.Free();
						
					}
					else {
						//
						// Found potential deletion.
						//

					}
				}
			}
			alignments.erase(queryIt);
		}
	}
	outFile.close();
	
}
