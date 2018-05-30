#include <iostream>
#include <string>
#include <sstream>
#include <map>

using namespace std;
#include "InversionAlign.h"

#include "FASTAReader.h"
#include "FASTASequence.h"

ifstream samIn;
vector<FASTASequence> genome;
map<string, FASTASequence*> chromMap;
MappingSemaphores semaphores;
ofstream tableOut;
int seqIndex;


int main(int argc, char* argv[]) {
	string queryFileName = argv[1];
	string genomeFileName = argv[2];
	string outputTable = argv[3];
	int argi = 4;
	int k = 15;
	while (argi < argc) {
		if (strcmp(argv[argi], "-k") == 0) {
			k = atoi(argv[++argi]);
		}
		++argi;
	}
	tableOut.open(outputTable.c_str());
	//
	// Since many alignments will be performed against the genome, 
	// read it all at once, rather than seeking.
	//

	FASTAReader reader;
	reader.Initialize(genomeFileName);
	FASTASequence target;
	genomeReader.GetNext(target);
	reader.Initialize(queryFileName);
	FASTASequene query;


	vector<StrandedFragment> fragmentSet;
	StoreMatches(query, target, k, fragmentSet);

	int f;
	for (f = 0; f < fragmentSet.size(); f++) {
		

	}
}

