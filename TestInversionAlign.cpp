#include "InversionAlign.h"
#include "FASTAReader.h"
#include "FASTASequence.h"


int main(int argc, char* argv[]) {

	string queryName, targetName;
	queryName = argv[1];
	targetName = argv[2];


	FASTAReader reader;
	reader.Initialize(queryName);
	FASTASequence query, target;
	reader.GetNext(query);
	reader.Initialize(targetName);
	reader.GetNext(target);
	int coords[4];
	vector<int> boxes;
	vector<StrandedFragment> fragments;
	InvParameters params;
	InversionAlign(query, target, coords, fragments, boxes, params);

}
