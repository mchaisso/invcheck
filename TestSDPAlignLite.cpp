#include "InversionAlign.h"
#include "FASTAReader.h"
#include "FASTASequence.h"


int main(int argc, char* argv[]) {

	string queryName, targetName;
	if (argc < 3) {
		cout << "usage: testsdpl query target k" << endl;
		exit(1);
	}
	queryName = argv[1];
	targetName = argv[2];
	int k = atoi(argv[3]);

	FASTAReader reader;
	reader.Initialize(queryName);
	FASTASequence query, target;
	reader.GetNext(query);
	reader.Initialize(targetName);
	reader.GetNext(target);
	vector<int> coords;
	InvParameters params;
	cout << SDPAlignLite(query, target, params) << endl;;
}
