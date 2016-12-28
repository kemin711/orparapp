#include <fstream>
#include <cstring>
#include <cstdlib>

#include "match.h"

using namespace std;

/* file: selortho.cpp
 * input: blastp -m 8 output
 * This is a blast post processing program, it creates new 
 * filtered results for loading into a relational table.
 * output: unique query x target pair consolidated
 *
 * top hits not selected yet. 
 * */

void usage();
string nameOutFile(const string &inf, string ext);

int main(int argc, char* argv[]) {
	if (argc < 2) {
		usage();
	}
	/* default values */
	float scoreCut=0.9, identityCut=0.95, lengthCut=0.7;
	// score>scoreCut or (identity>identityCut && length>lengthCut)
	string infile = "test.blp";
	int i=1;
	while (i < argc) {
		if (!strcmp(argv[i], "-sc")) scoreCut = atof(argv[++i]);
		else if (!strcmp(argv[i], "-lc")) lengthCut = atof(argv[++i]);
		else infile = argv[i];
		i++;
	}
	ifstream IN(infile.c_str());
	if (IN.fail()) {
		cerr << infile << " open failed\n";
		exit(1);
	}
	string rawfile = nameOutFile(infile, "raw");
	ofstream OU(rawfile.c_str());
	string selfile = nameOutFile(infile, "sel");  // selected results
	ofstream SEL(selfile.c_str());
	cerr << "output fields:\n" << hit::fields << endl;
	//cerr << "query\ttarget\taverage_identity\tsum_match_length\tsum_score\tselratio\n";
	string nextQuery, nextTarget;

	IN >> nextQuery >> nextTarget;
	int count=0;
	try {
		while (!IN.eof()) {    // loop through each query
			count++;
			if (count%200 == 0) { cerr << '.'; }
			hit topHit(nextQuery, nextTarget, IN);
			float maxIden = topHit.getIdentity();
			int maxLen = topHit.getLength();
			topHit.dumpRaw(OU);
			SEL << topHit << endl;
			while (topHit.sameQuery(nextQuery)) {
				hit ahit(nextQuery, nextTarget, IN);
				//cerr << ahit << endl;
				if (ahit.getIdentity() > maxIden) 
					maxIden = ahit.getIdentity();
				if (ahit.getLength() > maxLen) 
					maxLen = ahit.getLength();

				if (ahit.scoreAsBigAs(topHit, scoreCut)) {
					SEL << ahit << endl;
					ahit.dumpRaw(OU);
				}
				else if (ahit.getIdentity() > identityCut*maxIden && 
						ahit.getLength() > lengthCut*maxLen ) {
					//cout << "From next part ****\n";   // debug
					SEL << ahit << endl;
					ahit.dumpRaw(OU);
				}
			}
		}
	}
	catch(inputException) {
		cerr << "Input file has some problem at " << count <<endl;
		return 1;
	}
	catch(hit::end) {
		cerr << count << " queries done\n";
		cerr << "Selected result written to " << selfile << "\n";
		return 0;
	}

	return 0;
}

void usage() {
	cerr << "selortho -sc scoreCut (default=0.9) \n\tinfile\n";
	cerr << "\t-lc lengthCut (default=0.7)\n";
	exit(1);
}

string nameOutFile(const string &inf, string ext) {
	int i = inf.rfind('.');
	return (inf.substr(0, i+1) + ext);
}
