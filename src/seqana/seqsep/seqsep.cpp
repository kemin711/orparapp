#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <map>
#include <sstream>
#include <cstdlib>

using namespace std;

/*this program will take a gbformated sequence file
sequences from one organism such as Fugu rubripes or human will be 
separated from other sequences 
it will also output a list of gi, list of Acc files
remember there are two names for Fugu Fugu rubripes | Takigugu rubripes
I have used alias to solve this problem.  One extraction is needed
for Fugu sequences

Copy Right of Kemin Zhou
*/

void usage(const map<string, string> &m);

void getInput(char inf[]);
// name output files to store sequence of this and orther orgnanism
void nameFile(const char ipf[], string &s, string &nons, const string &name);

// converts a string to its lower case version
void toLower(string &s) { 
	for (size_t i=0; i<s.length(); i++)  s[i] = tolower(s[i]);
}

int main(int argc, char *argv[])
{
	char gbfile[51];
	string orgCommonName("fugu");  //default
	string org("Fugu rubripes"); //default
	//string mito("Mitochondrion " + org);
	int show = 1000;
	bool OUTGI = false;  //output gi and version information
	bool QUIET = false;  // No information to the terminal

	map<string, string> orgNameMap;        // for lazy typing by the user
	orgNameMap["human"] = "Homo sapiens";
	//orgNameMap["Human"] = "Homo sapiens";
	orgNameMap["fugu"] = "Fugu rubripes";
	//orgNameMap["Fugu"] = "Fugu rubripes";
	orgNameMap["mouse"] = "Mus musculus";
	orgNameMap["rat"] = "Rattus norvegicus";

	switch (argc) {
	case 1:
		usage(orgNameMap);
		exit(1);
		break;
	case 2:
		cout << "Will use default organims Fugu rubripes\n";
		strcpy(gbfile, argv[1]);  // use default organism Fugu rubripes
		break;
	default:
		int i = 1;
		while (i<argc) {
			if (!strcmp(argv[i], "-org")) {
				i++;
				orgCommonName = argv[i];
				toLower(orgCommonName);
				if (orgNameMap.find(argv[i]) != orgNameMap.end()) {
					org = orgNameMap[argv[i]];
				}
				else org = argv[i];  // supplied is a scientific name
			}
			else if (!strcmp(argv[i], "-s")) {
				show = atoi(argv[++i]);
			}
			else if (!strcmp(argv[i], "-og")) { OUTGI = true; }
			else if (!strcmp(argv[i], "-Q")) { QUIET = true; }
			else strcpy(gbfile, argv[i]);  // input file name
			i++;
		}
		break;
	}

	string orgalias; //, mitoalias("Mitochondrion ");

	// some organisms have more than one nmae, these are the ones
	if (org == "Fugu rubripes") {
		orgalias = "Takifugu rubripes"; //more than one name
	}
	else if (org == "Mus musculus") {
		orgalias = "Mus sp.";
	}
	size_t i=0;
	//remove space in common name for file extension
	while ( (i=orgCommonName.find(' ', i)) != string::npos) {
		orgCommonName[i] = '_';
		i++;
	}

	string thisSeq, otherSeq;
	nameFile(gbfile, thisSeq, otherSeq, orgCommonName);

	ifstream INF(gbfile);  //input file containing the GB files
	if (INF.fail()) {
		cout << "input file not found\n";
		exit(1);
	}
	ofstream THISORG(thisSeq.c_str());  
	//outfile containing sequences of this organism
	ofstream REST(otherSeq.c_str());  //all author objects
	ofstream GI, ACC;

	if (OUTGI) {
		GI.open("thisOrg.gi");
		ACC.open("otherOrg.acc");
	}

	int totalSeq = 0, thisOrgCnt =0, restOrgCnt=0;

	string ln, buff, acc, gi;
	//int gi;
	getline(INF, ln);
	while (ln.substr(0, 5) != "LOCUS") {
		cout << ln << endl; // to show the header information
		getline(INF, ln);
	}

	while (!INF.eof() && ln.substr(0, 5) == "LOCUS") {
		totalSeq++;
		if (totalSeq%show == 0 && !QUIET) {
			cerr << '#' << totalSeq << "\n";
			cerr << ln << endl;
		}

		ostringstream sbff; //out buffer
		sbff << ln << endl;  //ln contains the LOCUS line

		if (OUTGI) {
			getline(INF, ln);    //ln contains the first line of DEFINITION
			sbff << ln << endl;
			getline(INF, ln);
			while (ln.substr(0, 9) != "ACCESSION") {
				sbff << ln << endl;
				getline(INF, ln);
			}  //reads the rest of DEFINITION line if more than one lines
			sbff << ln << endl;
			acc=ln.substr(12);
			
			getline(INF, ln);
			while (ln.substr(0, 7) != "VERSION") {
				sbff << ln << endl;
				getline(INF, ln);
			}  //reads the second or more of the ACCESSION lines
			//very few records have ACESSION line longer than one line
			sbff << ln << endl;
			gi = ln.substr(ln.find("GI:", 12)+3);
		//	gi+=3;
		//	gi = atoi(ln.substr(gi).c_str());
		}

		getline(INF, ln);
		while (ln.substr(2, 8) != "ORGANISM") {
			sbff << ln << endl;
			getline(INF, ln);
		}
		unsigned int altIndex;  // alternative name index in the line
		if (ln.find(org) != string::npos) { // this organism
			thisOrgCnt++;
			THISORG << sbff.str();
			if (OUTGI) {
				GI << gi << endl;
				ACC << acc << endl;
			}
			while (!INF.eof() &&  ln.substr(0, 5) != "LOCUS") {
				THISORG << ln << endl;
				getline(INF, ln);
			}
		}
		else if (!orgalias.empty() && 
				(altIndex = ln.find(orgalias)) != string::npos) {
			/* replace the alias name with the normal name */
			ln.replace(altIndex, orgalias.length(), org);
			thisOrgCnt++;
			THISORG << sbff.str();
			if (OUTGI) {
				GI << gi << endl;
				ACC << acc << endl;
			}
			while (!INF.eof() && ln.substr(0, 5) != "LOCUS" ) {
				THISORG << ln << endl;
				getline(INF, ln);
			}
		}
		else { //sequence of other organism
			restOrgCnt++;
			REST << sbff.str();
			while (ln.substr(0, 5) != "LOCUS" && !INF.eof()) {
				REST << ln << endl;
				getline(INF, ln);
			}
		}
	}

	cout << "In genbank file: " << gbfile << "\n";
	cout << "Total " << totalSeq << " sequences\n";
	cout << "Total " << thisOrgCnt << " " << org << " sequences\n";
	cout << "Total " << restOrgCnt << " Non-" << org << " sequences\n";
	
	INF.close();
	THISORG.close();
	REST.close();
	if (OUTGI) {
		GI.close();
		ACC.close();
	}
	
	return 0;
}

void usage(const map<string, string> &m) {
	map<string, string>::const_iterator i = m.begin();
	string usageString("Seprates sequences of one organism from those of ");
	usageString.append("others\n");
	usageString.append("usage: seqsep infile -org [orgnaims name]\n");
	usageString.append("-s [num]: show Locus Line for every num sequence\n");
	usageString.append("-og for output all gi into one file\n");
	usageString.append("-Q for quiet, no screen output\n");

	cout << usageString << endl;
	cout << "You can use:\n";
	while (i != m.end()) {
		cout << i->first << " for " << i->second << endl;
		i++;
	}
}

void getInput(char inf[])
{
	cout << "Enter input file in GB format\n";
	cout << "no more than 50 characters\n";
	cin.getline(inf, 50);
}

/*
void nameOutfile(const char ipf[], char fs[], char nonfs[])
{
	char gbf[50]; 
	int i = 0;
	while (ipf[i] != '.' && ipf[i] != '\0') i++;
	strncpy(gbf, ipf, i); //discard .xxx extension
	gbf[i]='\0';
	strcpy(fs, gbf);
	strcat(fs, ".thisOrg");
	strcpy(nonfs, gbf);
	strcat(nonfs, ".otherOrg");
}
*/

void nameFile(const char ipf[], string &s, string &nons, const string &name) {
	const char *p = strchr(ipf, '.');
	if (p) {
		s = string(ipf, p-ipf);
	}
	else s = ipf;
	nons = s;
	s += '.';
	s += name;
	nons += ".other";
}

