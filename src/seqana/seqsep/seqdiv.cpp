#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iterator>
#include <cstdlib>
#include "nameoutfile.h"

//file: seqdiv.cpp

/** Separates GB sequence according to keywords in the lineage line
 * following the ORGANISM record 
 * mainly for separating sau and amp divisions from the vrt division
 *
 * Adding the ability to read from a config file example in seqdiv.conf
 * in this source code directory
 *
 * I will hardcode this into this source code, so that the caller don't
 * have to call the config file
*/

void getInput(char inf[]);
//void nameFile(const char ipf[], string &s, string &nons, const string &name);
void toLower(string &s) { 
	for (unsigned int i=0; i<s.length(); i++)  s[i] = tolower(s[i]);
}
map<string, vector<string> > readConf(const char cnf[]);

/** helper function to split string into vector of words
 */
vector<string> split(string str) {
	vector<string> tmp;
	istringstream IN(str);
	string wd;
	while (IN.eof()) {
		IN >> wd;
		tmp.push_back(wd);
	}
	return tmp;
}

int main(int argc, char *argv[])
{
	char gbfile[51] = "gbvrt.seq";
	vector<string> divkwd;  // keywords are case-sensitive
	vector<string> kwd;  // for reading into map

	/* loading the keyword */
	map<string, vector<string> > divkwdmap;
	kwd.push_back("Archosauria");
	kwd.push_back("Lepidosauria");
	kwd.push_back("Testudines");
	divkwdmap["sau"] = kwd;
	kwd.clear();
	kwd.push_back("Eutheria");
	divkwdmap["plc"]=kwd;
	kwd.clear();
	kwd.push_back("Amphibia");
	divkwdmap["amp"]=kwd;
	kwd.clear();
	kwd.push_back("Actinopterygii");
	kwd.push_back("Coelacanthiformes");
	kwd.push_back("Dipnoi");
	divkwdmap["bony"] = kwd;
	kwd.clear();
	kwd.push_back("Chondrichthyes");
	divkwdmap["cart"] = kwd;
	/////////////////////////////////
	
	string divName("sau");  // must be specified
	string configfile;

	string usage("seqdiv [-s integer] -f inputfile -d division_name keyword1 keyword2 ...\n");
	int i, show=5000;

	switch (argc) {
	case 1:
		cerr << usage;
		exit(1);
		break;
	default:
		i = 1;
		while (i<argc) {
			if (!strcmp(argv[i], "-f")) strcpy(gbfile,argv[++i]);
			else if (!strcmp(argv[i], "-s")) show = atoi(argv[++i]);
			else if (!strcmp(argv[i], "-d")) divName = argv[++i];
			else if (!strcmp(argv[i], "-r")) configfile = argv[++i];
			else divkwd.push_back(string(argv[i]));
			i++;
		}
		break;
	}
	if (!configfile.empty()) readConf(configfile.c_str());

	/* You needs only to specify the division name, the program will 
	 * automatically load the keywords
	 */
	if (divkwd.empty()) {
		if (divkwdmap.find(divName) != divkwdmap.end()) 
			divkwd = divkwdmap[divName];
		else {
			cerr << "division keyword not specified!\n";
			exit(1);
		}
	}
	cerr << "division: " << divName << " keywords: ";
	copy(divkwd.begin(), divkwd.end(), ostream_iterator<string>(cerr, " "));
	cerr << endl;

	string thisSeq, otherSeq;
	nameFile(gbfile, thisSeq, otherSeq, divName);

	ifstream INF(gbfile);  //input file containing the GB files
	if (INF.fail()) {
		cout << "input file: " << gbfile << " not found\n";
		exit(1);
	}
	ofstream THISORG(thisSeq.c_str());    // files with seq from thisdiv
	ofstream REST(otherSeq.c_str());  //all author objects

	int totalSeq = 0, totalF =0, totalNf=0;

	string ln, buff, acc;
	getline(INF, ln);
	while (ln.substr(0, 5) != "LOCUS") {
		cout << ln << endl; // to show the header information
		getline(INF, ln);
	}
	bool found=false;

	while (ln.substr(0, 5) == "LOCUS" && !INF.eof()) {
		totalSeq++;
		//if (totalSeq > 18000) break;  // for debug
		if (totalSeq%show == 0) {
			cout << '#' << totalSeq << ":\n";
			cout << ln << endl;
		}
		ostringstream sbff; //out buffer
		sbff << ln << endl;  //ln contains the LOCUS line

		getline(INF, ln);
		while (ln.substr(2, 8) != "ORGANISM") {
			sbff << ln << endl;
			getline(INF, ln);
		}
		sbff << ln << endl;  // Organism line into buffer

		string lineage;
		getline(INF,lineage);
		getline(INF,ln);
		while (ln[0] == ' ') {
			lineage.append("\n").append(ln);
			getline(INF, ln);
		}
		found = false;
		for (unsigned int j=0; j<divkwd.size(); j++) {
			if (lineage.find(divkwd[j]) != string::npos) {
				found = true;
				break;
			}
		}
		if (found) {
			totalF++;
			THISORG << sbff.str() << lineage << endl;
			while (ln.substr(0, 5) != "LOCUS" && !INF.eof()) {
				THISORG << ln << endl;
				getline(INF, ln);
			}
		}
		else {
			totalNf++;
			REST << sbff.str() << lineage << endl;
			while (ln.substr(0, 5) != "LOCUS" && !INF.eof()) {
				REST << ln << endl;
				getline(INF, ln);
			}
		}
	}

	cout << "In genbank file: " << gbfile << "\n";
	cout << "Total " << totalSeq << " sequences\n";
	cout << "Total " << totalF << " " << divName << " sequences\n";
	cout << "Total " << totalNf << " Non-" << divName << " sequences\n";
	
	INF.close();
	THISORG.close();
	REST.close();

	return 0;
}

void getInput(char inf[])
{
	cout << "Enter input file in GB format\n";
	cout << "no more than 50 characters\n";
	cin.getline(inf, 50);
}

/*
void nameFile(const char ipf[], string &s, string &nons, const string &name) {
	char *p = strchr(ipf, '.');
	if (p) {
		s = string(ipf, p-ipf);
	}
	else s = ipf;
	nons = s;
	s += '.';
	s += name;
	nons += ".other";
}
*/
map<string, vector<string> > readConf(const char cnf[]) {
	map<string, vector<string> > divkwd;
	ifstream IN(cnf);
	if (IN.fail()) {
		cerr << "Cannot read configure file: " << cnf << endl;
		exit(1);
	}
	string div, kwdstr;
	vector<string> kwdvec;

	getline(IN, div);
	while (!IN.eof()) {
		while (!IN.eof() && div[0] == '#') getline(IN, div);
		if (IN.eof()) break;
		getline(IN, kwdstr);
		kwdvec = split(kwdstr);
		divkwd[div] = kwdvec;
		getline(IN, div);
	}
	return divkwd;
}

