// file: picknewseq.cpp
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
#include <cstring>

#include "version.h"

using namespace std;

/** removes double quotes from beginning and ends*/
void chopEnds(char* &w);

void readAcc(ifstream &ins, set<version> &vset);
void save(ostringstream &bf, ifstream &ins, ofstream &ous);

const int LLEN = 83;  //line length of GB documents

/** Compares the ACC.version file from one file with the GenBank flat
 * file.  Picks new entries from the GenBank flat file
 * First go to ace and run the aql 
 * select s->version from s in class sequence
 * and dump the result as text
 */
int main(int argc, char* argv[]) 
{
	char accfile[50], gbfile[50];

	switch (argc) {
	case 1:
		cerr << "usage picknewseq -k ACC.version -g GenBank_flat_file\n";
		cerr << "version file format: \"acc#.version\"\n";
		return 1;
	default: 
		int i = 1;
		while (i<5) {
			if (!strcmp(argv[i], "-k")) {
				strcpy(accfile, argv[++i]);
			}
			else if (!strcmp(argv[i], "-g")) {
				strcpy(gbfile, argv[++i]);
			}
			else {
				cerr << "something may be wrong in command line\n";
				return 1;
			}
			i++;
		}
		break;
	}

	ifstream ACC(accfile);
	if (ACC.fail()) {
		cerr << "input file " << accfile << " open errro\n";
		exit(1);
	}
	ifstream GB(gbfile);
	if (GB.fail()) {
		cerr << "input gbfile " << gbfile << " open error\n";
		exit(1);
	}

//////making a file name for the new sequence///////////
	char newseqfile[55];
	int i=0;
	while (gbfile[i] != '.' && gbfile[i] != '\0') 
		newseqfile[i] = gbfile[i++];
	memcpy(newseqfile+i, "new.gb", 6);
	newseqfile[i+6]='\0';
///////////////////////////////////////////////////////////

	ofstream NGB(newseqfile);
	
	set<version> oldacc;  //those to be removed from the GB file
	set<version>::iterator sp;
	version vvv; //to store the version info of a GB record
	char line[LLEN+1];

	readAcc(ACC, oldacc);
	cout << "total " << oldacc.size() << " old acc\n";
	ACC.close();
	int newCount = 0;

	GB.getline(line, LLEN);
	while (strncmp(line, "LOCUS", 5)) GB.getline(line, LLEN);
	while (!strncmp(line, "LOCUS", 5) && (oldacc.size()>0)) {
		//temparily store the first few lines in Outbuff
		//cout << "    checking GenBank seq\n" << line << endl;

		//ostrstream Outbuff;      //cause memory leaks
		ostringstream tmpBuffer; //the new C++, no memory leak
		//Outbuff << line << endl;
		tmpBuffer << line << endl;
		while (strncmp(line, "VERSION", 7)) {
			GB.getline(line, LLEN);
			tmpBuffer << line << endl;
		}
		//tmpBuffer << endl;
		vvv.get(line+12);  //get info from version line
		                   //VERSION     AB000095.1  GI:2924600

		sp = oldacc.find(vvv);
		if (sp != oldacc.end()) {  //found in oldacc set
			if (vvv.newer(*sp)) {
				newCount++;
				save(tmpBuffer, GB, NGB);
				//save this one to NGB ofstream
				//cout <<  "new seq found for old ace:" << *sp << endl; //debug
				//cout << vvv << endl;
			}
			else {
				GB.getline(line, LLEN);
				while (strcmp(line, "//")) GB.getline(line, LLEN);
				//discard this entry from GB file
			}
			oldacc.erase(sp);  //both old and current version removed
		}
		else {  // not in the dump Acc.ver list
			newCount++;
			save(tmpBuffer, GB, NGB);
			//new entry not in old list, save this one to NGB ofstream
			//cout << "new entry not in ACEDB:" << vvv << endl;
		}
		GB.getline(line, LLEN);
		while (strncmp(line, "LOCUS", 5) && !GB.eof()) 
			GB.getline(line, LLEN); //in case there are empty lines 
			//between entries
	}
	if (oldacc.size() == 0 && !GB.eof()) { 
	//everything remaining in GB is new relative the ACC given
	//this will make the program run faster, because no need to 
	//check for conditions
		while (!GB.eof()) {
			NGB << line << endl;
			GB.getline(line, LLEN);
		}
	}
	
	char uncheckFile[70];
	strcpy(uncheckFile, gbfile);
	strcat(uncheckFile, "unck.acc");
	ofstream UNCHECK(uncheckFile);
	if (oldacc.size() > 0) { //output this reduced set as unchecked
		cout << oldacc.size() << " Acc.ver entries are not checked\n";
		cout << "they are written to the file: "
		<< uncheckFile << endl;
		sp = oldacc.begin();
		while (sp != oldacc.end()) {
			UNCHECK << *sp << endl;
			sp++;
		}
	}
	else {
		UNCHECK << "\n";
		cout << "nothing to check in the old acc list\a\n";
	}
	cout << "total new sequences: " << newCount << endl;
	//stream book keeping
	GB.close();
	NGB.close();

	return 0;
}

void chopEnds(char* &w)
{
	if (w[0] == '\"') {
		w++;  //to get rid of the first "
		char *p = w + strlen(w) - 1;
		while (isspace(*p)) p--;  // get rid of trailing space
		if (*p == '\"') *p = '\0';
		else *(p+1) = '\0';
	}
	else { // it is not double quoted, nothing needs to be done
	}
}

void readAcc(ifstream &ins, set<version> &vset)
{
	//this function should be able to read two type of input files
	//the one dumped from ACEDB has "xxxxx.n"\t format
	/* a sample version file:
	"AB000095.1"	
	"AB000099.1"	
	"AB000114.1"	
	"AB000115.1"	
	"AB000220.1"	
	"AB000221.1"	
	"AB000263.1"	
	"AB000267.1"	
	"AB000268.1"	
	*/
	//the one ouputed from this program has no double quotes and tab

	char ln[LLEN+1];
	char *p;
	version vv;
	//set<version>::iterator vp;

	ins.getline(ln, LLEN);
	if (ln[0] == '\"') { //from ACEDB tablemaker, double quoted
		while (!ins.eof()) {
			p=ln;
			chopEnds(p);        // removes double quotes 
			vv.quickGet(p);
			vset.insert(vv);
			ins.getline(ln, LLEN);
			//cout << vset.size() << endl;
		}
	}
	else {
		while (!ins.eof()) {
			vv.quickGet(ln);
			vset.insert(vv);
			ins.getline(ln, LLEN);
		}
	}
}

void save(ostringstream &bf, ifstream &ins, ofstream &ous)
{
	ous << bf.str();
	char ln[LLEN+1];
	ins.getline(ln, LLEN);
	while (strcmp(ln, "//")) {
		ous << ln << endl;
		ins.getline(ln, LLEN);
	}
	ous << "//\n";
}
