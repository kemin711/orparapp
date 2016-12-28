#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

/*this program will take a dbEST sequence file
sequences from one organism such as Fugu rubripes or human will be 
separated from other sequences 
*/

using namespace std;

void getInput(char inf[]);
//void nameFile(const char ipf[], string &s, string &nons, const string &name);

int main(int argc, char *argv[])
{
	char infile[51], organism[100]="Homo sapiens"; 
  const char dir[] = "/data/data7/inputs/dbEST/2001-02-13";
  const int LINE = 90;  //length of line
	//string org("Homo sapiiens"); //default
	//int show = 1000;

	switch (argc) {
	case 1:
		getInput(infile);
		break;
	case 2:
		strcpy(infile, argv[1]);
		break;
	default:
		int i = 1;
		while (i<argc) {
			if (!strcmp(argv[i], "--org") || !strcmp(argv[i], "-o")) {
				i++;
				strcpy(organism, argv[i]);
			}
			else strcpy(infile, argv[i]);
			i++;
		}
		break;
	}

  /// setting up input and output files  /////////////
  char infileFullName[150];
  strcpy(infileFullName, dir);
  strcat(infileFullName, "/");
  strcat(infileFullName, infile);

	ifstream INF(infileFullName);  //input file containing the EST files
	if (INF.fail()) {
		cerr << "Input file \"" << infileFullName << "\" not found\n";
		return 1;
	}

  char outfileThisOrg[150], outfileOtherOrg[150];
  strcpy(outfileThisOrg, infile);
  strcat(outfileThisOrg, ".human");
  strcpy(outfileOtherOrg, infile);
  strcat(outfileOtherOrg, ".other");
	ofstream THISORG(outfileThisOrg);  
	ofstream REST(outfileOtherOrg);  

	int seqCnt = 0, thisSeqCnt = 0, otherSeqCnt=0;

	//string ln, buff, acc;
  char ln[LINE+1], *p;

	INF.getline(ln, LINE);

  while (strncmp(ln, "IDENTIFIERS", 11)) INF.getline(ln, LINE);
  //while (ln.substr(0, 11) != "IDENTIFIERS") getline(INF, ln);

  while (!INF.eof() && !strncmp(ln, "IDENTIFIERS", 11) ) {
    seqCnt++;
    ostringstream sbff; //out buffer
    sbff << endl;
    while (strncmp(ln, "Organism", 8) && !INF.eof()) {
      sbff << ln << endl;  //ln contains the LOCUS line
      INF.getline(ln, LINE);
    }
    sbff << ends;
    p = ln + 16;
    if (!strcmp(p, "Homo sapiens")) {
      thisSeqCnt++;
      THISORG << sbff.str();
      while (strcmp(ln, "||") && !INF.eof()) {
        THISORG << ln << endl;
        INF.getline(ln, LINE);
      }
      THISORG << ln << endl;
    }
    else { //sequence of other organism
      otherSeqCnt++;
      REST << sbff.str();
      while (strcmp(ln, "||") && !INF.eof()) {
        REST << ln << endl;
        INF.getline(ln, LINE);
      }
      REST << ln << endl;
    }
    INF.getline(ln, LINE);
    while (!INF.eof() && strncmp(ln, "IDENTIFIERS", 11)) 
      INF.getline(ln, LINE);
	}

	cout << "Total " << seqCnt << " sequences\n";
	cout << "Total " << thisSeqCnt << " " << organism << " sequences\n";
	cout << "Total " << otherSeqCnt << " Non-" << organism << " sequences\n";
	
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

