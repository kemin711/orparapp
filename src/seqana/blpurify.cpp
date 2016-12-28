#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include "alnsummary.h"

#define LINE 90

using namespace std;


//usage blpurify infile
string getQuery(char ln[]);
int resultType(const char ln[]);
void nameOutfile(const char ipf[], const char ext[], char opf[]);

int main(int argc, char *argv[]) 
{
	char infile[50], outfile[55];

	switch (argc) {
		case 1:
			cout << "no input file specified\n";
			exit(1);
			break;
		case 2:
			if (strlen(argv[1]) > 49) {
				cout << "file name too long.  Should be < 49 char\n";
				exit(1); 
			}
			else strcpy(infile, argv[1]);
			break;
		default:
			cout << "too many arguments\n";
			break;
	}
	ifstream INF(infile);
	if (INF.fail()) {
		cout << "input file: " << infile << " open failer\n";
		exit(1);
	}
	nameOutfile(infile, "hi", outfile);
	ofstream OUF(outfile);

	char line[LINE+1];
	string buff, query;
	int scoreStart;

	INF.getline(line, LINE);
	while (!strncmp(line, "BLAST", 5)) {
		buff = line;
		buff += "\n";
		INF.getline(line, LINE);
		while (strncmp(line, "Query=", 6)) INF.getline(line, LINE); 
		//discard Reference: ....
		query = getQuery(line);
		buff += line;
		buff += "\n";
		INF.getline(line, LINE);
		int hit = resultType(line);
		while (!hit) {
			buff += line;
			buff += "\n";
			INF.getline(line, LINE);
			hit = resultType(line);
		}
		if (hit == 2) { //only this one contains hits
			OUF << buff;
			OUF << line << endl << endl;
			scoreStart = strstr(line, "(bit)") - line;
			INF.getline(line, LINE); //reads the empty line
			INF.getline(line, LINE); //reads the first summary line
			alnsummary summln;
			while (strlen(line)>2) {
				summln.read(line, scoreStart);
				cout << summln.subj << summln.s << summln.e << summln.n << endl;
				INF.getline(line, LINE);
			}
		}
		else {
			INF.getline(line, LINE);
			while (!INF.eof() && strncmp(line, "BLAST", 5)) 
				INF.getline(line, LINE);
		}
	}
	return 0;
}

string getQuery(char ln[])
{
	string tmp;
	char *p = strchr(ln, ' ');
	char *pp = strchr(p, ' ');
	tmp = p;
	if (pp) tmp = tmp.substr(0, pp-p);
	return tmp;
}

int resultType(const char ln[])
{
	//return 1, 2, 3 if line is nohit etc or producing ....
	if (strstr(ln, "No hits found")) return 1;
	else if (strstr(ln, "Sequences producing")) return 2;
	else if (strstr(ln, "done  Database")) return 3;
	else return 0;
}

void nameOutfile(const char ipf[], const char ext[], char opf[])
{ 
   char tmp[50];
   int i = 0;
   while (ipf[i] != '.' && ipf[i] != '\0') i++;
   strncpy(tmp, ipf, i); //discard .gb extension
   tmp[i]='\0';

   strcpy(opf, tmp);
	strcpy(opf, ".");
	strcat(opf, ext);
}

