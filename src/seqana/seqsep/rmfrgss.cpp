#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cstring>

/* This function removes all GSS sequence from Fugu database in the GB format
input file should be named fugu.gb
the gi number for GSS and nonGSS will be outputed
*/
//extern long _ftype;   // for macintoch text file compatibility

using namespace std;

int main()
{
	ifstream infs("fugu.gb");  //input file
	if (infs.fail()) {
		cout << "you must name your input file \"fugu.gb" << endl;
		exit(1);
	}
//	_ftype = 'TEXT';      // on macitoch machines
	ofstream outfs("fugunogss.gb");  //result file without GSS entries
	ofstream logs("logfile");  // first two lines of saved records
	ofstream rmlog("removedseqs");  //not included sequences
	ofstream badlog("badlog.rec");  //exceptions
	ofstream gssgi("GSS.gi");       //GSS gi numbers 
	//ofstream nongssgi("nonGSS.gi"); //nonGSS gi numbers, not collected at this version
	
	char line[130], locus[130];
	int gi, total;
	unsigned int numgb = 0, numgss = 0, numbad = 0;
	const char gssdef[] = "GSS sequence,";  
	infs.getline(locus, 129);  //reading the locus line
	infs.getline(line, 129);   //reading the definition line
	cout << "Removing GSS sequences ...........\n";
	while (!infs.eof())  {
		if (strstr(line, gssdef)) {
			numgss++;
			rmlog << locus << endl;
			infs >> line;
			while (strcmp(line, "NID")) infs >> line;
			infs.ignore(10);  //trying to get the gi number from "NID         g3270896"
			infs >> gi;
			gssgi << gi << "\t"; 
			if (numgss%20 == 0) gssgi << endl;  //output 20 gi in one line
			infs.getline(line, 129);
			while (strcmp(line, "//"))  //get rid of the rest of the file.   
				infs.getline(line, 129);
		}
		else  {  //not GSS seq
			numgb++;
			outfs << locus << endl << line << endl;
			//cout << locus << endl << line << endl;  //for debug 
			logs << locus << endl << line << endl;
			infs.getline(line, 129);
			while ( strncmp(line, "//", 2) )  {
				outfs << line << endl;
				infs.getline(line, 129);
			}
			outfs << "//\n";
		}
		total = numgss + numgb;
		if ( total%1000 == 0 ) cout << total << " processed\n";  //for moniter 
		if (!infs.eof()) infs.getline(locus, 129);  
		//read Lucus line only if not reaching the end yet
		
		while ( !strstr(locus, "LOCUS") && !infs.eof()) {
			badlog << locus << endl;
			numbad++;
			infs.getline(locus, 129);
			cout << "Not retrieved sequence from GB:\n";
			cout << locus << endl;
		}//read to the first LOCUS line
		if (!infs.eof()) 	infs.getline(line, 129);
		
	}
	infs.close();
	outfs.close();
	logs << "total " << numgb << " sequences\n";
	cout << "Total number of Non-GSS sequences processed: " 
			<< numgb << endl;
	logs.close();
	rmlog << "Total number of GSS sequences " << numgss << endl;
	rmlog.close();
	badlog << "Total bad sequences " << numbad << endl;
	badlog.close();
	return EXIT_SUCCESS;
}

