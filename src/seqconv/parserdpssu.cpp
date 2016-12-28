#include <iostream>
#include <fstream>

#include <strformat.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;

void extractHeader(const string &fasFile, const string &newfile, const string &header, unsigned int lencut);
string getStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if ( i != string::npos) {
      stem = infile.substr(0, i);
   }
   return stem;
}

void usage() {
   cerr << "parserdpssu -h header -o shortid.fas inputfile\n"
      << "parses small subunit including 16S and 18S ribisomal rRNA file\n";
   exit(1);
}


int main(int argc, char* argv[]) {
   int i=1;
   string inputfile, outputfile, headerfile;
   //string inputfile="Bacteria16S";
   //string outputfile="Bacteria16Ssid.fas";
   //string headerfile="Bacteria16Sheader.tab";
   unsigned int lengthcut=700;
   while (i < argc) {
      if (!strcmp(argv[i], "-t")) {
         headerfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outputfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-l")) {
         lengthcut = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
         usage();
      }
      else {
         inputfile = argv[i];
      }
      ++i;
   }
   if (inputfile.empty()) usage();
   if (outputfile.empty()) {
      outputfile = getStem(inputfile) + "cut.fas";
   }
   if (headerfile.empty()) {
      headerfile = getStem(inputfile) + "header.tab";
   }
   extractHeader(inputfile, outputfile, headerfile, lengthcut);
}

/**
 * return species, strain pair
 */
pair<string, string> breakupHeader(const string &hd) {
   string species, strain;
   string::size_type i = hd.find('\t');
   if (i == string::npos) {
      throw runtime_error("no tab");
   }
   species = hd.substr(0, i);
   //string lineaqge = hd.substr(i+1); // don't use this info yet
   i = species.find("; ");
   if (i != string::npos) {
      strain = species.substr(i+2);
      species = species.substr(0, i);
   }
   return make_pair(species, strain);
}

void extractHeader(const string &fasFile, const string &newfile, const string &header, unsigned int lencut) {
   DNA seq;
   ifstream inf(fasFile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + fasFile);
   }
   ofstream ouf(newfile);
   ofstream oufh(header);
   int cnt = 0;
   while (seq.read(inf)) {
      if (seq.getTitle().find("synthetic") != string::npos
            || seq.length() < lencut) continue;
      pair<string, string> taxon = breakupHeader(seq.getTitle());
      oufh << seq.getName() << '\t' << taxon.first << '\t' << taxon.second << endl;
      ouf << seq;
      ++cnt;
   }
   cout << cnt << " fast sequences longer than " << lencut 
      << " written to " << newfile << endl;
}
