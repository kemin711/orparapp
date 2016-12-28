#include <iostream>
#include <fstream>

#include <bioseq.h>
#include <strformat.h>

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
   cerr << "parsepatriacssu -h header -o shortid.fas inputfile\n"
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
      outputfile = getStem(inputfile) + "sid.fas";
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
   string::size_type i = hd.find('[');
   string::size_type j = hd.rfind(']');
   if (i == string::npos) {
      throw runtime_error("square bracket: " + hd );
   }
   species = hd.substr(i+1, j-i-1);
   vector<string> word = split(species, ' ');
   if (word.size() > 2) {
      if (isupper(word[0][0]) && islower(word[1][0])) {
         string newspecies = word[0] + " " + word[1];
         strain = species.substr(newspecies.length()+1);
         species = newspecies;
      }
   }
   return make_pair(species, strain);
}

bool notSSU(const string &hd) {
   if (hd.find("tRNA") != string::npos 
         || hd.find("Large Subunit") != string::npos) {
      return true;
   }
   if (hd.find("Small Subunit") != string::npos
         || hd.find("16S ribosomal") != string::npos)
      return false;
   // if it is not labeled as SSU then it is something else.
   // cerr << hd << " not sure it is SSU\n";
   return false;
}

string extractId(const string &oldid) {
   string newid=oldid;
   if (oldid[oldid.size()-1] == '|')
      newid=oldid.substr(0, oldid.size()-1);
   string::size_type i = newid.rfind('|');
   if (i == string::npos) 
      throw runtime_error("PATRIC id format wrong: " + oldid);
   return newid.substr(i+1);
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
            || seq.length() < lencut
            || notSSU(seq.getTitle())) {
         continue;
      }
      seq.setName(extractId(seq.getName()));
      pair<string, string> taxon = breakupHeader(seq.getTitle());
      oufh << seq.getName() << '\t' << taxon.first << '\t' << taxon.second << endl;
      ouf << seq;
      ++cnt;
   }
   cout << cnt << " fast sequences longer than " << lencut 
      << " written to " << newfile << endl;
}

