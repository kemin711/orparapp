#include <iostream>
#include <fstream>

#include <strformat.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;

void extractHeader(const string &fasFile, const string &newfile, const string &header);
string getStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if ( i != string::npos) {
      stem = infile.substr(0, i);
   }
   return stem;
}

void usage() {
   cerr << "parsencbi16s -h header -o shortid.fas inputfile\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   int i=1;
   string inputfile, outputfile, headerfile;
   //string inputfile="Bacteria16S";
   //string outputfile="Bacteria16Ssid.fas";
   //string headerfile="Bacteria16Sheader.tab";
   while (i < argc) {
      if (!strcmp(argv[i], "-t")) {
         headerfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outputfile = argv[++i];
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
   extractHeader(inputfile, outputfile, headerfile);
}

string getAccession(const string &nrid) {
   string::size_type i = nrid.find("ref|");
   if (i == string::npos) {
      throw runtime_error("nr id format wrong" + nrid);
   }
   string acc=nrid.substr(i+4);
   if (acc[acc.length()-1] == '|') {
      acc.resize(acc.length()-1);
   }
   return acc;
}

string cleanName(const string &input) {
   string name=input;
   if (input[0] == '[') {
      name=input.substr(1);
      name.erase(name.find(']'), 1);
   }
   return name;
}

bool endWithQuotedSpecies(const string &sp) {
   string::size_type i=sp.find('(');
   string::size_type j=sp.rfind(')');
   if (j == sp.length()-1) {
      if (isScientificSpeciesName(sp.substr(i+1,j-i-1))) {
         return true;
      }
   }
   return false;
}

bool isSymbion(const string &species, const string &strain) {
   if (species.find(" sp. ") != string::npos) {
      string::size_type i=species.find('(');
      string::size_type j=species.rfind(')');
      if (isScientificSpeciesName(species.substr(i+1,j-i-1)) 
            && (isupper(strain) || wc(strain) == 1)) {
         return true;
      }
   }
   else if (endWithQuotedSpecies(species)) {
      return true;
   }
   return false;
}

tuple<string, string, bool> breakupHeader(const string &hd) {
   string tmp=hd;
   bool partial=false;
   string::size_type i = hd.rfind(", complete sequence");
   if (i != string::npos) {
      tmp = hd.substr(0, i);
   }
   else if ((i = hd.rfind(", cpmplete sequence")) != string::npos) {
      // spelling error
      tmp = hd.substr(0, i);
   }
   else {
      i = hd.rfind(", partial sequence");
      if (i != string::npos) {
         tmp = hd.substr(0,i);
         partial = true;
      }
      //else {
         //throw runtime_error("failed to parse header " + hd);
      //   cerr << "failed to find partial/complete label in:\n" 
      //      << hd << " |END" << endl;
      //}
   }
   i = tmp.find("16S ribosomal RNA");
   if (i != string::npos) {  // useless lebel
      if (i+17 == tmp.length()) {
         // ends with this, remove it
         tmp.resize(tmp.length()-17);
      }
      else {
         tmp.erase(i, 17);
      }
   }
   tmp = tmp.substr(0, i-1);
   i = tmp.find(" strain ");
   string strain;
   if (i != string::npos) {
      strain = tmp.substr(i + 8);
      tmp = tmp.substr(0, i);
      if ((i=tmp.rfind(" str. ")) != string::npos) { // duplicated strain, remove 
         tmp = tmp.substr(0, i);
      }
      if (strain.substr(0,2) == ": ") {
         strain = strain.substr(2);
      }
   }
   if (strain.empty()) {
      i = tmp.find(" isolate ");
      strain = tmp.substr(i+9);
      tmp = tmp.substr(0, i);
   }
   // Mycoplasma haemofelis str. Langford 1 strain Langford 1 16S ribosomal RNA, complete sequence
   // now trimmmed to |-strain->END
   string strstrain = "str. " + strain;
   if (endwith(tmp, strstrain)) {
      tmp.resize(tmp.length() - strain.length() - 6);
   }
   else if (endwith(tmp, strain)) {
      tmp.resize(tmp.length() - strain.length() - 1);
   }
   else if (isSymbion(tmp, strain)) {
      // Blattabacterium sp. (Blaberus giganteus) strain BGIGA
      //cerr << tmp << " is symbion\n";
   }
   else {
      if (!isScientificSpeciesName(tmp) && !endwith(tmp, " sp.")) {
         string clean = cleanName(tmp);
         if (!isScientificSpeciesName(cleanName(clean))) {
            string tailwd = lastword(clean);
            if (tailwd.rfind(strain) != string::npos) {
               strain = tailwd;
               clean.resize(clean.length() - tailwd.length());
               tmp = clean;
            }
            else {
               //cerr << hd << endl;
               //cerr << "no redundant strain in title: " << tmp << endl;
            }
         }
         else tmp=clean;
      }
   }
   // no extract strains for Symbion
   //cout << hd << endl;
   //cout << tmp << " | " << strain << " | " << partial << endl;
   return make_tuple(tmp, strain, partial);
}

void extractHeader(const string &fasFile, const string &newfile, const string &header) {
   bioseq seq;
   ifstream inf(fasFile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + fasFile);
   }
   ofstream ouf(newfile);
   ofstream oufh(header);
   while (seq.read(inf)) {
      //cout << seq.getName() << endl
      //   << seq.getTitle() << endl;
      seq.setName(getAccession(seq.getName()));
      //cout << seq.getName() << endl;
      tuple<string, string, bool> taxon = breakupHeader(seq.getTitle());
      oufh << seq.getName() << '\t' << get<0>(taxon) << '\t' << get<1>(taxon) << endl;
      ouf << seq;
   }
   cout << "fast file with shorter id written to " << newfile << endl;
}

