#include <iostream>
#include <fstream>

#include <strformat.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;

void extractHeader(const string &fasFile, const string &newfile, const string &header,
      unsigned int lencut);
string getStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if ( i != string::npos) {
      stem = infile.substr(0, i);
   }
   return stem;
}

void usage() {
   cerr << "parsencbissu -h header -o shortid.fas inputfile\n"
      << "parses small subunit including 16S and 18S ribisomal rRNA file\n"
      << "Most sequences are full length or close to full length so\n"
      << "there is no need for length filtering\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   int i=1;
   string inputfile, outputfile, headerfile;
   //inputfile="/home/users/zhouk15/overflow/refseq/hmp/hmp_ncbi.fasta";
   inputfile="/home/kzhou/work/metag/refseq/hmp_ncbi.fasta";
   //string outputfile="Bacteria16Ssid.fas";
   //string headerfile="Bacteria16Sheader.tab";
   unsigned int lengthcut = 700;
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

string extractTag(const string &str) {
   static vector<string> tags = {"gb", "emb", "dbj" };
   string id;
   for (unsigned i = 0; i<tags.size(); ++i) {
      string::size_type x = str.find(tags[i] + "|");
      if (x != string::npos) {
         id = str.substr(x+tags[i].length()+1);
         string::size_type y = id.find('|');
         if (y != string::npos) {
            id = id.substr(0, y);
         }
         break;
      }
   }
   if (id.empty()) {
      throw runtime_error("Failed to extract GenBank id from " + str);
   }
   return id;
}

string getAccession(const string &nrid) {
   string acc;
   string::size_type i = nrid.find("ref|");
   if (i != string::npos) {
      acc=nrid.substr(i+4);
      if (acc[acc.length()-1] == '|') {
         acc.resize(acc.length()-1);
      }
   }
   else { // not refseq, normal GB id
      acc = extractTag(nrid);
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

void removeSLabel(string &title) {
   static vector<string> word = {" 18S ribosomal RNA gene", " 18S ribosomal RNA", " 16S ribosomal RNA gene", " 16S ribosomal RNA", " 16S rRNA gene, ", "gene for 16S rRNA" };
   bool found = false;
   // remove from the end
   for (unsigned i = 0; i < word.size(); ++i) {
      if (endwith(title, word[i])) {
         title.resize(title.length() - word[i].length());
         found = true;
         break;
      }
   }
   //if (endwith(title, ",")) {
   //   title.resize(title.length()-1);
   //}
   // remove from the middle
   if (!found) {
      for (unsigned i = 0; i<word.size(); ++i) {
         string::size_type x = title.rfind(word[i]);
         if (x != string::npos) {
            title.erase(x, word[i].length());
            found = true;
            break;
         }
      }
   }
   //if (!found) {
   //   cerr << title << " has no ribosomal RNA lable\n";
   //}
}

tuple<string, string, bool> breakupHeader(const string &hd) {
   string tmp=hd;
   bool partial=false;
   string::size_type i;
   // spelling error, get partial or complete information
   if (endwith(hd, ", complete sequence") || endwith(hd, ", cpmplete sequence")) {
      tmp.resize(tmp.length() - 19);
   }
   else if ((i = tmp.rfind(" partial 16S rRNA gene, ")) != string::npos) {
      tmp.erase(i, 23);
      partial = true;
   }
   else if (endwith(tmp, ", partial sequence")) {
      tmp.resize(tmp.length()-18);
      partial = true;
   }
   else if ((i = tmp.rfind(", partial sequence, ")) != string::npos) {
      tmp.erase(i, 19);
      partial = true;
   }
   else if ((i=tmp.rfind(" partial 16S ribosomal RNA gene, ")) != string::npos) {
      tmp.erase(i, 32);
      partial = true;
   }
   else {
      cerr << "no paritial or complete label found: " << tmp << endl;
   }
   removeSLabel(tmp);
   string strain;
   i = tmp.find(" type strain ");
   if (i != string::npos) {
      strain = tmp.substr(i+13);
      tmp = tmp.substr(0, i);
   }
   if (strain.empty()) {
      i = tmp.find(" strain ");
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
   }
   if (strain.empty()) {
      i = tmp.find(" isolate ");
      if (i != string::npos) {
         strain = tmp.substr(i+9);
         tmp = tmp.substr(0, i);
      }
   }
   if (strain.empty()) {
      i = tmp.find(" strain: ");
      if (i != string::npos) {
         strain = tmp.substr(i+9);
         tmp = tmp.substr(0, i);
      }
   }
   if (strain.empty()) {
      vector<string> word = split(tmp, ' ');
      if (word.size() > 2) {
         string candidateSP = word[0] + " " + word[1];
         if (isScientificSpeciesName(candidateSP)) {
            strain = tmp.substr(candidateSP.length() + 1);
            tmp = candidateSP;
         }
      }
   }
   if (strain.empty()) {
      i = tmp.find(" clone ");
      if (i != string::npos) {
         strain = tmp.substr(i+7);
         tmp = tmp.substr(0, i);
      }
   }
   if ((i=tmp.find(" ATCC ")) != string::npos) {
      if (strain.empty()) {
         strain = tmp.substr(i+1);
      }
      else {
         strain += tmp.substr(i);
      }
      tmp = tmp.substr(0, i);
   }

   // Mycoplasma haemofelis str. Langford 1 strain Langford 1 16S ribosomal RNA, complete sequence
   // now trimmmed to |-strain->END
   if (!strain.empty()) {
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
            if (!isScientificSpeciesName(clean)) {
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
   }
   // no extract strains for Symbion
   //cout << hd << endl;
   //cout << tmp << " | " << strain << " | " << partial << endl;
   return make_tuple(tmp, strain, partial);
}

void extractHeader(const string &fasFile, const string &newfile, const string &header,
      unsigned int lencut) 
{
   bioseq seq;
   ifstream inf(fasFile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + fasFile);
   }
   ofstream ouf(newfile);
   ofstream oufh(header);
   while (seq.read(inf)) {
      //cout << seq.getName() << endl << seq.getTitle() << endl;
      if (seq.length() < lencut) continue;
      seq.setName(getAccession(seq.getName()));
      //cout << seq.getName() << endl;
      tuple<string, string, bool> taxon = breakupHeader(seq.getTitle());
      oufh << seq.getName() << '\t' << get<0>(taxon) << '\t' << get<1>(taxon) << endl;
      ouf << seq;
   }
   cout << "fast file with shorter id written to " << newfile << endl;
}

