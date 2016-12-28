#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <functional>
#include <set>

#include "bioseq.h"
#include "strformat.h"

using namespace std;
using namespace orpara;

void byWord(const string &oufile, const string &file, const string &wd, const unsigned int length);
list<DNA> byAmbiguous(const string &file, const unsigned int length);

void usage() {
   cerr << "fastadiscard -n file.fas remove sequences with N from file.fas\n"
      << "Options:\n"
      << "   -w word in the title of the sequence\n"
      << "   -i inputfile any fasta file\n"
      << "   -n remove sequences with N in them.\n"
      << "   -o outputfile\n"
      << "   -l lengthcut  discard length below this one\n"
      << "   --byid file containing a list of ids. One sequence id per line\n"
      << "     idline can also be: >id more title words\n"
      << " Right now these operation are not combined\n"
      << " Noe: the GC cutoff must be quoted, and shell use > \n"
      << "      for special purposes\n";
   exit(1);
}

void writeSequences(const list<DNA> &seq, const string &file);
set<string> readIds(const string &idfile);
void discardIds(const string &oufile, const string &infile, const string &idfile, const unsigned int lencut);

string getFileStem(const string &fname) {
   string stem=fname;
   string::size_type x = fname.rfind('/');
   if (x != string::npos) {
      stem = fname.substr(x+1);
   }
   x = stem.rfind('.');
   if (x != string::npos) 
      stem = stem.substr(0, x);
   return stem;
}

string fillSpace(const string &word) {
   string result = word;
   for (size_t i=0; i<word.size(); ++i) {
      if (word[i] == ' ') {
         result[i] = '_';
      }
   }
   return result;
}

int main(int argc, char* argv[]) {
   //string word="Mycobacterium";
   string word, inputFile;
   //inputFile="/remote/DataAnalysis/metag/refseq/silva123.fas";
   unsigned int lengthcut=500;
   bool removeN=false;
   string outputFile, idlistFile;
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-w")) { word=argv[++i]; }
      else if (!strcmp(argv[i], "-i")) { inputFile=argv[++i]; }
      else if (!strcmp(argv[i], "-n")) { removeN=true; }
      else if (!strcmp(argv[i], "-l")) { lengthcut=atoi(argv[++i]); }
      else if (!strcmp(argv[i], "-o")) { outputFile=argv[++i]; }
      else if (!strcmp(argv[i], "--byid")) { idlistFile=argv[++i]; }
      else if (!strcmp(argv[i], "--help")) { usage(); }
      else {
         inputFile=argv[i];
      }
      ++i;
   }

   if (inputFile.empty()) usage();

   if (!idlistFile.empty()) {
      if (outputFile.empty()) {
         outputFile = getFileStem(inputFile) + "_clean.fas";
      }
      discardIds(outputFile, inputFile, idlistFile, lengthcut);
   }

   if (!word.empty()) {
      if (outputFile.empty()) {
         outputFile=getFileStem(inputFile) + "_" + fillSpace(word) + ".fas";
      }
      byWord(outputFile, inputFile, word, lengthcut);
   }
   if (removeN) {
      if (outputFile.empty()) {
         outputFile = getFileStem(inputFile) + "noN.fas";
      }
      list<DNA> unambiguous=byAmbiguous(inputFile, lengthcut);
      writeSequences(unambiguous, outputFile);
   }

   return 0;
}

void writeSequences(const list<DNA> &seq, const string &file) {
   ofstream ouf(file);
   for (auto it=seq.begin(); it != seq.end(); ++it) {
      ouf << *it;
   }
   cout << seq.size() << " sequences written to " << file << endl;
}

list<DNA> byAmbiguous(const string &file, const unsigned int length) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + file);
   }
   cerr << "input file: " << file << endl;
   list<DNA> result;
   int cnt=0;
   int every=10;
   while (dna.read(inf)) {
      ++cnt;
      if (dna.length() > length && dna.noAmbiguous()) {
         result.push_back(dna);
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
   return result;
}

void byWord(const string &oufile, const string &file, const string &wd, const unsigned int length) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + file);
   }
   cerr << "input file: " << file << endl;
   string lowerwd=getLower(wd);
   int cnt=0;
   int every=10;
   ofstream ouf(oufile);
   while (dna.read(inf)) {
      ++cnt;
      if (dna.length() < length) continue;
      if (getLower(dna.getTitle()).find(lowerwd) == string::npos) {
         ouf << dna;
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
}

//>JP829389.1.1014 Bacteria;Firmicutes;Bacilli;Bacillales;Paenibacillaceae;Cohnella;Triticum aestivum (bread wheat)
//>JW028727.1.1036 Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;Triticum aestivum (bread wheat)
//>HX130966.2.728 Bacteria;Actinobacteria;Actinobacteria;Micrococcales;Microbacteriaceae;Amnibacterium;Triticum aestivum (bread wheat)
//>HP624559.7.961 Eukaryota;Opisthokonta;Holozoa;Metazoa (Animalia);Eumetazoa;Bilateria;Arthropoda;Hexapoda;Insecta;Pterygota;Neoptera;Hemiptera;Triticum aestivum (bread wheat)
//
//Id file can be just fastaa headers without any processing.
set<string> readIds(const string &idfile) {
   ifstream inf(idfile);
   string line;
   string::size_type s;
   set<string> result;
   while (!getline(inf, line).eof()) {
      if (line[0] == '>') line = line.substr(1);
      s=line.find(' ');
      if (s != string::npos) {
         result.insert(line.substr(0, s));
      }
      else result.insert(line);
   }
   cerr << result.size() << " ids read\n";
   return result;
}

void discardIds(const string &oufile, const string &infile, const string &idfile, const unsigned int lencut) {
   set<string> ids = readIds(idfile);
   DNA dna;
   ifstream inf(infile);
   if (inf.fail()) {
      throw runtime_error("cannot read " + infile);
   }
   ofstream ouf(oufile);
   if (ouf.fail()) throw runtime_error("failed to open file for write: " + oufile);
   int cnt=0;
   int every=100;
   while (dna.read(inf)) {
      ++cnt;
      if (dna.length() < lencut) continue;
      if (ids.find(dna.getName()) == ids.end()) {
         ouf << dna;
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << endl;
         every *= 2;
      }
   }
}

