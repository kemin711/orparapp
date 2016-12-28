#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <functional>
#include <string>
#include <set>
#include <cstring>

#include "bioseq.h"

using namespace std;
using namespace orpara;

void usage() {
   cerr << "fastapickpack -i idname_file -n seqfilename_file -o Bacgenomic\n"
      << "Options:\n"
      << "   -i file name containing the sequence ids wanted\n"
      << "   -n FILE containing the input file name\n"
      << "   -o output file name prefix.  Multiple file will produced\n"
      << "      file size is limited to 4 billion bases\n"
      << "   -h or --help print this message\n";
   exit(1);
}

set<string> readSet(const string &fname);
vector<string> readVector(const string &fname);
void pick(const string &fname, set<string> &ids, list<DNA> &store, unsigned long int &nb);
void writeSequences(const list<DNA> &seq, const string &file);
void work(const vector<string> &files, set<string> &ids, const string &prefix);

/**
 * pick from a file containing a list of sequence ids
 * and files matching certain pattern or all files in 
 * a directory.  Then write the squences into filess
 * with certain limits to the total number of bases
 * in each file. Say 4 billion bases per file.
 */
int main(int argc, char* argv[]) {
   // work in the current directory
   // one id per line or separate by space
   string inputidsFile="reprefseqgid.tab";
   string infnameFile="seqfilenames.txt";
   // output prefix
   string outprefix="bactgenomic";
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-i")) { inputidsFile=argv[++i]; }
      else if (!strcmp(argv[i], "-n")) { infnameFile=argv[++i]; }
      else if (!strcmp(argv[i], "-o")) { outprefix=argv[++i]; }
      else if (!strcmp(argv[i], "-h")) { usage(); }
      else {
         inputidsFile=argv[i];
      }
      ++i;
   }

   vector<string> inputFiles=readVector(infnameFile);
   set<string> inputIds=readSet(inputidsFile);
   work(inputFiles, inputIds, outprefix);

   return 0;
}

void work(const vector<string> &files, set<string> &ids, const string &prefix) {
   static unsigned long int limit = 4000000000;
   unsigned long int nb=0;
   list<DNA> store;
   int bucket=1;
   for (size_t i=0; i<files.size(); ++i) {
      cerr << "working on " << files[i] << endl;
      pick(files[i], ids, store, nb);
      if (nb > limit) {
         string oufname = prefix + to_string(bucket++) + ".fas";
         writeSequences(store, oufname);
         store.clear();
         nb=0;
      }
   }
}

void pick(const string &fname, set<string> &ids, list<DNA> &store, unsigned long int &nb) {
   ifstream inf(fname);
   DNA dna;
   set<string>::iterator it;
   while (dna.read(inf)) {
      it = ids.find(dna.getName());
      if (it != ids.end()) {
         nb += dna.length();
         store.push_back(dna);
         ids.erase(it);
      }
   }
}

set<string> readSet(const string &fname) {
   set<string> res;
   ifstream inf(fname);
   string item;
   while (inf >> item) {
      res.insert(item);
   }
   cout << res.size() << " items stored in set\n";
   return res;
}

vector<string> readVector(const string &fname) {
   vector<string> res;
   ifstream inf(fname);
   string item;
   while (inf >> item) {
      res.push_back(item);
   }
   cout << res.size() << " items stored in vector\n";
   return res;
}

void writeSequences(const list<DNA> &seq, const string &file) {
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("failed to open file: " + file);
   }
   for (auto it=seq.begin(); it != seq.end(); ++it) {
      ouf << *it;
   }
   cout << seq.size() << " sequences written to " << file << endl;
}

