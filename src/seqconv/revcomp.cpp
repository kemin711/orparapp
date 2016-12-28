#include <iostream>
#include <fstream>

#include <bioseq.h>

using namespace std;
using namespace orpara;

string makeOutputFileName(const string &inf);
void usage() {
   cerr << "Usage: revcomp -o outfile input.fasta\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   int i = 1;
   string infile, outfile;
   while (i < argc) {
      if (string(argv[i]) == "-o") {
         outfile = string(argv[++i]);
      }
      else {
         infile = string(argv[i]);
      }
      ++i;
   }
   if (infile.empty()) usage();

   DNA dna;
   dna.read(infile);
   dna.revcomp();
   if (outfile.empty()) {
      outfile = makeOutputFileName(infile);
   }
   ofstream of(outfile.c_str());
   of << dna;
   cout << "reverse completed sequence written to file: " << outfile << endl;

   return 0;
}

string makeOutputFileName(const string &inf) {
   string::size_type i;
   string ofile=inf;
   if ((i=inf.rfind('.')) != string::npos) {
      ofile = inf.substr(0, i) + "rc.fas";
   }
   return ofile;
}
