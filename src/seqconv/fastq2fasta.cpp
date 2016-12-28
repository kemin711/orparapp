#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>

#include <fastq.h>
#include <stddev.h>

using namespace std;
using namespace orpara;

void usage() {
   cerr << "fastq2fasta seq.fastq\n"
      << "Options\n"
      << "   -i inputFastq file\n"
      << "   -o output_file_name\n"
      << "   -h or --help display help message\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   string infile, outfile;
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
         usage();
      else {
         if (argv[i][0] == '-') {
            cerr << "wrong options: " << argv[i] << endl;
            usage();
         }
         infile=argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            ++i;
            outfile = argv[i];
         }
      }
      ++i;
   }

   if (infile.empty()) {
      usage();
   }

   unsigned int j;
   if (outfile.empty()) {
      j = infile.find('.');
      if (j != string::npos) outfile = infile.substr(0,j);
      else outfile = infile;
      outfile += ".fas";
   }

   ifstream inf(infile.c_str());
   ofstream ouf(outfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      return 1;
   }
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
   }
   cout << "Converting fastq: " << infile << " into fasta: " << outfile << endl;

   stddev avgstd;
   Fastq fasq;
   j = 0;
   while (fasq.read(inf)) {
      fasq.writeFasta(ouf);
      ++j;
   }
   cout << j << " fastq sequences converted\n";

   return 0;
}

