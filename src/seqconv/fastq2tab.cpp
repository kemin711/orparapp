#include <fstream>
#include <iostream>
#include <cstdlib>
#include <map>

#include <fastq.h>

using namespace std;
using namespace orpara;

void usage() {
   cerr << "Usage: fastq2tab input-fastqfile outputfile\n"
      << "Options:\n"
      << "  -i inputFastqFile input fastq file\n"
      << "  -o outputFile table format\n"
      << "  -s output sequence. This is a flag, default not output sequence\n"
      << "Convert a fastq file into a table format for loading\n"
      << "into relational database. You must have very large disk\n"
      << "for your database installation if you want to store\n"
      << "the actual sequence\n";
   exit(1);
}

string buildOutputName(const string &infile, const string &suffix);

int main(int argc, char* argv[]) {
   string infile, outfile;
   string lenfile; //("length.dis");
   bool outputSeq=false;
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "-s") outputSeq = true;
      else if (string(argv[i]) == "--help") usage();
      else {
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
   if (outfile.empty()) {
      outfile = buildOutputName(infile, ".tab");
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
   cout << "Converting fastq file: " << infile << " to tabular format\n";
   Fastq fasq;
   while (fasq.read(inf)) {
      ouf << fasq.getName() << "\t" << fasq.getTitle()
         << "\t" << fasq.length();
      if (outputSeq) ouf << fasq.getSequence();
      ouf << endl;
   }
   inf.close();
   ouf.close();
   cout << "tabular file: " << outfile << endl;

   return 0;
}

string buildOutputName(const string &infile, const string &suffix) {
   string oufile;
   string::size_type j = infile.rfind('.');
   if (j != string::npos) oufile = infile.substr(0,j);
   else oufile = infile;
   oufile += suffix;
   return oufile;
}

