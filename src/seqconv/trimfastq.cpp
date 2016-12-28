#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>

#include <fastq.h>
#include <stddev.h>
#include <strformat.h>

using namespace std;
using namespace orpara;

void usage();
/**
 * @return the oputput file name with a pattern infile.stem_trimWC.fastq.
 *         where W and C are two integers of window and cutoff values.
 */
string nameOutputFile(const string &infile, const int w, const int c);
string nameStatFile(const string &infile, const int w, const int c);
string getFileStem(const string &infile);
void writeStatHeader(ostream &ous) {
   ous << "state\taverage\tstddev\tn\n";
}
void writeStatValue(ostream &ous, const stddev &stat, const string &msg) {
   ous << msg << '\t' << stat.getMean() << '\t' << stat.getStd() << '\t' << stat.getCount() << endl;
}

int main(int argc, char* argv[]) {
   //string infile("bar4.fastq");
   //string outfile("bar4trim.fastq");
   string infile, outfile;
   int window=6, cutoff=20;
   unsigned int lengthcut = 17;
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "-w") window = atoi(argv[++i]);
      else if (string(argv[i]) == "-c") cutoff = atoi(argv[++i]);
      else if (string(argv[i]) == "-l") lengthcut = atoi(argv[++i]);
      else if (string(argv[i]) == "--help" || argv[i][0]=='?') {
         usage(); return 1;
      }
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
      return 1;
   }
   if (outfile.empty()) {
      outfile = nameOutputFile(infile, window, cutoff);
   }

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      return 1;
   }
   // output file stream
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writing " << endl;
      return 1;
   }
   cout << "Using w c: " << window << " " << cutoff << endl;
   stddev avgstdBefore, avgstdAfter;
   Fastq fasq;
   int count=0;
   int goodcnt = 0;
   while (fasq.read(inf)) {
      ++count;
      avgstdBefore(fasq.length());
      //cerr << count << " th seq\n";
      /*
      if (fasq.trimLowq()) {
         ouf << "\n**Trimmed\n";
      }
      */
      fasq.trimLowq(window, cutoff);
      if (fasq.length() > lengthcut) {
         ++goodcnt;
         fasq.write(ouf);
         avgstdAfter(fasq.length());
      }
      //if (count > 100) break;
   }
   string statFile = nameStatFile(infile, window, cutoff);

   ofstream sout(statFile.c_str());
   sout << goodcnt << " good seq out of " << count << " output written to " << outfile << endl;
   writeStatHeader(sout);
   writeStatValue(sout, avgstdBefore, "before");
   writeStatValue(sout, avgstdAfter, "after");
   sout.close();
   
   cout << goodcnt << " good seq out of " << count << " output written to " << outfile << endl;

   return 0;
}

void usage() {
   cout << "Usage: trimfastq --help will print help messasge\n"
      << "  trimfastq ? will also print help message\n"
      << "  You can also do this: trimfastq infile outfile\n"
      << "  Options:\n"
      << "      -i input fastq file\n"
      << "      -o output fastq file. If not given will be constructed by this program\n"
      << "      -w trimming window parameter. Default 5. Try 4-7 \n"
      << "      -c trimming cut parameter. Default 20. Try 18-30. w=5,c=20 gave best results for one batch.\n"
      << "      -l length of fastq to save. Default 17\n"
      << " Except for -i, all other options are optional. The input file \n" 
      << " Can also be given directly without the -i options\n";
}



/**
 * @return the trim output file name for a given input file name
 *    and parameters used.
 */
string nameOutputFile(const string &infile, const int w, const int c) {
   return getFileStem(infile) + "_trim" + itos(w) + itos(c) + ".fastq";
};

/** 
 * given input fastq file name
 * @return the stat file name.
 */
string nameStatFile(const string &infile, const int w, const int c) {
   string stem = getFileStem(infile);
   return stem + "_trim" + itos(w) + itos(c) + ".stat";
};

string getFileStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      stem = infile.substr(0, i);
   }
   //cerr << "file stem of " << infile << " is " << stem << endl;
   return stem;
}
