#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <string>
#include <iostream>
#include <fstream>
#include <map>

#include <fastq.h>
#include <bioseq.h>
#include <dynalnt.h>

using namespace orpara;

bool isSimple(const bioseq &bs);
bool isArtifact(const string &seq);
string nameOutfile(const string &infile);
string nameStatfile(const string &infile);
void printLengthStat(ostream &ous, const map<int,int> &lstat);
/**
 * Write the length stat to a file
 */
void fileStat(const string &lendisfile, const map<int,int> &lstat);
void correctQV(const string &bamf, const string &fasqf, const string &statfile, 
      const int lengthCutoff, const int shift);
void convert(const string &bamf, const string &fasqf, const string &statfile, const int lengthCutoff);

void usage() {
   cerr << "bam2cleanfastq inputfile.bam outputfile.fastq\n"
      << "  Or you can use command line options:\n"
      << "  -i inputfile.bam \n"
      << "  -o outputfile.fastq if not given this program will make one for you.\n"
      << "  -c lengthcutoff default 50\n"
      << "  -l length distribution file. All lengths seen in the input file.\n"
      << " There are two ways to specify the input and output file\n";
}

int main(int argc, char* argv[]) {
   // read command line arguments
   string inputfile, outputfile, statfile;
   //statfile = "length.dis";
   int lengthCutoff = 50;
   if (argc == 1) {
      usage();
      return 1;
   }
   int i = 1;
   bool correctQ=false;
   while (i < argc) {
      if (!strcmp(argv[i], "-i")) { inputfile = string(argv[++i]); }
      else if (!strcmp(argv[i], "-o")) { outputfile = string(argv[++i]); }
      else if (!strcmp(argv[i], "-c")) { lengthCutoff = atoi(argv[++i]); }
      else if (!strcmp(argv[i], "-l")) { statfile = argv[++i]; }
      else if (!strcmp(argv[i], "--correct")) { correctQ=true; }
      else if (!strcmp(argv[i], "--help")) { 
         usage();
         return 1;
      }
      else {
         inputfile = string(argv[i]);
         if (i+1 < argc && argv[i+1][0] != '-') {
            outputfile = string(argv[++i]);
         }
      }
      ++i;
   }

   if (inputfile.empty()) {
      usage();
      return 1;
   }

   if (outputfile.empty()) {
      outputfile = nameOutfile(inputfile);
   }
   if (statfile.empty()) {
      statfile = nameStatfile(inputfile);
   }
   if (correctQ) correctQV(inputfile, outputfile, statfile, lengthCutoff, -33);
   else convert(inputfile, outputfile, statfile, lengthCutoff);

   return 0;
}

void convert(const string &bamf, const string &fasqf, const string &statfile, const int lengthCutoff) {
   BamTools::BamAlignment bam;
   BamTools::BamReader breader;
   if (!breader.Open(bamf)) {
      cerr << "Failed to open " << bamf << endl;
      exit(1);
   }
   ofstream ouf(fasqf.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << fasqf << " for writing!\n";
      exit(1);
   }
   int shortCnt=0; 
   int simpleCnt=0;
   int goodCnt=0;
   int every=10;
   int i=1;
   map<int, int> lengthStat;
   map<string, int> ids;
   while (breader.GetNextAlignment(bam)) {
      if (i % every == 0) {
         cout << "working on " << i << endl;
         every *= 2;
      }
      if (bam.Length < lengthCutoff) {
         ++shortCnt;
         ++lengthStat[bam.Length];
      }
      else if (isArtifact(bam.QueryBases)) {
         ++simpleCnt;
      }
      else { // write to fastq
         ++ids[bam.Name];
         if (ids[bam.Name] > 1) {
            cerr << "duplicated sequence id: " << bam.Name << " discarded\n";
         }
         else {
            Fastq fasq(bam.Name, bam.QueryBases, bam.Qualities);
            ouf << fasq;
            ++goodCnt;
            ++lengthStat[bam.Length];
         }
      }
      ++i;
   }

   cout << shortCnt << " short sequences\n" << simpleCnt << " artifacts\n"
      << goodCnt << " good\n";
   fileStat(statfile, lengthStat);
}

// only do simple length cutoff
// correction ty 33
// checks for duplicated names
void correctQV(const string &bamf, const string &fasqf, const string &statfile, 
      const int lengthCutoff, const int shift) {
   BamTools::BamAlignment bam;
   BamTools::BamReader breader;
   if (!breader.Open(bamf)) {
      cerr << "Failed to open " << bamf << endl;
      exit(1);
   }
   ofstream ouf(fasqf.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open output fastq file " << fasqf << " for writing!\n";
      exit(1);
   }

   int shortCnt=0; 
   int simpleCnt=0;
   int goodCnt=0;
   int every=10;
   int i=1;
   map<int, int> lengthStat;
   map<string, int> ids;
   while (breader.GetNextAlignment(bam)) {
      if (i % every == 0) {
         cout << "working on " << i << endl;
         every *= 2;
      }
      if (bam.Length < lengthCutoff) {
         ++shortCnt;
         ++lengthStat[bam.Length];
      }
      else if (isArtifact(bam.QueryBases)) {
         ++simpleCnt;
      }
      else { // write to fastq
         ++ids[bam.Name];
         if (ids[bam.Name] > 1) {
            cerr << "duplicated sequence id: " << bam.Name << " discarded\n";
         }
         else {
            Fastq fasq(bam.Name, bam.QueryBases, bam.Qualities);
            fasq.shiftQuality(shift);
            ouf << fasq;
            ++goodCnt;
            ++lengthStat[bam.Length];
         }
      }
      ++i;
   }
   cerr << "min score: " << Fastq::getMinScore << " max score: " << Fastq::getMaxScore() << endl;
   cout << shortCnt << " short sequences\n" << simpleCnt << " artifacts\n"
      << goodCnt << " good\n";
   fileStat(statfile, lengthStat);
}

void printLengthStat(ostream &ous, const map<int,int> &lstat) {
   ous << "length\tcount\n";
   for (auto itr=lstat.begin(); itr != lstat.end(); ++itr) {
      ous << itr->first << "\t" << itr->second << endl;
   }
}

void fileStat(const string &lendisfile, const map<int,int> &lstat) {
   ofstream outf(lendisfile.c_str());
   if (outf.fail()) {
      cerr << "Failed to open " << lendisfile << " for writing!\n";
   }
   //outf << "Length satistic of sequences (excluding artifacts):\n";
   printLengthStat(outf, lstat);
   cerr << "Length satistic of sequences (excluding artifacts) written to file: "
      << lendisfile << endl;
   outf.close();
}

// if repeat of 10 or more is greater than 50%
// then this sequence is considered as junk.
bool isArtifact(const string &seq) {
   unsigned int i = 0;
   unsigned int lencut = 10;
   unsigned int repeatSum = 0;
   while (i < seq.length()) {
      unsigned int j = i+1;
      while (j<seq.length() && seq[j] == seq[i]) ++j;
      if (j-i > lencut) repeatSum += (j-i);
      i = j+1;
   }
   if (repeatSum / (float)seq.length() > 0.5) {
      cout << seq << " is artifact!\n";
      return true;
   }
   return false;
}


bool isSimple(const bioseq &bs) {
   map<char, double> freq = bs.getFrequency();
   if (freq.size() <= 2) return true;
   map<char, double>::const_iterator it = freq.begin();
   while (it != freq.end()) {
      if (it->second > 0.85) return true;
      ++it;
   }
   return false;
}
  
string nameOutfile(const string &infile) {
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      return infile.substr(0, i) + ".fastq";
   }
   return infile + ".fastq";
}

string nameStatfile(const string &infile) {
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      return infile.substr(0, i) + "_length.dis";
   }
   return infile + "_length.dis";
}

