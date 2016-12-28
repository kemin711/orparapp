#include <iostream>
#include <string>
#include <fstream>
#include <stddev.h>
#include <cstring>
#include <exception>
#include <glob.h>
#include <vector>

#include <bioseq.h>

using namespace std;
using namespace orpara;

// given a fasta file, it will eliminate identical sequences

void usage() {
   cerr << "fastaleninfo seqfile.fas\n"
      << " Generate length summary information for the input file.\n"
      << " --list fileName list the length of the fasta file\n";
   exit(1);
}

/**
 * produce seqid => seqlength table
 * from a fasta file.
 */
void listLength(const string &inf, const string &ouf);
/**
 * @return number of sequences processed from inf
 */
int listLength(const string &inf, ostream &ous);
void listLengthMultiple(const vector<string> &infiles, const string &ouf);
vector<string> getInputFiles(const string &pat);

int main(int argc, char* argv[]) {
   string infasfile, outfile, pattern;
   int i = 1;
   while (i < argc) {
      if (string(argv[i]) == "--help") {
         usage();
      }
      else if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      else if (!strcmp(argv[i], "--list")) {
         pattern=argv[++i];
      }
      else {
         infasfile=string(argv[i]);
      }
      ++i;
   }

   if (infasfile.empty() && pattern.empty()) {
      cerr << "You must provided either a input file or a patern through --list option\n";
      usage();
   }

   // work on multiple input files
   if (!pattern.empty()) {
      vector<string> files = getInputFiles(pattern);
      if (outfile.empty()) {
         outfile = "sequence_length.tab";
      }
      listLengthMultiple(files, outfile);
      return 0;
   }

   if (outfile.empty()) {
      outfile = infasfile.substr(0, infasfile.rfind('.')) + ".lenstat";
   }
   ofstream ofs(outfile.c_str());
   if (ofs.fail()) {
      cerr << "Failed to open outfile " << outfile << endl;
      return 1;
   }

   stddev as; // average and standard deviation
   bioseq bs;
   ifstream ifs(infasfile);
   if (ifs.fail()) {
      cerr << "Failed to open infasfile: " << infasfile << endl;
      return 1;
   }
   while (bs.read(ifs)) {
      as(bs.length());
   }
   ofs << as << endl;
   cout << as << endl;
   cout << "sequence length summary information written to "
      << outfile << endl;

   return 0;
}

vector<string> getInputFiles(const string &pat) {
   cerr << "get files matching " << pat << endl;
   glob_t fobj;
   glob(pat.c_str(), 0, NULL, &fobj);
   cerr << "Found " << fobj.gl_pathc << " files from pattern: " << pat << endl;
   vector<string> result;
   for(unsigned int i=0; i < fobj.gl_pathc; ++i){
      result.push_back(string(fobj.gl_pathv[i]));
   }
   globfree(&fobj);
   return result;
}

void listLength(const string &inf, const string &ouf) {
   ofstream ofs(ouf);
   listLength(inf, ofs);
}

int listLength(const string &inf, ostream &ous) {
   bioseq bs;
   ifstream ifs(inf);
   if (ifs.fail()) {
      cerr << "Failed to open infasfile: " << inf << endl;
      throw runtime_error("cannot open file: " + inf);
   }
   cout << "processing " << inf << endl;
   int cnt=0;
   while (bs.read(ifs)) {
      ous << bs.getName() << '\t' << bs.length() << endl;
      ++cnt;
   }
   cout << inf << " has " << cnt << " sequences\n";
   return cnt;
}

void listLengthMultiple(const vector<string> &infiles, const string &ouf) {
   ofstream ofs(ouf);
   int totalseq=0;
   for (int i=0; i<infiles.size(); ++i) {
      totalseq += listLength(infiles[i], ofs);
   }
   cout << totalseq << " total sequences processed\n";
}

