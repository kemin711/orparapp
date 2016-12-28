#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <string.h>

#include <fastq.h>
#include <stddev.h>

using namespace std;
using namespace orpara;


/**
 * collect length and average q value
 */
map<string, pair<int, double> > loadSeqinfo(const string &fname);
void joinInfo(const map<string, pair<int,double> > &oldq, 
      const map<string, pair<int,double> > &newq, const string &fname);
string shortenMovieName(const string &mvn);
void usage() {
   cerr << "usage comparefastq -n newfile -o oldfile\n"
      << "  or comparefastq newfile oldfile\n";
}
int main(int argc, char* argv[]) {
   //string newfile="new/A_onepass_ccs.fastq";
   //string oldfile="old/m150911_011314_42232_c100853722550000001823189501241694_s1_p0.ccs.fastq";
   string newfile, oldfile;

   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-n")) {
         newfile=argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         oldfile=argv[++i];
      }
      else if (i+1<argc) {
         newfile=argv[i];
         oldfile=argv[++i];
      }
      ++i;
   }
   if (newfile.empty()) {
      usage();
      return 1;
   }

   map<string, pair<int,double> > newfq = loadSeqinfo(newfile);
   map<string, pair<int,double> > oldfq = loadSeqinfo(oldfile);
   cerr << "fastq length and average QV loaded from two files\n";
   joinInfo(oldfq, newfq, "lenqvcompare.tab");

   return 0;
}

string shortenMovieName(const string &mvn) {
   size_t i = mvn.find('/');
   size_t j = mvn.rfind('/');
   return "ccs" + mvn.substr(i+1, j-i-1);
}

void joinInfo(const map<string, pair<int,double> > &oldq, 
      const map<string, pair<int,double> > &newq, const string &fname) {
   ofstream ouf(fname.c_str());
   // old id old_seqlen old_QV new_len new_QV
   ouf << "seqid\told_len\told_qv\tnew_len\tnew_qv\n";
   map<string, pair<int,double> >::const_iterator it, newit;
   it=oldq.begin();
   int missed=0;
   while (it != oldq.end()) {
      ouf << shortenMovieName(it->first) << "\t" << it->second.first << "\t" << it->second.second;
      newit = newq.find(it->first);
      if (newit == newq.end()) { // not found in new
         ouf << "\t0\t0\n";
         ++missed;
      }
      else {
         ouf << "\t" << newit->second.first << "\t" << newit->second.second << endl;
      }
      ++it;
   }
   cerr << missed << " old fastq not found in new fastq\n";
   // in new but not found in old
   newit = newq.begin();
   missed=0;
   while (newit != newq.end()) {
      it = oldq.find(newit->first);
      if (it == oldq.end()) {
         ouf << shortenMovieName(newit->first) << "\t0\t0\t" << newit->second.first 
            << "\t" << newit->second.second << endl;
         ++missed;
      }
      ++newit;
   }
   cerr << missed << " new fastq not found in old fastq\n";
   ouf.close();
   cerr << "result written to " << fname << endl;
}

map<string, pair<int, double> > loadSeqinfo(const string &fname) {
   ifstream inf(fname.c_str());
   if (!inf) {
      cerr << "Failed to open " << fname << endl;
      exit(1);
   }
   Fastq fastq;
   map<string, pair<int, double> > lenq;
   while (fastq.read(inf)) {
      stddev stat;
      int *qv = fastq.getQuality();
      for (size_t i=0; i<fastq.length(); ++i) {
         stat(qv[i]);
      }
      lenq[fastq.getName()]=make_pair(fastq.length(), stat.getMean() - Fastq::getConverter());
   }
   cerr << lenq.size() << " sequences processed\n";
   return lenq;
}
