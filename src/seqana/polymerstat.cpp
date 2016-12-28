#include <iostream>
#include <string>
#include <fstream>
#include <map>

#include <bioseq.h>

using namespace std;
using namespace orpara;


void findRepeat(const bioseq& seq, map<char, map<int,int> > &res);
int main(int argc, char* argv[]) {
   //string fasfile="/path/to/any/fasta/file.fasta";
   string fasfile;
   int i = 1;
   while (i<argc) {
      fasfile=string(argv[i]);
      ++i;
   }

   ifstream inf(fasfile.c_str());
   if (!inf) {
      cerr << "Failed to open " << fasfile << " for reading\n";
   }
   bioseq fseq;
   string header;
   cerr << "processing " << fasfile << endl;
   i=0;
   map<char, map<int,int> > repeatCount;
   int every=100;
   size_t numbase = 0;
   while (fseq.read(inf, header)) {
      if (i % every == 0) {
         cerr << "working on " << i+1 << " " << fseq.getName() << endl;
         every *= 2;
      }
      findRepeat(fseq, repeatCount);
      numbase += fseq.length();
      ++i;
   }
   inf.close();
   cout << i << " sequences processed\n";
   cout << "repeat information " << endl;
   ofstream ouf("refseq16Snr_repeatinfo.txt");
   ouf << "derived from " << fasfile 
      << " total " << numbase << " bases" << endl;
   for (map<char, map<int,int> >::const_iterator b=repeatCount.begin(); b != repeatCount.end(); ++b) {
      ouf << b->first <<  endl;
      cout << b->first <<  endl;
      for (auto c=b->second.begin(); c != b->second.end(); ++c) {
         cout << c->first << " " << c->second 
            << " " << double(c->second)*c->first/numbase << endl;
         ouf << c->first << " " << c->second 
            << " " << double(c->second)*c->first/numbase << endl;
      }
      ouf << endl;
      cout << endl;
   }
   ouf.close();
   return 0;
}

void findRepeat(const bioseq& seq, map<char, map<int, int> > &res) {
   size_t i=0;
   while (i<seq.length()) {
      char r = seq[i];
      ++i;
      int rlen=1;
      while (r == seq[i]) {
         ++rlen;
         ++i;
      }
      ++res[r][rlen];
   }
}

