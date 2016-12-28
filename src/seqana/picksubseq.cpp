#include <iostream>
#include <fstream>
#include <string>

#include <dynalnt.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;


void alnreftoall(const DNA &ref, const string &libfile, const string &outfile) {
   ifstream ifs(libfile.c_str());
   if (ifs.fail()) {
      exit(1);
   }
   ofstream ofs(outfile.c_str());
   if (ofs.fail()) {
      exit(1);
   }

   SimpleScoreMethod simplesm(10, -9, -19, -3);
   Dynaln<SimpleScoreMethod> align(simplesm);
   align.setSeq1(ref);

   DNA libseq;
   string header;
   int cntambiguous=0;
   int cnt=0;
   while (libseq.read(ifs, header)) {
      align.setSeq2(libseq);
      align.runlocal();
      //align.printAlign(cout);
      DNA sub = libseq.subseq(align.bottomBeginIndex()+1, align.bottomEndIndex()+1);
      //cout << sub << endl;
      if (sub.ambiguous()) {
         cout << "sequence has ambiguous base, ignored\n";
         ++cntambiguous;
      }
      else {
         ++cnt;
         align.printAlign(cout);
         DNA sub = libseq.subseq(align.bottomBeginIndex()+1, align.bottomEndIndex()+1);
         sub.appendTitle("identity to ref: " + to_string(align.getIdentity()));
         ofs << sub;
      }
   }
   cout << cnt << " noambiguous " << cntambiguous << " ambiguous\n";
}

int main(int argc, char* argv[]) {
   string refseqfile("hcvref.fas"), libseqfile("hcvblasthits.fas");
   string libsubfile("hcvamplicon.lib.fas");
   DNA refseq;
   refseq.read(refseqfile);
   alnreftoall(refseq, libseqfile, libsubfile);
}


