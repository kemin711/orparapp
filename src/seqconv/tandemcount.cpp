#include <iostream>
#include <fstream>

#include <bioseq.h>

using namespace std;
using namespace orpara;

class TandemBase {
   public:
      TandemBase() : base('?'), count(0) { }
      TandemBase(char b) : base(b), count(1) { }
      TandemBase(char b, int c) : base(b), count(c) { }
      string toString() const { return string(count, base); }
      friend ostream& operator<<(ostream &ous, const TandemBase &tb) {
         ous << (char)toupper(tb.base) << '\t' << tb.count; return ous; }
      // base first the count
      /*
      bool operator<(const TandemBase &tb) const {
         if (base < tb.base) return true;
         else return count < tb.count; }
      */
      // count first then base
      bool operator<(const TandemBase &tb) const {
         if (count < tb.count) return true;
         else if (count > tb.count) return false;
         else return base < tb.base; 
      }

   private:
      char base;
      int count;
};

string makeOutputFileName(const string &inf);
void countRepeats(const string &seq, map<TandemBase, int> &tc);
void writeStore(const map<TandemBase, int> &tc, ostream &ous);

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

   DNA dna;
   dna.read(infile);
   map<TandemBase, int> tstore;
   countRepeats(dna.toString(), tstore);
   if (outfile.empty()) {
      outfile = makeOutputFileName(infile);
   }
   ofstream of(outfile.c_str());
   writeStore(tstore, of);
   writeStore(tstore, cout);
   cout << "Tandem repeats written to file: " << outfile << endl;

   return 0;
}

void writeStore(const map<TandemBase, int> &tc, ostream &ous) {
   map<TandemBase, int>::const_iterator it = tc.begin();
   while (it != tc.end()) {
      ous << it->first << "\t" << it->second << endl;
      ++it;
   }
}

void countRepeats(const string &seq, map<TandemBase, int> &tc) {
   unsigned int i = 0;
   while (i < seq.length()) {
      char B = seq[i];
      int count = 1;
      unsigned int j = i + 1;
      if (j == seq.length()) break;
      while (j < seq.length() && seq[j] == B) {
         ++count;
         ++j;
      }
      if (count > 3) {
         ++tc[TandemBase(B, count)];
      }
      i = j + 1;
   }
}

string makeOutputFileName(const string &inf) {
   string::size_type i;
   string ofile=inf;
   if ((i=inf.rfind('.')) != string::npos) {
      ofile = inf.substr(0, i) + "tandem.list";
   }
   return ofile;
}
