#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <cstdlib>

using namespace std;


void cleanfasta(const string &file, const string &ofile, const string &prefix);
void trimN(string &seq);
// use composition to see sequence is trash or not
bool isSimple(const string &seq);

int main(int argc, char* argv[]) {
   int i=1;
   string infile, oufile;
   string prefix;
   while (i<argc) {
      if (!strcmp(argv[i], "-o")) {
         oufile=argv[++i];
      }
      else if (!strcmp(argv[i], "-p")) {
         prefix=argv[++i];
      }
      else {
         infile=argv[i];
      }
      ++i;
   }
   if (oufile.empty()) oufile=infile + ".clean";
   cleanfasta(infile,oufile, prefix);
   return 0;
}

void cleanfasta(const string &file, const string &ofile, const string &prefix) {
   ifstream ins(file.c_str());
   ofstream ous(ofile.c_str());
   if (ins.fail()) {
      cerr << "Failed to open " << file << endl;
      exit(1);
   }
   string line;
   getline(ins, line);
   int unsigned id=1;
   while (!ins.eof() && line[0] == '>') {
      getline(ins,line);
      string seq;
      while (!ins.eof() && !line.empty() && line[0] != '>') {
         seq += line;
         getline(ins, line);
      }
      trimN(seq);
      if (seq.length() > 24 && !isSimple(seq)) {
         ous << '>';
         if (!prefix.empty()) ous << prefix;
         ous << id++ << endl << seq << endl;
      }
   }
}

void trimN(string &seq) {
   string::size_type i=0;
   while (i<seq.length() && seq[i] == 'N') {
      ++i;
   }
   if (i == seq.length() || seq.length()-i < 26) {
      seq.clear();
      return;
   }
   seq=seq.substr(i);
   i=seq.length()-1;
   while (i>0 && seq[i] == 'N') --i;
   if (i < 25) { 
      seq.clear();
      return;
   }
   i=seq.find('N');
   if (i == string::npos) return;
   while (seq.length() > 25 && i != string::npos) {
      if (i > seq.length()/2) {
         seq.resize(i);
      }
      else {
         seq=seq.substr(i+1);
      }
      i=seq.find('N');
   }
   if (seq.length() < 26) seq.clear();
}

bool isSimple(const string &seq) {
   // A, C, G, T, N
   int B[6];
   for (int i=0; i<seq.length(); ++i) {
      if (seq[i] == 'A' || seq[i] == 'a') {
         ++B[0];
      }
      else if (seq[i] == 'C' || seq[i] == 'c') {
         ++B[1];
      }
      else if (seq[i] == 'G' || seq[i] == 'g') {
         ++B[2];
      }
      else if (seq[i] == 'T' || seq[i] == 't') {
         ++B[3];
      }
      else if (seq[i] == 'N' || seq[i] == 'n') {
         ++B[4];
      }
      else {
         cerr << "new base: " << seq[i] << " not N or 4 canonical bases\n";
      }
   }
   int sum=B[0]+B[1]+B[2]+B[3]+B[4];
   static const float cutoff=0.9;
   if (float(B[0])/sum > cutoff || float(B[1])/sum > cutoff
         || float(B[2])/sum > cutoff || float(B[3])/sum > cutoff 
         || float(B[4])/sum > cutoff)  {
      cerr << seq << " is simple\n";
      return true;
   }
   return false;
}
