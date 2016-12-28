#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <cstring>

#include <bioseq.h>

using namespace std;
using namespace orpara;

// given two fasta files, it will eliminate identical sequences
// and combine the two. Always eliminate from the second file
// so there is a directionality

void usage() {
   cerr << "combineseq seqfile1.fas seqfile2.fas\n"
      << " you can give a sequence in fasta format, it will\n"
      << " eliminate all identical sequences from the input file\n"
      << " Then it will write into a new file with extension .uniq.fas\n"
      << "Optons\n"
      << "   -o outputFile the unique sequences\n"
      << "   -1 inputFile1 the input fasta file that may have dupicated sequences\n"
      << "   -2 inputFile2 the input fasta file that should not have dupicated sequences\n"
      << "   -d fileName headers for identical sequences\n";
   exit(1);
}

/**
 * suff must contain . if you want it
 * input.abc after thisfunction("input.abc", ".def")
 * @return stem + suff, for example input.def
 */
string makeOutputFile(const string &file1, const string &file2, const string &suff) {
   string outfile = file1;
   string::size_type i = outfile.rfind('/');
   if (i != string::npos) {
      outfile = outfile.substr(i+1);
   } // get rid of path info
   string tmp = file2;
   i = file2.rfind('/');
   if (i != string::npos) {
      tmp = file2.substr(i+1);
   }
   if (outfile.rfind('.') != string::npos) {
      outfile = outfile.substr(0, outfile.rfind('.'));
   }
   if (tmp.rfind('.') != string::npos) {
      outfile += "_" + tmp.substr(0, tmp.rfind('.'));
   }
   else outfile += "_" + tmp;
   outfile += suff;
   return outfile;
}
/**
 * @param input is the input sequence file name.
 *   Sequences in this file should be unique.
 *   If not unique, then duplicated entries will be
 *   discarded.
 */
set<bioseq> slurpFasta(const string &input);

int main(int argc, char* argv[]) {
   string infile1, infile2, outfile, dupheaderfile;
   int i = 1;
   while (i < argc) {
      if (string(argv[i]) == "--help" || string(argv[i]) == "-h") {
         usage();
      }
      else if (!strcmp(argv[i], "-1")) {
         infile1 = argv[++i];
      }
      else if (!strcmp(argv[i], "-2")) {
         infile2 = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-d")) {
         dupheaderfile = argv[++i];
      }
      else {
         if (i+1 > argc) {
            usage();
         }
         infile1=string(argv[i]);
         infile2=string(argv[++i]);
      }
      ++i;
   }
   if (infile1.empty() || infile2.empty()) usage();

   if (outfile.empty()) {
      outfile = makeOutputFile(infile1, infile2, ".comb.fas");
   }
   ofstream ofs(outfile);
   if (ofs.fail()) {
      throw runtime_error("Failed to open outfile " + outfile);
   }
   if (dupheaderfile.empty()) {
      dupheaderfile = makeOutputFile(infile1, infile2, "_duph.txt");
   }
   ofstream odu(dupheaderfile);
   if (odu.fail()) throw runtime_error("cannot open " + dupheaderfile);
   odu << "Only the first sequence id used to collect duplicated headers\n";

   set<bioseq> uniqueseq = slurpFasta(infile1);
   bioseq bs;
   ifstream ifs(infile2);
   if (ifs.fail()) {
      cerr << "Failed to open second fasta file: " << infile2 << endl;
      return 1;
   }
   int cnt=0;
   int dupcnt=0;
   map<string, vector<string> > defline;
   while (bs.read(ifs)) {
      auto itr = uniqueseq.find(bs);
      if (itr != uniqueseq.end()) {
         auto mit = defline.find(itr->getName());
         if (mit == defline.end()) {
            // the tile in the set will not be duplicated here
            vector<string> hh = { bs.getTitle() };
            defline.insert(make_pair(itr->getName(), hh));
         }
         else {
            mit->second.push_back(bs.getTitle());
         }
         ++dupcnt;
      }
      else {
         ofs << bs;
      }
      ++cnt;
   }
   // write sequences from the first file to output
   map<string, vector<string> >::const_iterator it;
   for (auto itr=uniqueseq.begin(); itr != uniqueseq.end(); ++itr) {
      it = defline.find(itr->getName());
      if (it != defline.end()) {
         odu << itr->getName() << endl << itr->getTitle() << "\nmultiple headers:\n";
         for (unsigned i=0; i < (it->second).size(); ++i) {
            odu << (it->second)[i] << endl;
         }
         odu << string(70, '-') << endl;
      }
      ofs << *itr;
   }
   cout << uniqueseq.size() << " unique sequences from " << infile1
      << endl
      << cnt << " from " << infile2 << " " << dupcnt << " duplicated\n" 
      << " combined sequences to " << outfile << "\n"
      << " duplicated header written to " << dupheaderfile << endl;
   return 0;
}

set<bioseq> slurpFasta(const string &input) {
   set<bioseq> uniq;
   bioseq bs;
   ifstream ifs(input);
   if (ifs.fail()) {
      throw runtime_error("Failed to open infasfile: " + input);
   }
   int cnt=0;
   while (bs.read(ifs)) {
      uniq.insert(bs);
      ++cnt;
   }
   cerr << cnt << " sequences from " << input << endl;
   return uniq;
}
