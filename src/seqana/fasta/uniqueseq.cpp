#include <iostream>
#include <string>
#include <fstream>
#include <map>
#include <set>
#include <cstring>

#include <bioseq.h>

using namespace std;
using namespace orpara;

// given a fasta file, it will eliminate identical sequences

void usage() {
   cerr << "uniqueseq seqfile.fas\n"
      << " you can give a sequence in fasta format, it will\n"
      << " eliminate all identical sequences from the input file\n"
      << " Then it will write into a new file with extension .uniq.fas\n"
      << "Optons\n"
      << "   -o outputFile the unique sequences\n"
      << "   -i inputFile the input fasta file that may have dupicated sequences\n"
      << "   -d fileName headers for identical sequences\n"
      << "   -dupid flag to produce duplicated ids in tabular format\n";
   exit(1);
}

/**
 * suff must contain . if you want it
 * input.abc after thisfunction("input.abc", ".def")
 * @return stem + suff, for example input.def
 */
string switchSuffix(const string &file, const string &suff) {
   string outfile = file;
   if (file.rfind('.') != string::npos) {
      outfile = file.substr(0, file.rfind('.'));
   }
   outfile += suff;
   return outfile;
}

void duplicatedIds(const string &input);

int main(int argc, char* argv[]) {
   string infasfile, outfile, dupheaderfile;
   bool listDupid = false;
   int i = 1;
   while (i < argc) {
      if (string(argv[i]) == "--help" || string(argv[i]) == "-h") {
         usage();
      }
      else if (!strcmp(argv[i], "-i")) {
         infasfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-d")) {
         dupheaderfile = argv[++i];
      }
      else if (!strcmp(argv[i], "--dupid")) {
         listDupid = true;
      }
      else {
         infasfile=string(argv[i]);
      }
      ++i;
   }
   if (infasfile.empty()) usage();

   if (outfile.empty()) {
      outfile = switchSuffix(infasfile, ".uniq.fas");
   }
   ofstream ofs(outfile);
   if (ofs.fail()) {
      throw runtime_error("Failed to open outfile " + outfile);
   }
   if (dupheaderfile.empty()) {
      dupheaderfile = switchSuffix(infasfile, "_duph.txt");
   }
   ofstream odu(dupheaderfile);
   if (odu.fail()) throw runtime_error("cannot open " + dupheaderfile);
   odu << "Only the first sequence id used to collect duplicated headers\n";

   set<bioseq> uniqueseq;
   bioseq bs;
   ifstream ifs(infasfile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open infasfile: " << infasfile << endl;
      return 1;
   }
   int cnt=0;
   // store duplicated headers
   map<string, vector<string> > defline;
   while (bs.read(ifs)) {
      auto itr = uniqueseq.find(bs);
      if (itr != uniqueseq.end()) {
         auto mit = defline.find(itr->getName());
         if (mit == defline.end()) {
            // the tile in the set will not be duplicated here
            vector<string> hh = { bs.getName() + ": " + bs.getTitle() };
            defline.insert(make_pair(itr->getName(), hh));
         }
         else {
            mit->second.push_back(bs.getName() + ": " + bs.getTitle());
         }
      }
      else {
         uniqueseq.insert(bs);
      }
      ++cnt;
   }
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
   cout << uniqueseq.size() << " unique sequences out of " << cnt 
      << " written to " << outfile << "\n"
      << " duplicated header written to " << dupheaderfile << endl;

   if (listDupid) {
      duplicatedIds(infasfile);
   }
   return 0;
}

void duplicatedIds(const string &input) {
   set<bioseq> uniqueseq;
   bioseq bs;
   ifstream ifs(input);
   if (ifs.fail()) {
      throw runtime_error("Failed to open infasfile: " + input);
   }
   string dupidfile = switchSuffix(input, "_dupid.tab");
   ofstream odu(dupidfile);
   if (odu.fail()) throw runtime_error("cannot open " + dupidfile);
   odu << "Only the first sequence id used to collect duplicated headers\n";
   int cnt=0;
   // store duplicated ids
   map<string, vector<string> > dupids;
   while (bs.read(ifs)) {
      auto itr = uniqueseq.find(bs);
      if (itr != uniqueseq.end()) {
         auto mit = dupids.find(itr->getName());
         if (mit == dupids.end()) {
            vector<string> dd = { bs.getName() };
            dupids.insert(make_pair(itr->getName(), dd));
         }
         else {
            mit->second.push_back(bs.getName());
         }
      }
      else {
         uniqueseq.insert(bs);
      }
      ++cnt;
   }
   map<string, vector<string> >::const_iterator it;
   for (it=dupids.begin(); it != dupids.end(); ++it) {
      for (unsigned i=0; i < (it->second).size(); ++i) {
         odu << it->first << '\t' << (it->second)[i] << endl;
      }
   }
   cerr << "duplicated ids written to " << dupidfile << endl;
}

