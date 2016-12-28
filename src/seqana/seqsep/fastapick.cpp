#include <iostream>
#include <fstream>
#include <list>
#include <stdexcept>
#include <functional>
#include <string>

#include "bioseq.h"
#include "strformat.h"

using namespace std;
using namespace orpara;

list<DNA> byWord(const string &file, const string &wd, const unsigned int length);

template<class C>
list<DNA> byGCContent(const string &file, const C &cmp, unsigned int lencut) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot open file: " + file);
   }
   list<DNA> result;
   while (dna.read(inf)) {
      if (dna.length() > lencut && cmp(dna.GCContent())) {
         result.push_back(dna);
      }
   }
   return result;
}
void filterByLength(const string &inputf, unsigned int length, const string &outputf);
list<DNA> byWordList(const string &file, const set<string> &wd, const unsigned int length);
/**
 * One phrase per line
 * Should be reading phrases.
 */
set<string> readWords(const string &file);
/**
 * @param lencut length cutoff for the input sequence.
 * @param wlfile is the file containing all words to be picked
 *     one or more words per line. The whole line is considered as
 *     one "word" in the search process.
 */
void pickAllWords(const string &inputf, string &outputf, const string &wlfile, 
      unsigned int lencut);
/**
 * fetch one DNA object by unique id
 * @return only one object
 */
DNA byId(const string &file, const string &id);
void pickByIdList(const string &inputf, string &outputf, const string &idfile,
      ios_base::openmode fmode);
set<string> readId(const string &infile);
/**
 * @param idlist is a string of ids separated by comma
 */
set<string> splitId(const string &idlist);
string makeListOutfileName(const string &ifname, const string &lfname);

void usage() {
   cerr << "fastapick --gc '>0.65' -i inputfasta -o outputfile\n"
      << "Options:\n"
      << "   -w word in the title of the sequence\n"
      << "   -W FILE containing phrases one per line to be used as query\n"
      << "      to search the title of fasta sequence\n"
      << "   --gc > or <FLOAT number. Pick sequence with GC content\n"
      << "     above or below certain value. [0-1].\n"
      << "   -l LENGTH_CUTOFF default 1400\n"
      << "   -o OUTPUT_FILE name, if not given will make one for you\n"
      << "   --id seqid, unique identifier for a sequence\n"
      << "   --idlist file_of_seqids, one or more id per line, one or more lines.\n"
      << "        or a list of ids separate by comma. Single or double quote\n"
      << "        provide a hint that it is a list of ids not file name\n"
      << "   --append flag to append to output file\n"
      << "   -h or --help print this message\n"
      << " Length cutoff can be combined with other operation with AND logic.\n"
      << " Note: the GC cutoff must be quoted, and shell use > \n"
      << "      for special purposes\n";
   exit(1);
}

void writeSequences(const list<DNA> &seq, const string &file);

int main(int argc, char* argv[]) {
   //string word="Mycobacterium";
   string word, id;
   string inputFile;
   float gc=0;
   bool gcbelow=true;
   unsigned int lengthcut=1400;
   string outputFile;
   string wordlist, idlist;
   // for testing
   /*
   idlist="Bac563_1114,Bac563_1115,Bac563_1116";
   inputFile="/isilon/Data/inf_dis/metag/mashref/genomes/bacteria.563.1.genomic_sis.fas";
   outputFile="xxyytest.fas";
   */

   ios_base::openmode writeMode=ios_base::out;
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-w")) { word=argv[++i]; }
      else if (!strcmp(argv[i], "-W")) { wordlist=argv[++i]; }
      else if (!strcmp(argv[i], "-i")) { inputFile=argv[++i]; }
      else if (!strcmp(argv[i], "-l")) { lengthcut=atoi(argv[++i]); }
      else if (!strcmp(argv[i], "-o")) { outputFile=argv[++i]; }
      else if (!strcmp(argv[i], "--id")) { id=argv[++i]; }
      else if (!strcmp(argv[i], "--idlist")) { idlist=argv[++i]; }
      else if (!strcmp(argv[i], "--append")) { writeMode=ios_base::app; }
      else if (!strcmp(argv[i], "--gc")) {
         if (argv[i+1][0] == '>') {
            gc=atof(argv[i+1]+1);
            gcbelow=false;
         }
         else if (argv[i+1][0] == '<') {
            gc=atof(argv[i+1]+1);
         }
         else {
            usage();
         }
         ++i;
      }
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")
            || !strcmp(argv[i], "?")) 
      {
         usage();
      }
      else if (argv[i][0] == '-') {
         cerr << "option " << argv[i] << " not recognized\n";
         usage(); return 1;
      }
      else {
         inputFile=argv[i];
      }
      ++i;
   }
   if (inputFile.empty()) {
      cerr << "you must provide an inputfile\n";
      usage();
      return 1;
   }
   if (!id.empty()) {
      DNA result = byId(inputFile, id);
      if (result.empty()) {
         cerr << "did not find sequence with id: " 
            << id << " in file: " << inputFile << endl;
         return 1;
      }
      if (outputFile.empty()) {
         outputFile=inputFile;
         string::size_type x = outputFile.rfind('.');
         if (x != string::npos) {
            outputFile = outputFile.substr(0,x);
         }
         outputFile += id + ".fas";
      }
      ofstream ouf(outputFile);
      if (ouf.fail()) {
         cerr << "Failed to open " << outputFile << " for writing\n";
         return 1;
      }
      ouf << result << endl;
      cerr << "result written to " << outputFile << endl;
      return 0;
   }
   if (!idlist.empty()) {
      pickByIdList(inputFile, outputFile, idlist, writeMode);
      return 0;
   }

   if (!word.empty()) {
      outputFile=inputFile;
      string::size_type x = inputFile.rfind('/');
      if (x != string::npos) {
         outputFile = inputFile.substr(x+1);
      }
      string suffix=word;
      for (size_t s=0; s<suffix.size(); ++s) {
         if (suffix[s] == ' ')
            suffix.replace(s, 1, 1, '_');
      }
      outputFile += "." + suffix;
      list<DNA> myselect = byWord(inputFile, word, lengthcut);
      cout << myselect.size() << " " << word << " sequences\n";
      writeSequences(myselect, outputFile);
   }

   list<DNA> GCextreme;
   string gcfile;
   if (gc > 0) {
      cout << "GC cutoff: " << gc << endl;
      if (gcbelow) {
         GCextreme=byGCContent(inputFile, bind2nd(less<float>(), gc), lengthcut);
         cerr << GCextreme.size() << " gc < " << gc << endl;
         gcfile=inputFile + ".GCLow";
      }
      else {
         GCextreme=byGCContent(inputFile, bind2nd(greater<float>(), gc), lengthcut);
         cerr << GCextreme.size() << " gc > " << gc << endl;
         gcfile=inputFile + ".GCHight";
      }
      writeSequences(GCextreme, gcfile); 
      cerr << "extreme GC sequences written to " << gcfile << endl;
   }
   if (!wordlist.empty()) {
      pickAllWords(inputFile, outputFile, wordlist, lengthcut);
   }
   if (word.empty() && gc < 0.01 && wordlist.empty()) {
      // do length filtering only
      if (outputFile.empty()) {
         outputFile=inputFile;
         string::size_type x=outputFile.rfind('.');
         if (x != string::npos) {
            outputFile = outputFile.substr(0,x);
         }
         outputFile += ".lencut.fas";
      }
      filterByLength(inputFile, lengthcut, outputFile);
   }

   return 0;
}

set<string> readWords(const string &file) {
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("failed to open " + file);
   }
   set<string> words;
   string wd;
   getline(inf,wd);
   while (!inf.eof()) {
      words.insert(getLower(wd));
      getline(inf,wd);
   }
   cerr << words.size() << " words\n";
   return words;
}

void filterByLength(const string &inputf, unsigned int length, const string &outputf) {
   ifstream inf(inputf);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + inputf);
   }
   ofstream ouf(outputf);
   if (ouf.fail()) {
      throw runtime_error("cannot write to file: " + outputf);
   }
   DNA dna;
   int cnt=0;
   int p=0;
   int every=10;
   while (dna.read(inf)) {
      ++cnt;
      if (dna.length() > length) {
         ouf << dna;
         ++p;
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
   cout << p << " picked out of " << cnt << " sequences\n";
}

void writeSequences(const list<DNA> &seq, const string &file) {
   ofstream ouf(file);
   if (ouf.fail()) {
      throw runtime_error("failed to open file: " + file);
   }
   for (auto it=seq.begin(); it != seq.end(); ++it) {
      ouf << *it;
   }
   cout << seq.size() << " sequences written to " << file << endl;
}

list<DNA> byWord(const string &file, const string &wd, const unsigned int length) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + file);
   }
   cerr << "input file: " << file << endl;
   list<DNA> result;
   string lowerwd=getLower(wd);
   int cnt=0;
   int every=10;
   while (dna.read(inf)) {
      ++cnt;
      if (dna.length() > length && getLower(dna.getTitle()).find(lowerwd) != string::npos) {
         result.push_back(dna);
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
   return result;
}

// by single id
DNA byId(const string &file, const string &id) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + file);
   }
   cerr << "input file: " << file << endl;
   int every=10;
   int cnt=0;
   while (dna.read(inf)) {
      ++cnt;
      if (dna.getName() == id) {
         return dna;
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
   dna.clear();
   return dna;
}

bool titleHasWord(const string &title, const set<string> &wds) {
   for (auto it=wds.begin(); it != wds.end(); ++it) {
      if (title.find(*it) != string::npos) return true;
   }
   return false;
}

list<DNA> byWordList(const string &file, const set<string> &wd, const unsigned int length) {
   DNA dna;
   ifstream inf(file);
   if (inf.fail()) {
      throw runtime_error("cannot read file: " + file);
   }
   cerr << "input file: " << file << endl;
   list<DNA> result;
   int cnt=0;
   int every=10;
   while (dna.read(inf)) {
      ++cnt;
      string title=getLower(dna.getTitle());
      if (dna.length() > length && titleHasWord(getLower(dna.getTitle()), wd)) {
         result.push_back(dna);
      }
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " ...\n";
         every *= 2;
      }
   }
   return result;
}

void pickAllWords(const string &inputf, string &outputf, const string &wlfile, 
      unsigned int lencut) {
   if (outputf.empty()) {
      if (inputf.rfind('.') != string::npos) {
         outputf=inputf.substr(0, inputf.rfind('.'));
      }
      else outputf = inputf;
      outputf += ".bywdlist.fas";
   }
   set<string> wlist=readWords(wlfile);
   list<DNA> candidate = byWordList(inputf, wlist, lencut);
   writeSequences(candidate, outputf);
}

set<string> readId(const string &infile) {
   ifstream inf(infile);
   if (inf.fail()) throw runtime_error("failed to open " + infile);
   string id;
   set<string> res;
   inf >> id;
   while (!inf.eof()) {
      res.insert(id);
      inf >> id;
   }
   cerr << res.size() << " ids read from " << infile << endl;
   cerr << "last one " << id << endl;
   return res;
}

set<string> splitId(const string &idlist) {
   string tmp=idlist;
   if (idlist[0] == '\'' || idlist[0] == '"') {
      tmp = idlist.substr(1);
      if (tmp[tmp.size()-1] == '\'' || tmp[tmp.size()-1] == '"')
         tmp=tmp.substr(0,tmp.size()-1);
   }
   string::size_type i,j;
   i=0;
   set<string> res;
   while (i<tmp.size()) {
      j = tmp.find(',', i);
      if (j == string::npos) {
         res.insert(tmp.substr(i));
         break;
      }
      else {
         res.insert(tmp.substr(i, j-i));
         i=j+1;
      }
   }
   /*
   for (auto it=res.begin(); it != res.end(); ++it) {
      cerr << *it << endl;
   }
   cerr << res.size() << " ids from input string\n";
   */
   return res;
}

string makeListOutfileName(const string &ifname, const string &lfname) {
   // make up a convenient output file name for the user
   string outputf = ifname;
   string::size_type i = outputf.rfind('.');
   if (i != string::npos) {
      outputf = outputf.substr(0,i);
   }
   i = lfname.rfind('.');
   outputf += "_byid_";
   if (i == string::npos) {
      outputf += lfname + ".fas";
   }
   else {
      outputf += lfname.substr(0,i) + ".fas";
   }
   return outputf;
}

void pickByIdList(const string &inputf, string &outputf, const string &idfile,
      ios_base::openmode fmode) {
   cerr << "picking up sequences from " << inputf << endl;
   set<string> idlist;
   if (idfile[0] == '\'' || idfile[0] == '"' || idfile.find(',') != string::npos) {
      cerr << "list of ids provieded: " << idfile << endl;
      // id provided directly not through file
      idlist = splitId(idfile);
   }
   else {
      cerr << "idlist file provided: " << idfile << endl;
      if (outputf.empty()) {
         outputf=makeListOutfileName(inputf, idfile);
      }
      idlist=readId(idfile);
   }
   set<string>::const_iterator it;
   ifstream inf(inputf);
   if (inf.fail()) {
      throw runtime_error("failed to open " + inputf);
   }
   ofstream ouf(outputf, fmode);
   DNA dna;
   unsigned int cnt=0;
   cerr << idlist.size() << " ids to fetch\n";
   while (dna.read(inf)) {
      it = idlist.find(dna.getName());
      if (it != idlist.end()) {
         //cerr << cnt << " found " << dna.getName() << endl;
         ouf << dna;
         ++cnt;
         if (cnt == idlist.size()) {
            cerr << "all " << cnt << " candidates found\n";
            break;
         }
      }
   }
   cerr << cnt << " sequences written to " << outputf << endl;
}

