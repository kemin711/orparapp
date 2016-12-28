#include <iostream>
#include <string>
#include <fstream>
#include <set>
#include <functional>

#include <strformat.h>
#include <dynalnt.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;


// given a fasta file, it will eliminate identical sequences

string trimVersion(const string &seqid) {
   size_t x = seqid.find('.');
   if (x != string::npos) return seqid.substr(0, x);
   else return seqid;
}

class PidNameCleanup : public unary_function<string,string> {
   public:
      string operator()(const string& oldname) {
         if (oldname.substr(0, 9) == "UniRef50_") {
            return oldname.substr(9);
         }
         //else if (oldname.substr(0,8) == "BMSPROT:") {
         else if (oldname.substr(0,8) == "GCGPROT:") {
            size_t dot = oldname.rfind('.');
            if (dot != string::npos) {
               return oldname.substr(8, dot-8);
            }
            else {
               cerr <<  "protein id has no version: " << oldname << endl;
               exit(1);
            }
         }
         else if (oldname.substr(0,3) == "gi|") {
            vector<string> elem = split(oldname, '|');
            return trimVersion(elem[3]);
         }
         else {
            cerr << "unexpected old pid name: " << oldname << endl;
            exit(1);
         }
      }
};

class RefseqNameCleanup : public unary_function<string,string> {
   public:
      string operator()(const string& oldname) {
         if (oldname.substr(0, 8) == "REFSEQP:") {
            size_t dot = oldname.rfind('.');
            if (dot == string::npos) {
               cerr << oldname << " has no dot at end!\n";
               exit(1);
            }
            else {
               return oldname.substr(8, dot-8);
            }
         }
         else if (oldname.substr(0,3) == "gi|") {
            vector<string> elem = split(oldname, '|');
            return trimVersion(elem[3]);
         }
         else {
            cerr << "unexpected old refseq name: " << oldname << endl;
            exit(1);
         }
      }
};

// read fasta 
map<string, Protein> collectRefseqpFasta(const string &fastaFile, map<string,int> &list) {
   ifstream ifs(fastaFile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open " << fastaFile << endl;
      exit(1);
   }
   map<string, Protein> tmp;
   Protein prt;
   string header;
   //int limit = 1000;
   int cnt = 0;
   int every=10;
   RefseqNameCleanup cleaner;
   while (prt.read(ifs, header)) {
      if ((cnt+1) % every == 0) {
         cout << "working on " << cnt+1 << endl;
         every *= 2;
      }
      //if (cnt > limit) break;
      string newname = cleaner(prt.getName());
      //cout << prt.getName() << " " << newname << endl;
      if (list.find(newname) != list.end()) {
         ++list[newname];
         tmp.insert(make_pair(newname, prt));
      }
      //else {
      //   cout << newname << " not in our wanted list\n";
      //}
      ++cnt;
   }
   cout << cnt << " sequences in " << fastaFile << endl;
   cout << tmp.size() << " found out of " << list.size() << endl;
   if (tmp.size() < cnt) {
      string missedFile;
      size_t x = fastaFile.rfind('/');
      if (x != string::npos) {
         missedFile = fastaFile.substr(x+1);
      }
      else missedFile = fastaFile;

      x=missedFile.rfind('.');
      if (x == string::npos) { 
         missedFile = missedFile + ".missed";
      }
      else {
         missedFile = missedFile.substr(0, missedFile.rfind('.')) + ".missed";
      }
      ofstream ofs(missedFile.c_str());
      cout << "the following ids are not found in " << fastaFile << endl;
      for (auto itr=list.begin(); itr != list.end(); ++itr) {
         if (itr->second == 0)
            ofs << itr->first << endl;
      }
      cout << "missed entries written to " << missedFile << endl;
   }
   return tmp;
}
template<class trimmer>
map<string, Protein> collectFasta(const string &fastaFile, map<string,int> &list, 
      trimmer stripper) {

   ifstream ifs(fastaFile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open " << fastaFile << endl;
      exit(1);
   }
   map<string, Protein> tmp;
   Protein prt;
   string header;
   //int limit = 1000;
   int cnt = 0;
   int every=10;
   //RefseqNameCleanup cleaner;
   while (prt.read(ifs, header)) {
      if ((cnt+1) % every == 0) {
         cout << "working on " << cnt+1 << endl;
         every *= 2;
      }
      //if (cnt > limit) break;
      string newname = stripper(prt.getName());
      //RefseqNameCleanup reftrimmer;
      //string newname = reftrimmer(prt.getName());
      //cout << prt.getName() << " " << newname << endl;
      if (list.find(newname) != list.end()) {
         ++list[newname];
         tmp.insert(make_pair(newname, prt));
      }
      //else {
      //   cout << newname << " not in our wanted list\n";
      //}
      ++cnt;
   }
   cout << cnt << " sequences in " << fastaFile << endl;
   cout << tmp.size() << " found out of " << list.size() << endl;
   if (tmp.size() < cnt) {
      string missedFile;
      size_t x = fastaFile.rfind('/');
      if (x != string::npos) {
         missedFile = fastaFile.substr(x+1);
      }
      else missedFile = fastaFile;

      if (missedFile.rfind('.') == string::npos) { 
         missedFile = missedFile + ".missed";
      }
      else {
         missedFile = missedFile.substr(0, missedFile.rfind('.')) + ".missed";
      }
      ofstream ofs(missedFile.c_str());
      cout << "the following ids are not found in " << fastaFile << endl;
      for (auto itr=list.begin(); itr != list.end(); ++itr) {
         if (itr->second == 0)
            ofs << itr->first << endl;
      }
      cout << "missed entries written to " << missedFile << endl;
   }
   return tmp;
}


//set<string> readList(const string &listFile) {
map<string,int> readList(const string &listFile) {
   ifstream inf(listFile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open file: " << listFile << endl;
      exit(1);
   }
   map<string, int> tmp;
   string id;
   int cnt=0;
   inf >> id;
   while (!inf.eof()) {
      ++cnt;
      tmp.insert(make_pair(id, 0));
      inf >> id;
   }
   cout << cnt << " id " << tmp.size() << " unique id\n";
   return tmp;
}

void alignMapping(const string &fname, const map<string, Protein> &pidbuff,
      const map<string, Protein> &refbuff) {
   ifstream ifs(fname.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open file: " << fname << endl;
      exit(1);
   }
   string resultFile("validMapping.tab");
   string logfile("invalidMapping.log");
   ofstream ofres(resultFile.c_str());
   ofstream oflog(logfile.c_str());

   Dynaln<ProteinScoreMethod> aligner;
   string line;
   getline(ifs, line);
   map<string, Protein>::const_iterator p1, p2;
   int cnt = 0;
   int bottomPos;
   int every=10;
   while (!ifs.eof()) {
      //if (cnt > 50) break;
      ++cnt;
      if (cnt % every == 0) {
         cout << "working on " << cnt << endl;
         every *= 2;
      }
      //cout << line << endl;
      getline(ifs, line);
      vector<string> row = split(line, '\t');
      p1 = pidbuff.find(row[0]);
      if (p1 == pidbuff.end()) {
         oflog << row[0] << " protein id not found\n";
         continue;
      }
      vector<string> refids = split(row[1], '|');
      for (int i=0; i<refids.size(); ++i) {
         p2 = refbuff.find(refids[i]);
         if (p2 != refbuff.end()) {
            Protein prt1(p1->second);
            Protein prt2(p2->second);
            //cout << prt1 << endl;
            //cout << prt2 << endl;
            aligner.setSeq(prt1, prt2);
            //aligner.runlocal();
            aligner.runglobal();
            pair<char, char> matched = aligner.getResiduesByTopPosition(stoi(row[3]), bottomPos);
            char bottomchar = (char)toupper(matched.second);
            if (matched.first != row[2][0]) {
               oflog << row[0] << " " << refids[i] << " " << row[2] << " " << row[3] 
                  << " " << bottomchar << " " << bottomPos << endl;
               oflog << "Top residue " << (char)matched.first << " different from input!\n";
               aligner.printAlign(oflog, 80);
            }
            else if (bottomchar == '-') {
               oflog << row[0] << " " << refids[i] << " " << row[2] << " " << row[3] 
                  << " " << bottomchar << " " << bottomPos << endl;
               oflog << "There is no bottom residue corresponding to top phosphorylation site\n";
               aligner.printAlign(oflog, 80);
            }
            else if (bottomchar != 'S' && bottomchar != 'T' && bottomchar != 'Y'
                  && bottomchar != 'D') {
               oflog << row[0] << " " << refids[i] << " " << row[2] << " " << row[3] 
                  << " " << bottomchar << " " << bottomPos << endl;
               oflog << "bottom residue not animal phosphorylation site\n";
               aligner.printAlign(oflog, 80);
            }
            else {
               ofres << row[0] << '\t' << refids[i] << '\t' 
                  << bottomchar << bottomPos << endl;
               //cout << "bottom residue: " << matched.second << " position: " << bottomPos << endl;
            }
         }
         else {
            oflog << row[1] << " refid not found\n";
         }
      }
   }
   cout << cnt << " rows processed. result written to " << resultFile << endl;
}

void writeBuffer(map<string, Protein> &buff, const string &file) {
   ofstream ofs(file.c_str());
   if (ofs.fail()) {
      cerr << "Failed to open " << file << " for writing\n";
      exit(1);
   }
   for (auto itr=buff.begin(); itr != buff.end(); ++itr) {
      (itr->second).setName(itr->first);
      ofs << itr->second;
   }
   cerr << "fasta in memeory written to file: " << file << endl;
}

void collectFasta(map<string, Protein> &pidbuffer, map<string, Protein> &refseqbuffer,
      const string &collectedpidFile, const string &collectedrefFile) {
   string rootdir="/ngs/ngs12/stefanRequest/PSP";
   string pidfile="substrate_pid.list";
   string reffile="substrate_ref.list";

   // protein sequence store
   string seqstore="/net/minerva/gcgblast";
   string uniref50File="uniref50.fasta";
   string refseqpFile="refseqp"; // reference sequence protein

   cout << "reading pid ...\n";
   map<string,int> pids = readList(rootdir + "/" + pidfile);
   cout << "reading ref ...\n";
   map<string,int> refids = readList(rootdir + "/" + reffile);

   RefseqNameCleanup ref_trimmer;
   PidNameCleanup pid_trimmer;
   // read refseq protein
   //map<string, Protein> refseqs = collectRefseqpFasta(seqstore + "/" + refseqpFile, refids);
   cout << "buffering ref sequuences ...\n";
   //map<string, Protein> refseqbuffer = collectFasta(seqstore + "/" + refseqpFile, 
   refseqbuffer = collectFasta(seqstore + "/" + refseqpFile, refids, ref_trimmer);
   cout << "buffering from hand downloaded refseq, missed from seqstoore ..\n";
   map<string, Protein> tmpbuffer = collectFasta(rootdir + "/refmissed.fas",
         refids, ref_trimmer);
   refseqbuffer.insert(tmpbuffer.begin(), tmpbuffer.end());

   cout << "buffering uniprot sequuences ...\n";
   //map<string, Protein> pidbuffer = collectFasta(seqstore + "/" + uniref50File, 
   pidbuffer = collectFasta(seqstore + "/" + uniref50File, pids, pid_trimmer);
   cout << "buffering swiss_prot sequuences ...\n";
   tmpbuffer = collectFasta(seqstore + "/swiss_prot", pids, pid_trimmer);
   pidbuffer.insert(tmpbuffer.begin(), tmpbuffer.end());
   cout << "buffering hand downloaded, missed in seqstoere ...\n";
   tmpbuffer = collectFasta(rootdir + "/pidmissed.fas", pids, pid_trimmer);
   pidbuffer.insert(tmpbuffer.begin(), tmpbuffer.end());

   // output collected sequences to file
   //string collectedpidFile = "pid.collect.fas";
   //string collectedrefFile = "ref.collect.fas";
   writeBuffer(pidbuffer, collectedpidFile);
   writeBuffer(refseqbuffer, collectedrefFile);
}

map<string, Protein> loadFasta(const string &fasfile) {
   ifstream ifs(fasfile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open " << fasfile << endl;
      exit(1);
   }
   map<string, Protein> tmp;
   Protein prt;
   string header;
   int cnt = 0;
   int every=10;
   while (prt.read(ifs, header)) {
      if ((cnt+1) % every == 0) {
         cout << "working on " << cnt+1 << endl;
         every *= 2;
      }
      tmp.insert(make_pair(prt.getName(), prt));
      ++cnt;
   }
   cout << cnt << " sequences in " << fasfile << endl;
   return tmp;
}

int main(int argc, char* argv[]) {
   string collectedpidFile = "pid.collect.fas";
   string collectedrefFile = "ref.collect.fas";
   ifstream testinf(collectedpidFile.c_str());
   map<string, Protein> pidbuffer, refseqbuffer;
   if (testinf.good()) {
      cout << "Loading from collection ...\n";
      pidbuffer = loadFasta(collectedpidFile);
      refseqbuffer = loadFasta(collectedrefFile);
   }
   else {
      cout << "Collecting from seqstore ...\n";
      collectFasta(pidbuffer, refseqbuffer, collectedpidFile, collectedrefFile);
   }

   string pid2reffile = "kinase_substrate_pid2ref.tab";
   alignMapping(pid2reffile, pidbuffer, refseqbuffer);


/*
   string infasfile("amplicon_nrlib.fas");
   if (argc > 1) {
      infasfile = argv[1];
   }
   string outfile = infasfile.substr(0, infasfile.rfind('.')) 
      + ".uniq.fas";
   ofstream ofs(outfile.c_str());
   if (ofs.fail()) {
      cerr << "Failed to open outfile " << outfile << endl;
      return 1;
   }

   set<bioseq> uniqueseq;
   bioseq bs;
   string header;
   ifstream ifs(infasfile.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open infasfile: " << infasfile << endl;
      return 1;
   }
   int cnt=0;
   while (bs.read(ifs, header)) {
      uniqueseq.insert(bs);
      ++cnt;
   }
   for (auto itr=uniqueseq.begin(); itr != uniqueseq.end(); ++itr) {
      ofs << *itr;
   }
   cout << uniqueseq.size() << " unique sequences out of " << cnt << "\n";
   */

   return 0;
}
