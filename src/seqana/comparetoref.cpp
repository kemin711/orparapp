#include <iostream>
#include <string>
#include <map>
#include <stack>
#include <iterator>

#include <alninfo.h>
#include <strformat.h>
#include <dynalnt.h>
#include <fastq.h>
#include "genoanalyzer.h"

#define DEBUG

using namespace std;
using namespace orpara;


/**
 * Simple method to compare reads to reference.
 * This algorithm tolerates large errors and gaps.
 * If the error rates are 2% or less, the BWA
 * algorithms should be used.
 */

//mutex mapupdate_mutex;
string nameAlninfoFile(const string &infile);
void printMap(const map<string, int> &info, ostream &ous);

/**
 * Global variable use by tabulateCodon for efficiency
 * Start 0-based index of the codon of interest in the reference seq.
 */
void displayAlninfo(const map<AlignInfo, int> &res, ostream &ous);
void aligntoref(const string &inputfile, const string &alninfofile, const string &refseqfile, const pair<float, float>& alncut, const vector<int> &codonIdx);
/**
 * read all fastq into a single container as DNA objects
 * for more complicated job, we need to use Fastq
 */
//vector<DNA> slurpFastq(const string &file, int lengthcut=150);
void aligntoref2(const string &inputfile, const string &alninfofile, 
      const string &refseqfile, const vector<int> &codonIdx, 
      int seqlencut, int readDirection, int refframe);
void aligntorefPair(const string &inputfile, const string &inputfile2,
      const string &alninfofile, const string &refseqfile, const vector<int> &codonIdx,
      int seqlencut);
void assignTargetPosition(vector<int> &pos, const string &posstr);
vector<int> readTargetIndex(const string &file);
map<string,string> readConfig(const string &file);
bool fileExists(const string &file);

void usage() {
   cerr << "comparetoref -r refseq.fas read.fastq result.info\n"
        << " You first need to convert bam files to fastq files.\n"
        << " Please use bam2cleanfastq program\n"
        << " This program will also generate length distribution\n"
        << "Options:\n"
        << "    -r <reference sequence fasta file>\n"
        << "    -i <input fastq file>\n"
        << "    -o <output information file> for statistics and analysis.\n"
        << "    -c identitycut coveragecut default 0.85 0.86\n"
        << "        nolonger using the refseq horizontal coverage\n"
        << "    -p codon position in the refseq nucleotide position\n"
        << "       in the format: 51,57,60,63 no space allowed\n"
        << "       The index is 0-based.\n"
        << "    -x <filename> target index file for the linked codon positions\n"
        << "       the index is the same 0-based as specified by -p\n"
        << "    --no-link flag to not count linked codons\n"
        << "    --pair forwardFile backwardFile for paired end reads.\n"
        << "       this will also set the sequence length cutoff to 45 nt.\n"
        << "    --libfile reference library file\n"
        << "    --conf config_file configuration file name\n"
        << "    -f frame 0,1,2.  Reading frame in the refseq. default 0.\n"
        << "\n";
}

int main(int argc, char* argv[]) {
   bool threadrun = false;
   string infile, alninfoFile, refseqfile, infile2;
   // for debug run
   //infile="test.fastq";
   //refseqfile="/d1/work/virology/bar4consensus.fas";
   // for testing!
   //infile="/ngs/ngs16/zhouke/virology/testrun/bar1.fastq";
   //infile="/ngs/ngs16/zhouke/virology/test.fastq";
   //infile="/ng14/zhouke/virology/testrun/bar10.fastq";
   //string inputdir="/ng14/zhouke/virology/ddl/run1/test5B/";
   //infile=inputdir + "test10S1R1.fastq";
   //infile2=inputdir + "test10S1R2.fastq";

   //pair<float,float> alncutoff = make_pair(0.85, 0.86);
   int fastqLengthCut = 180; // this should be 1/2 of refseq
   // in NS5A ORF 60=>91(L31), 256=>277(Y93H)
   //vector<int> targetCodonIndex{51, 57, 60, 63, 129, 141, 153, 243, 246};
   //0-based index in the amplicon
   vector<int> targetCodonIndex{60, 246}; // will be loaded from an external file
   string refposStr, targetIndexFile, libfile;
   //libfile = "/ngs/ngs16/zhouke/virology/hcvseq/amplicon_ulib";
   //for testing DDL long amplicon and shotgun short reads
   //libfile= "/ngs/ngs16/zhouke/virology/hcvseq/fdaref/hcv3NS5Alib.fas";
   string configFile="comparetoref.conf";
   int readDirection=0;
   int frame=0;
   int i = 1;
   while (i < argc) {
      if (!strcmp(argv[i], "-i")) infile = string(argv[++i]);
      else if (!strcmp(argv[i], "-o")) alninfoFile = string(argv[++i]);
      else if (!strcmp(argv[i], "-r")) refseqfile = string(argv[++i]);
      else if (!strcmp(argv[i], "-t")) threadrun = true;
      else if (!strcmp(argv[i], "-l")) fastqLengthCut = atoi(argv[++i]);
      else if (!strcmp(argv[i], "-f")) frame = atoi(argv[++i]);
      //else if (!strcmp(argv[i], "-c")) {
      //   float idencut = atof(argv[++i]);
      //   float covcut = atof(argv[++i]);
      //   alncutoff = make_pair(idencut, covcut);
      //}
      else if (!strcmp(argv[i], "-p")) refposStr = string(argv[++i]);
      else if (!strcmp(argv[i], "--help")) {
         usage(); return 1;
      }
      else if (!strcmp(argv[i], "-x")) targetIndexFile = string(argv[++i]);
      else if (!strcmp(argv[i], "--no-link")) targetCodonIndex.clear();
      else if (!strcmp(argv[i], "--reads-direction")) 
         readDirection=atoi(argv[++i]);
      else if (!strcmp(argv[i], "--pair")) {
         infile=string(argv[++i]);
         infile2=string(argv[++i]);
         fastqLengthCut=45;
      }
      else if (!strcmp(argv[i], "--libfile")) {
         libfile=string(argv[++i]);
      }
      else if (!strcmp(argv[i], "--conf")) {
         configFile=string(argv[++i]);
      }
      else {
         infile = string(argv[i]);
         if (i+1 < argc && argv[i+1][0] != '-') 
            alninfoFile = string(argv[++i]);
      }
      ++i;
   }
   cerr << "using configure file: " << configFile << endl;

   // configuration file does overwrite command line
   if (fileExists(configFile)) {
      map<string,string> configProperty=readConfig(configFile);
      if (libfile.empty())
         libfile=configProperty["reflibdir"] + "/" + configProperty["ref"];
      targetCodonIndex.clear();
      if (configProperty["targetIndex"] != "none") {
         vector<string> tmp=split(configProperty["targetIndex"], ',');
         for (size_t x=0; x<tmp.size(); ++x) {
            if (tmp[x].empty()) continue;
            targetCodonIndex.push_back(stoi(tmp[x]));
         }
      }
      fastqLengthCut=stoi(configProperty["lengthcut"]);
      string workdir;
      if (configProperty.find("workdirectory") != configProperty.end()) {
         if (configProperty["workdirectory"] == ".") {
            workdir.clear();
         }
         else {
            workdir=configProperty["workdirectory"];
         }
      }
      if (configProperty.find("frame") != configProperty.end()) {
         frame=atoi(configProperty["frame"].c_str());
      }
      if (configProperty.find("seqfiles") != configProperty.end()) {
         size_t x=configProperty["seqfiles"].find(',');
         if ( x != string::npos) {
            infile=configProperty["seqfiles"].substr(0,x);
            infile2=configProperty["seqfiles"].substr(x+1);
            if (!workdir.empty()) {
               infile = workdir + "/" + infile;
               infile2 = workdir + "/" + infile2;
            }
         }
         else {
            infile=configProperty["seqfiles"];
            if (!workdir.empty()) {
               infile = workdir + "/" + infile;
            }
         }
      }
   }

   if (infile.empty()) {
      cerr << "you must provide an input fastq file!\n";
      usage(); return 1;
   }
   if (!targetIndexFile.empty()) {
      targetCodonIndex=readTargetIndex(targetIndexFile);
   }
   /* this can be empty now.
   if (refseqfile.empty()) {
      cerr << "You must provide a reference file!\n";
      usage();
      return 1;
   }
   */
   // if user override the default positions
   if (!refposStr.empty()) assignTargetPosition(targetCodonIndex, refposStr);
   //copy(targetCodonIndex.begin(), targetCodonIndex.end(), ostream_iterator<int>(cerr, ", "));
   //cerr << endl << "after\n";

   if (alninfoFile.empty()) {
      if (!infile.empty()) {
         alninfoFile = nameAlninfoFile(infile);
      }
      else {
         usage();
         return 1;
      }
   }
   // for running on europa
   //GenotypeAnalyzer::setReflibFile("/ngs/ngs16/zhouke/virology/hcvseq/amplicon_ulib");

   //if (threadrun) {
   //   threadedAlign(infile, alninfoFile, refseqfile);
   //}
   //else {
      //aligntoref(infile, alninfoFile, refseqfile, alncutoff, targetCodonIndex);
   //}
   if (infile2.empty()) { // unpaired reads
      if (libfile.empty()) 
         libfile = "/ngs/ngs16/zhouke/virology/hcvseq/amplicon_ulib";
      GenotypeAnalyzer::setReflibFile(libfile);
      aligntoref2(infile, alninfoFile, refseqfile, targetCodonIndex, fastqLengthCut, readDirection, frame);
   }
   else { // paired ends
      //if (libfile.empty())
         //libfile= "/ngs/ngs16/zhouke/virology/hcvseq/fdaref/hcv3NS5Alib.fas";
         //libfile= "/ngs/ngs16/zhouke/virology/hcvseq/fdaref/hcv3NS5Blib.fas";
      GenotypeAnalyzer::setReflibFile(libfile);
      //targetCodonIndex.clear();
      //fastqLengthCut = 45;
      aligntorefPair(infile, infile2, alninfoFile, refseqfile, targetCodonIndex, fastqLengthCut);
   }

   return 0;
}

string getPairedFastqFileStem(const string &fname) {
   string stem = fname.substr(0, fname.rfind('.'));
   unsigned int i;
   if ((i=stem.rfind("R1")) != string::npos) {
      stem.erase(i, 2);
   }
   else if ((i=stem.rfind("R2")) != string::npos) {
      stem.erase(i, 2);
   }
   return stem;
}

string getReflibFileStem(const string &fname) {
   string::size_type i = fname.rfind('/');
   string stem;
   if (i != string::npos) {
      stem=fname.substr(i+1);
   }
   i = stem.rfind('.');
   if (i != string::npos) {
      stem = stem.substr(0, i);
   }
   return stem;
}

string buildStem(const string &pairedq, const string &libf) {
   return getPairedFastqFileStem(pairedq) + "_" + getReflibFileStem(libf);
}

// version for paired ends
void aligntorefPair(const string &inputfile, const string &inputfile2,
      const string &alninfofile, const string &refseqfile, const vector<int> &codonIdx,
      int seqlencut) {
   GenotypeAnalyzerPair analyzer;
   if (!codonIdx.empty()) {
      analyzer.setGenotypeIndex(codonIdx);
   }
   // reads direction is both by default
   analyzer.setBothDirections();
   if (!refseqfile.empty()) {
      DNA ref;
      ref.read(refseqfile);
      analyzer.setRefseq(ref); // force a particular refsequence.
   }
   // use the same cutoff for both length and alngnment
   GenotypeAnalyzer::setAlnlengthCut(seqlencut);
   analyzer.suckupReadPair(inputfile, inputfile2);
   string stem=buildStem(inputfile, analyzer.getReflibFile());
   string logFile = stem + ".log";
   ofstream olog(logFile.c_str());
   if (olog.fail()) {
      cerr << "Failed to open " << logFile << endl;
      exit(1);
   }
#ifdef DEBUG
   string alnfile = stem + ".aln";
   analyzer.setAlignOutputFile(alnfile);
#endif
   
   //analyzer.consume(reads);
   map<string,int> matchinfo = analyzer.consumePair();
   //analyzer.reorderResult();
   olog << analyzer.computeBestHaplotype();
   // will reset the reference sequence
   //analyzer.computeBestHaplotype();
   printMap(matchinfo, olog);
   cout << "Doing the second round of analysis ...\n";
   analyzer.fillgapYes();
   matchinfo = analyzer.consumePair();
   olog << "number of sequences passed align cutoff: " 
      << matchinfo["goodcnt"] << endl
      << "number of " << analyzer.getExpandedAlnlengthCutoff()
      << " nt or longer failed cutoff: " 
      << matchinfo["failcnt"] << endl; 
   printMap(matchinfo, olog);
   analyzer.printAlninfo(alninfofile);

   // output files

   string codonfile = stem + "codon.txt";
   string genofile = stem + "_genotype.txt";
   string confile = stem + "consensus.txt";
   string genoaafile = stem + "_genotypeaa.txt";
   string basefile = stem + "_base.txt";
   string qualfile = stem + "_qual.tab";
   //if (!codonIdx.empty()) {
      analyzer.reorderResult();
      analyzer.printResult(codonfile, genofile, confile, genoaafile, basefile, qualfile);
   //}
   //analyzer.computeBestHaplotype();
   olog << "second round: " << analyzer.computeBestHaplotype();
   olog << "identity cutoff: " << GenotypeAnalyzer::getIdentityCutoff()
      << " alnlen cutoff: " << GenotypeAnalyzer::getAlnlenCutoff() << endl;
   string donefile = inputfile.substr(0, inputfile.rfind('.') + 1) + "done";
   ofstream odone(donefile.c_str());
   if (odone.fail()) {
      cerr << "Failed to open " << donefile << " for writing\n";
      exit(1);
   }
   cerr << "genotype analysis done for " << inputfile 
      << " and " << inputfile2 << "\n";
}


// version 2
// TODO: Should be able to control whether to do one or two rounds
// of mapping
void aligntoref2(const string &inputfile, const string &alninfofile, 
      const string &refseqfile, const vector<int> &codonIdx,
      int seqlencut, int readDirection, int refframe) {
   GenotypeAnalyzer analyzer;
   if (!codonIdx.empty()) {
      analyzer.setGenotypeIndex(codonIdx);
   }
   if (readDirection > 0) {
      analyzer.setReadDirection(readDirection);
   }
   analyzer.setFrame(refframe);
   if (!refseqfile.empty()) {
      DNA ref;
      ref.read(refseqfile);
      analyzer.setRefseq(ref); // force a particular refsequence.
   }
   // use the same cutoff for both length and alngnment
   GenotypeAnalyzer::setAlnlengthCut(seqlencut);
   analyzer.suckupReads(inputfile);
   string logFile = inputfile.substr(0, inputfile.rfind('.') + 1) + "log";
   ofstream olog(logFile.c_str());
   if (olog.fail()) {
      cerr << "Failed to open " << logFile << endl;
      exit(1);
   }
#ifdef DEBUG
   string alnfile = inputfile.substr(0, inputfile.rfind('.') + 1) + "aln";
   analyzer.setAlignOutputFile(alnfile);
#endif
   analyzer.consumeForward();
   olog << analyzer.computeBestHaplotype();
   // will reset the reference sequence
   //analyzer.computeBestHaplotype();
   cout << "Doing the second round of analysis ...\n";
   analyzer.fillgapYes();
   pair<int,int> goodfail = analyzer.consumeForward();
   olog << "number of sequences passed align cutoff: " << goodfail.first << endl
      << "number of " << analyzer.getExpandedAlnlengthCutoff()
         << " nt or longer failed cutoff: " << goodfail.second << endl; 
   analyzer.printAlninfo(alninfofile);

   // output files
   string stem=inputfile.substr(0, inputfile.rfind('.'));
   string codonfile = stem  + "codon.txt";
   string genofile = stem + "_genotype.txt";
   string confile = stem + "consensus.txt";
   string genoaafile = stem + "_genotypeaa.txt";
   string basefile = stem + "_base.txt";
   string qualityfile = stem + "_qual.tab";
   if (!codonIdx.empty()) {
      analyzer.reorderResult();
      analyzer.printResult(codonfile, genofile, confile, genoaafile, basefile, qualityfile);
   }
   //analyzer.computeBestHaplotype();
   olog << "second round: " << analyzer.computeBestHaplotype();
   olog << "identity cutoff: " << GenotypeAnalyzer::getIdentityCutoff()
      << " alnlen cutoff: " << GenotypeAnalyzer::getAlnlenCutoff() << endl;
   string donefile = inputfile.substr(0, inputfile.rfind('.') + 1) + "done";
   ofstream odone(donefile.c_str());
   if (odone.fail()) {
      cerr << "Failed to open " << donefile << " for writing\n";
      exit(1);
   }
   cerr << "genotype analysis done for " << inputfile << "\n";
}

/* have the same problem as when using template class as parameter direclty.
 * I cannot hide the syntax error for some reason.
 *
void threadedAlign(const string &infile, const string &alninfoFile, const string &refseqfile) {
   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << infile << endl;
      return;
   }
   thread t1, t2, t3;
   map<AlignInfo, int> result1, result2, result3, result4;
   t1 = thread(mapref, inf, refseqfile, 1, result1);
   t2 = thread(mapref, inf, refseqfile, 2, result2);
   t3 = thread(mapref, inf, refseqfile, 3, result3);
   mapref(inf, refseqfile, 4, result4);
   t1.join();
   t2.join();
   t3.join();
   cout << "all jobs done, combine results\n";
   result4.insert(result1.begin(), result1.end());
   result4.insert(result2.begin(), result2.end());
   result4.insert(result3.begin(), result3.end());

   ofstream ouf("threadResult.tab");
   displayAlninfo(result4, ouf);
   ouf.close();
   cerr << "result written to threadResult.tab\n";
}

void mapref(istream &ins, const string &refile, int ith, map<AlignInfo, int> &result) {
   DNA ref;
   ref.read(refile);
   SimpleScoreMethod scoreMethod(10, -9, -11, -3);
   Dynaln<SimpleScoreMethod> align(scoreMethod);
   align.setSeq1(ref);

   Fastq fasq;
   int cnt=0;
   int every=10;
   while (fasq.read(ins) && cnt % ith == 0) {
      ++cnt;
      DNA raw(fasq.getName(), fasq.getSequence());
      align.setSeq2(raw);
      align.runlocal();
      ++result[AlignInfo(align.getScore(), floor(align.getIdentity()*100)/100, 
            floor(align.getCov1()*100)/100)];
      if (cnt % every == 0) {
         cerr << "working on " << cnt << endl
            << "   identity: " << floor(align.getIdentity()*100)/100 << " coverage: " 
            << floor(100*align.getCov1())/100 << endl;
         every *= 2;
      }
   }
   cout << cnt << " sequences processed\n";
}
*/
/*
void doAlign(Dynaln<SimpleScoreMethod> &aligner, const DNA *seq, map<AlignInfo, int> &res) {
   aligner.setSeq2(*seq);
   aligner.runlocal();
   mapupdate_mutex.lock();
   ++res[AlignInfo(aligner.getScore(), roundPercent(aligner.getIdentity()), 
               roundPercent(aligner.getCov1()),
               seq->length())];
   mapupdate_mutex.unlock();
}
*/

void foo(Dynaln<SimpleScoreMethod> &a) {
   // for testing only
}

/**
 * The threaded version is not working because of bug in gnu 4.8.1
 * I have also seen other poepole posting the same bug in 4.8.2.
 * I may have to use 4.9 in the future.
 */
/*
void aligntorefThread(const string &inputfile, const string &alninfofile, const string &refseqfile) {
   DNA ref;
   ref.read(refseqfile);

   Fastq fasq;
   ifstream inf(inputfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << inputfile << endl;
      return;
   }

   SimpleScoreMethod scoreMethod(10, -9, -11, -3);
   //cerr << "score method made\n";
   Dynaln<SimpleScoreMethod> align(scoreMethod);
   Dynaln<SimpleScoreMethod> align1(scoreMethod);
   Dynaln<SimpleScoreMethod> align2(scoreMethod);
   Dynaln<SimpleScoreMethod> align3(scoreMethod);
   //cerr << "aligner intialized\n";
   align.setSeq1(ref);
   cerr << "refseq set on aligner\n";

   // the test code is fine
   thread newthread(foo, &align);
   newthread.join();

   map<AlignInfo, int> result;
   int cnt=0;
   int every=10;
   thread trd1, trd2, trd3;
   DNA *raw1, *raw2, *raw3, *raw;
   while (fasq.read(inf)) {
      ++cnt;
      if (cnt % every == 0) {
         cerr << "working on " << cnt << endl;
         every *= 2;
      }
      DNA *raw = new DNA(fasq.getName(), fasq.getSequence());
      if (cnt % 1 == 0) {
         raw1 = new DNA(fasq.getName(), fasq.getSequence());
         trd1 = thread(doAlign, align1, raw1, result);
      }
      else if (cnt % 2 == 0) {
         raw2 = new DNA(fasq.getName(), fasq.getSequence());
         trd2 = thread(doAlign, align1, raw2, result);
      }
      else if (cnt % 3 == 0) {
         raw3 = new DNA(fasq.getName(), fasq.getSequence());
         trd3 = thread(doAlign, align1, raw3, result);
      }
      else { // main thread do some work
         raw = new DNA(fasq.getName(), fasq.getSequence());
         align.setSeq2(*raw);
         align.runlocal();

         trd1.join();
         trd2.join();
         trd3.join();

         ++result[AlignInfo(align.getScore(), floor(align.getIdentity()*100)/100, 
               floor(align.getCov1()*100)/100)];
         delete raw1; delete raw2; delete raw3; delete raw;
         raw1=raw2=raw3=0;
      }
   }
   // deal with end condition
   if (raw1 == 0) { // done
   }
   else if (raw2 == 0) { // raw1 done
      trd1.join();
      delete raw1;
   }
   else if (raw3 == 0) {
      trd1.join(); trd2.join();
      delete raw1; delete raw2;
   }
   else if (raw == 0) {
      trd1.join(); trd2.join(); trd3.join();
      delete raw1; delete raw2; delete raw3;
   }

   cout << cnt << " sequences processed\n";

   ofstream ouf(alninfofile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open file: " << alninfofile << " for output!\n";
      exit(1);
   }
   displayAlninfo(result, ouf);
}
*/

// non threaded version. The compiler is not friendly to
// threading with template classes at this point
/*
void aligntoref(const string &inputfile, const string &alninfofile, 
      const string &refseqfile, const pair<float, float> &alncut,
      const vector<int> &codonIdx) {
   DNA ref;
   ref.read(refseqfile);

   Fastq fasq;
   ifstream inf(inputfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input file: " << inputfile << endl;
      return;
   }
#ifdef DEBUG
   string alnfile = inputfile.substr(0, inputfile.rfind('.') + 1) + "aln";
   ofstream oaln(alnfile.c_str());
   if (oaln.fail()) {
      cerr << "Failed to open " << alnfile << " for writting!\n";
   }
#endif

   // is producing double gaps, gap score not low enough
   //SimpleScoreMethod scoreMethod(10, -9, -11, -3);
   SimpleScoreMethod scoreMethod(10, -9, -19, -3);
   Dynaln<SimpleScoreMethod> align(scoreMethod);
   align.setSeq1(ref);

   map<AlignInfo, int> alnsummary;
   // codons with missing base or gaps will be considered 
   // uncertain.
   GenotypeAnalyzer analyzer(ref.getSequence(), codonIdx);

   int cnt=0;
   int goodcnt=0;
   int every=10;
   while (fasq.read(inf)) {
      ++cnt;
      DNA raw(fasq.getName(), fasq.getSequence());
      align.setSeq2(raw);
      align.runlocal();
      if (align.getIdentity() > alncut.first && align.getCov1() > alncut.second) {
         ++goodcnt;
      }
      //cerr << "working on read: " << fasq.getName() << endl;
      if (!codonIdx.empty() && align.getIdentity() > 0.85 && align.getSeq1AlignedLength() > 200) {
         analyzer.accumulate(align.getTopAln(), align.getBottomAln(), (unsigned int)(align.topBeginIndex()));
      }
      ++alnsummary[AlignInfo(align.getScore(), roundPercent(align.getIdentity()), 
            roundPercent(align.getCov1()), raw.length())];
#ifdef DEBUG
      align.printAlign(oaln, 80);
#endif
      if (cnt % every == 0) {
         cerr << "working on " << cnt << " " << fasq.getName()
            << ":  identity=" << roundPercent(align.getIdentity()) 
            << " coverage=" << roundPercent(align.getCov1()) << endl;
         every *= 2;
      }
   }
   cout << cnt << " sequences processed\n"
      << goodcnt << " alighments with identity > " << alncut.first
      << " ref_coverage > " << alncut.second << endl;

   ofstream ouf(alninfofile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open file: " << alninfofile << " for output!\n";
      exit(1);
   }
   displayAlninfo(alnsummary, ouf);
   inf.close();
   ouf.close();
#ifdef DEBUG
   oaln.close();
   cerr << "alignment file written to " << alnfile << endl;
#endif
   // output files
   string codonfile = inputfile.substr(0, inputfile.rfind('.')) + "codon.txt";
   string genofile = inputfile.substr(0, inputfile.rfind('.')) + "_genotype.txt";
   string confile = inputfile.substr(0, inputfile.rfind('.')) + "consensus.txt";
   string genoaafile = inputfile.substr(0, inputfile.rfind('.')) + "_genotypeaa.txt";
   string basefile = inputfile.substr(0, inputfile.rfind('.')) + "_base.txt";
   if (!codonIdx.empty()) {
      analyzer.reorderResult();
      analyzer.printResult(codonfile, genofile, confile, genoaafile, basefile);
   }
   cerr << "analysis done\n";
}
*/

void displayAlninfo(const map<AlignInfo, int> &res, ostream &ous) {
   ous << "score\tidentity\tcoverage\treadlen\tcount\n";
   /* C98 version
   map<AlignInfo, int>::const_iterator it = res.begin();
   while (it != res.end()) {
      ous << it->first << " => " << it->second << endl;
      ++it;
   }
   C2011 version
   */
   for (auto itr = res.cbegin(); itr != res.cend(); ++itr) {
      ous << itr->first << "\t" << itr->second << endl;
   }
}

/**
 * Info for statistics.
 */
string nameAlninfoFile(const string &infile) {
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      return infile.substr(0,i+1) + "info";
   }
   else {
      return infile + ".info";
   }
}

void printMap(const map<string, int> &info, ostream &ous) {
   for (auto ptr=info.begin(); ptr != info.end(); ++ptr) {
      ous << ptr->first << "\t" << ptr->second << endl;
   }
}

void assignTargetPosition(vector<int> &pos, const string &posstr) {
   //cerr << posstr << endl;
   if (!pos.empty()) pos.clear();
   vector<string> posv = split(posstr, ',');
   //copy(posv.begin(), posv.end(), ostream_iterator<string>(cerr, ", "));
   //cerr << endl;
   for (size_t i=0; i<posv.size(); ++i) {
      pos.push_back(atoi(posv[i].c_str()));
   }
}

/*
void doAlign(istream &ins, const string &refile, int ith, 
      map<AlignInfo, int> &result);
*/

/**
   in NS5A ORF 60=>91(L31), 256=>277(Y93H)
   vector<int> targetCodonIndex{51, 57, 60, 63, 129, 141, 153, 243, 246};
   The target index is the index in the amplicon which may be different
   from the index when counting from the start of the ORF.
   It is the DNA index.
   The line starting with # will be considered comment and will be discarded.
   The index is one or more numbers on the same or different lines.
for example:
-----
# this is comment
60 256
-----
*/

vector<int> readTargetIndex(const string &file) {
   ifstream inf(file.c_str());
   if (inf.fail()) {
      cerr << "Failed to open target index position file: " << file << endl;
      exit(1);
   }
   cout << "reading index from " << file << endl;
   string ln;
   getline(inf, ln);
   vector<int> result;
   while (!inf.eof()) {
      cout << "processing line: " << ln << endl;
      if (ln[0] != '#') {
         vector<string> tmp = dissect(ln);
         for (unsigned int i=0; i<tmp.size(); ++i) {
            if (tmp[i].empty()) continue;
            cout << stoi(tmp[i]) << " ";
            result.push_back(stoi(tmp[i]));
         }
      }
      getline(inf, ln);
   }
   cout << endl;
   inf.close();
   //cout << "numbers stored in file:\n";
   //copy(result.begin(), result.end(), ostream_iterator<int>(cout, " "));
   return result;
}


/*
vector<DNA> slurpFastq(const string &file, int lengthcut) {
   cerr << "Reading fastq into memory ...\n";
   ifstream inf(file.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << file << endl;
      exit(1);
   }
   vector<DNA> tmp;
   if (inf.fail()) {
      cerr << "Failed to open input file: " << file << endl;
      return tmp;
   }
   Fastq fasq;
   int cnt=0;
   int total=0;
   while (fasq.read(inf)) {
      ++total;
      if (fasq.length() > lengthcut) {
         ++cnt;
         // remove the @ from fastq names
         DNA raw(fasq.getName().substr(1), fasq.getSequence());
         tmp.push_back(raw);
      }
   }
   cerr << tmp.size() << " Fastq sequences longer than " << lengthcut
      << " pulled into memory out of " << total << " from file: " << file << "\n";
   return tmp;
}
*/

map<string,string> readConfig(const string &file) {
   //string conf="comparetoref.conf";
   cerr << "reading config file: " << file << endl;
   ifstream ifs(file.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open comparetoref configuration file: " << file << endl;
      exit(1);
   }
   map<string, string> confmap;
   string ln;
   getline(ifs, ln);
   size_t s;
   while (!ifs.eof()) {
      if (ln[0] != '#' && !ln.empty()) {
         //cerr << ln << endl;
         while ((s=ln.rfind(' ')) != string::npos) {
            ln.erase(s);
         }
         s=ln.find('=');
         confmap[ln.substr(0,s)]=ln.substr(s+1);
      }
      getline(ifs, ln);
   }
   // print out, for debug
   cerr << endl;
   for (auto ptr=confmap.begin(); ptr != confmap.end(); ++ptr) {
      cerr << ptr->first << " => " << ptr->second << endl;
   }
   return confmap;
}

bool fileExists(const string &file) {
   ifstream ifs(file.c_str());
   if (ifs.fail()) {
      return false;
   }
   return true;
}
