#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <map>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>

#include <bioseq.h>
#include <strformat.h>
#include <dynalnt.h>

#include "nucaln.h"
#include "MismatchMapping.h"

using namespace std;
using namespace orpara;

class QCInfo {
   public:
      QCInfo() : leftIdentity(0), rightIdentity(0), stemIdentity(0), stemLen(0) { }
      QCInfo(float left, float right, float center, int len) : leftIdentity(left), rightIdentity(right), stemIdentity(center), stemLen(len) { }
      QCInfo(const QCInfo &info) : leftIdentity(info.leftIdentity), rightIdentity(info.rightIdentity),
         stemIdentity(info.stemIdentity), stemLen(info.stemLen) { }
      void setLeftIdentity(float left) { leftIdentity=left; }
      void setRightIdentity(float right) { rightIdentity=right; }
      void setStemIdentity(float center) { stemIdentity=center; }
      void setStemLength(int len) { stemLen = len; }
      int getStemLength() const { return stemLen; }
      float getLeftIdentity() const { return leftIdentity; }
      float getRightIdentity() const { return rightIdentity; }
      float getStemIdentity() const { return stemIdentity; }
      QCInfo& operator=(const QCInfo& info);
      bool operator<(const QCInfo &other) const;
      friend ostream& operator<<(ostream &ous, const QCInfo &info);
      bool sameIdentity(const QCInfo &info) const { return leftIdentity == info.leftIdentity 
         && rightIdentity == info.rightIdentity && stemIdentity == info.stemIdentity; }
      ostream& printIdentity(ostream &ous) const {
         ous << leftIdentity << '\t' << rightIdentity << '\t' << stemIdentity;
         return ous;
      }

   private:
      float leftIdentity;
      float rightIdentity;
      float stemIdentity;
      int stemLen;
};

/**
 * @return the length of the stem loop alignment.
 * 0 if no stem loop found
 * Only update first,second if stem loop quality good enough
 */
template<class T>
int getStemSequence(const Dynaln<T> &aln, string &first, string &second);

QCInfo& QCInfo::operator=(const QCInfo &info) {
   if (this != &info) {
      leftIdentity = info.leftIdentity;
      rightIdentity = info.rightIdentity;
      stemIdentity = info.stemIdentity;
      stemLen = info.stemLen;
   }
   return *this;
}

bool QCInfo::operator<(const QCInfo &other) const { 
   if (leftIdentity < other.leftIdentity) return true;
   else if (leftIdentity > other.leftIdentity) return false;
   if (rightIdentity < other.rightIdentity) return true;
   else if (rightIdentity > other.rightIdentity) return false;
   if (stemIdentity < other.stemIdentity) return true;
   else if (stemIdentity > other.stemIdentity) return false;
   if (stemLen < other.stemLen) return true;
   else return false;
}

ostream& operator<<(ostream &ous, const QCInfo &info) {
   ous << info.leftIdentity << '\t' << info.rightIdentity << '\t' << info.stemIdentity
      << '\t' << info.stemLen;
   return ous;
}

// helper function put here for compiler to find
void showWarning(const DNA &seq, const vector<int> &pos, 
   const vector<char> &mk, const string &message, ostream &ouf) {
   string filler(seq.length(), '-');
   for (size_t i = 0; i < pos.size(); ++i) {
      filler.replace(pos[i], 1, string(1, mk[i]));
   }
   ouf << message << endl
      << seq.toString() << endl << filler << endl;
}
int updateRefcount(map<string, int> &rfc, const string &first, const string &second,
      ostream &ouf);

template<class T>
bool alignWith(Dynaln<T> &aln, const bioseq &seq, ostream &ous, float identityCut) {
   //cout << "aligning with " << name << "\n" << seq << endl;
   //DNA rawseq(name, seq); // when it is out of scope it is gone!
   aln.setSeq2(seq);
   aln.runlocal();

   //cout << "identity: " << aln.getIdentity() << endl;
   if (aln.getIdentity() > identityCut && aln.getCov1() > 0.65) {
      aln.printAlign(ous);
      return true;
   }
   return false;
}


template<class T>
bool alignWithNoOutput(Dynaln<T> &aln, const bioseq &seq, float identityCut) {
   aln.setSeq2(seq);
   aln.runlocal();
   if (aln.getIdentity() > identityCut && aln.getCov1() > 0.65) {
      return true;
   }
   return false;
}

/**
 * Fast version, hard code cutoff values.
 */
template<class T>
bool alignWithFast(Dynaln<T> &aln, const bioseq &seq) {
   aln.setSeq2(seq);
   aln.runlocal();

   if (aln.getIdentity() > 0.72 && aln.getCov1() > 0.4 && aln.getAlnlen() > 9) 
   {
      return true;
   }
   return false;
}

float round2(double fn) {
   return (float)(floor(fn * 100 + 0.5)/100);
}

template<class T1, class T2> void printMap(const map<T1, T2> &m, ostream &ous);
/**
 * read reference into map
 * File content like the following
 * AAAGAGATGCAAAGGATCCTC,TRCN0000010402
 * AAAGAGGCTCCTACCAACGAC,TRCN0000010300
 * AAATATGAATCTGCTCACTTC,TRCN0000040153
 * AAATGATGTGTGGGAGGTTAC,TRCN0000010181
 * ..... many lines
 */
void readReference(const string &file, map<string, int> &refcount);
int numdiff(const string &longstr, const string &shortstr);
int numdiffTail(const string &longstr, const string &shortstr);
void generateTrimmedFastaFile(const DNA &left, const DNA &right, const string &bamfile,
      const int leftFlankLength, const int rightFlankLength,
      const float identityCut);
int tallyStemLoop(const DNA &left, const DNA &right, const string &refile, 
      const string &inputfile, const string &outfile, float identityCut);

template<class T>
bool findStemLoopLeft(const DNA &seq, const int e, ostream &ouf, Dynaln<T> &aligner, 
      map<string, int> &rfc, QCInfo &qci, map<pair<string,string>,int> &doubleHit) {
   ouf << "Align to right only, Cannot find left margin!\n";
   if (e < 25) {
      ouf << "5' flanking too short!\n"
         << seq << endl;
      return false;
   }
   int bi = e - 48;
   if (bi < 0) {
      bi = 0;
   }
   int citmp = e-29;
   if (citmp < 0) citmp = 0;

   int ci = seq.locateSubsequence("CTCGAG", citmp, e);
   if (ci == -1) {
      ci = (e - 24); // set at imaginary location
      vector<int> marks;
      vector<char> markchar;
      marks.push_back(bi);
      markchar.push_back('?');
      marks.push_back(ci);
      markchar.push_back('?');
      marks.push_back(e);
      markchar.push_back('|');
      showWarning(seq, marks, markchar, "No Xho I Site", ouf);
      if (qci.getRightIdentity() < 0.81) {
         ouf << "Right identity too low for further consideration\n";
         return false;
      }
   }
   else { // with Xho I site
      ci += 3;
   }

   try {
      DNA top = seq.subsequence(bi, ci-bi);
      DNA bottom = seq.subsequence(ci, 24);
      bottom.revcomp();
      aligner.setSeq1(top);
      aligner.setSeq2(bottom);
      //aligner.runlocal();
      aligner.runglobal();
      ouf << "top seq:\n" << top.toString() << endl 
         << "bottom seq:\n" << bottom.toString() << endl;
      ouf << "alignment:\n";
      aligner.printAlign(ouf);
      if (aligner.getIdentity() < 0.78) {
         ouf << "poor stem align, ignored\n";
         return false;
      }

      string first, second;
      qci.setStemLength(getStemSequence(aligner, first, second));
      qci.setStemIdentity(round2(aligner.getIdentity()));
      if (qci.getStemLength() < 24) {
         ouf << "stem: " << qci.getStemLength() << " not long enough\n";
         return false;
      }
      else if (qci.getStemIdentity() > 0.79) {
         if (updateRefcount(rfc, first, second, ouf) == 2) 
            ++doubleHit[make_pair(first,second)];
      }
   }
   catch (AlnInputException &e) {
      cerr << e.what() << endl << "Failed findStemloopLeft\n";
      exit(1);
   }

   return true;
}


/**
 * @param aligner is the Dynaln used to align the two stems.
 * @return the length of the tail for informational purpose.
 * From the right side of 5'-primer match.
 * return the length of the tail.
 */
template<class T>
int findStemLoopRight(const DNA &seq, const int b, ostream &ouf, Dynaln<T> &aligner,
      map<string, int> &rfc, QCInfo &qci, map<pair<string, string>, int> &doubleHit) {
   int tailLen = seq.length() - b;

   ouf << "Align to left only, Stemloop no right margin!\n"
      << "seq.length: " << seq.length() << " insert b: " << b << endl;
   if (tailLen < 24 || qci.getLeftIdentity() < 0.78) {
      ouf << "insert length: " << tailLen << " too short or identity too low, not considered.\n";
      qci.setStemIdentity(0);
      qci.setStemLength(0);
      return tailLen;
   }

   int ci = seq.locateSubsequence("CTCGAG", b + 18);
   if (ci == -1) { // no Xho I site
      ci = (b + 24); // set at imaginary location
      vector<int> marks;
      vector<char> markchar;
      marks.push_back(b);
      markchar.push_back('|');
      marks.push_back(ci);
      markchar.push_back('?');
      showWarning(seq, marks, markchar, "No Xho I Site", ouf);
      if (qci.getLeftIdentity() < 0.85) {
         ouf << "Left Identity too low, no considering further\n";
         qci.setStemIdentity(0);
         qci.setStemLength(0);
         return tailLen;
      }
   }
   else { // has Xho I site
      ci += 3;
   }

   string first, second;
   DNA top = seq.subsequence(b, ci-b);
   if (ci >= seq.length()) {
      ouf << "Xho I site " << ci << " outside sequence. Only got half of stem loop!\n";
      ouf << top.toString() << endl;
      if (qci.getLeftIdentity() < 0.87) {
         ouf << "Left identity too low combined with Xho I side outside sequence\n";
         qci.setStemIdentity(0);
         qci.setStemLength(0);
         return tailLen;
      }
      first = top.substring(0, 24);
      updateRefcount(rfc, first, second, ouf);
      qci.setStemIdentity(0);
      qci.setStemLength(ci-b);
   }
   else { // get bottom seq and do the alignment
      DNA bottom = seq.subsequence(ci, 24);
      bottom.revcomp();
      if (bottom.length() < 4) {
         ouf << "Bottom sequence too short, not worth analysis\n";
         if (qci.getLeftIdentity() < 0.86) {
            qci.setStemIdentity(0);
            qci.setStemLength(0);
            return tailLen;
         }
         
         first=top.substring(0,24);
         updateRefcount(rfc, first, second, ouf); // these is no bottom stem
         qci.setStemIdentity(0);
         qci.setStemLength(ci-b);
      }
      else {
         try {
            aligner.setSeq1(top);
            aligner.setSeq2(bottom);
            //aligner.runlocal();
            aligner.runglobal();
            ouf << "top seq:\n" << top.toString() << endl 
               << "bottom seq:\n" << bottom.toString() << endl;
            ouf << "alignment:\n";
            aligner.printAlign(ouf);
            if (aligner.getIdentical() == 0) {
               // very poor sequence such as 
               // ACCAAAAAACAACCCCCCCCCCCC------------------------
               // ------------------------GGGTTTTTTTGGGGGGGGGGGGTT
               // The reads has only A and C
               ouf << "garabage alignment ignored\n";
               return tailLen;
            }
            else if (aligner.getNoTerminalGapIdentity() < 0.81) {
               ouf << "stem align too poor to be considered\n";
               return tailLen;
            }
            qci.setStemLength(getStemSequence(aligner, first, second));
            qci.setStemIdentity(round2(aligner.getIdentity()));
            if (qci.getStemLength() < 24) {
               ouf << "Stem: " << qci.getStemLength() << " not long enough\n";
            }
            else if (updateRefcount(rfc, first, second, ouf) == 2) 
               ++doubleHit[make_pair(first,second)];
         }
         catch (AlnInputException &e) {
            cerr << e.what() << endl << " failed findStemLoopRight\n";
            exit(1);
         }
      }
   }
   return tailLen;
}

/**
 * @return the length of the insert for informational purpose.
 */
template<class T>
int findStemLoop(const DNA &seq, const int b, const int e, ostream &ouf, Dynaln<T> &aligner,
      map<string, int> &rfc, QCInfo &qci, map<pair<string, string>, int> &doubleHit) {
   int insertLen = e - b + 1;
   if (insertLen < 1) {
      ouf << "primer dimer, ignored\n";
      qci.setStemIdentity(0);
      qci.setStemLength(0);
      return insertLen;
   }
   else if (insertLen < 24) {
      //return false;
      qci.setStemIdentity(0);
      qci.setStemLength(0);
      ouf << "insert too short ignored\n";
      return insertLen;
   }

   ouf << "finding stem loop structure ...\n"
      << "b e: " << b << " " << e << endl;
   int ci = seq.locateSubsequence("CTCGAG", b+18, e);
   // simple approach
   if (ci == -1) { // no Xho I site
      vector<int> marks;
      vector<char> markchar;
      ci = (b + e + 1)/2; // set at imaginary location
      marks.push_back(b);
      markchar.push_back('|');
      marks.push_back(e);
      markchar.push_back('|');
      marks.push_back(ci);
      markchar.push_back('?');
      showWarning(seq, marks, markchar, "No Xho I Site", ouf);
   }
   else ci += 3;

   string first, second;
   if (e-ci+1 == 0) {
      //cerr << "Second stem empty: " << b << " " << ci << " " << e << endl;
      vector<int> marks;
      vector<char> markchar;
      marks.push_back(b);
      markchar.push_back('|');
      marks.push_back(e);
      markchar.push_back('|');
      marks.push_back(ci);
      markchar.push_back('|');
      showWarning(seq, marks, markchar, "Second stem empty", ouf);
      ouf << seq.substring(b, ci-b) << endl;
      first = seq.substring(b, 24);
      updateRefcount(rfc, first, second, ouf); // no bottom stem seq
      qci.setStemIdentity(0);
      qci.setStemLength(e-b+1);
      return insertLen;
   }
   // Now ci is marking the start of the second half in ideal situation.
   //              ci
   //              |
   // b--21nt--CTC*GAG--21nt--e
   try {
      DNA top = seq.subsequence(b, ci-b);
      DNA bottom = seq.subsequence(ci, e-ci+1);
      //cerr << "top and bottom DNA:\n" << top << bottom << endl;
      bottom.revcomp();
      aligner.setSeq1(top);
      aligner.setSeq2(bottom);
      ouf << "top seq:\n" << top.toString() << endl 
         << "bottom seq:\n" << bottom.toString() << endl;
      int extra = bottom.length() - 24;
      if (extra > 0 && boost::starts_with(bottom.toString(), string(extra, 'A'))) {
         bottom = bottom.subsequence(extra);
         ouf << "bottom trimmed to " << bottom.toString() << endl;
      }
      ouf << "alignment:\n";
      aligner.runglobal();
      aligner.printAlign(ouf);

      qci.setStemIdentity(round2(aligner.getIdentity()));
      qci.setStemLength(getStemSequence(aligner, first, second));
      if (qci.getStemLength() < 24) {
         ouf << "stem: " << qci.getStemLength() << " not long enough\n";
      }
      else if (qci.getStemIdentity() > 0.78) {
         if (updateRefcount(rfc, first, second, ouf) == 2) 
            ++doubleHit[make_pair(first, second)];
      }
   }
   catch (AlnInputException &e) {
      cerr << "Bad imput sequence in findStemLoop() " << e.what() << endl;
      return insertLen;
   }
   return insertLen;
}

/**
 * @return the length of the stem loop alignment.
 * 0 if no stem loop found
 * Only update first,second if stem loop quality good enough
 */
template<class T>
int getStemSequence(const Dynaln<T> &aln, string &first, string &second) {
   Nucaln nucaln;

   //if (aln.getIdentity() > 0.88 && (!boost::ends_with(aln.getTopAln(), "CTC")
   if (!boost::ends_with(aln.getTopAln(), "CTC") || !boost::ends_with(aln.getBottomAln(), "CTC")) {
      if (aln.getNoTerminalGapIdentity() < 0.82) {
         first.clear();
         second.clear();
         //cout << "Junk aln discard\n";
         return 0;
      }
      else {
         //cout << "Not ending with CTC with potential fix for loop\n";
         //   aln.printAlign(cout);

         if (boost::ends_with(aln.getTopAln(), "CTC") && boost::ends_with(aln.getBottomAln(), "CCC")) {
            string tmp = aln.getBottomAln();
            tmp.replace(tmp.length()-2, 1, "T");
            nucaln = Nucaln(aln.getTopAln(), tmp);
         }
         else if (boost::ends_with(aln.getBottomAln(), "CTC") && boost::ends_with(aln.getTopAln(), "CCC")) {
            string tmp = aln.getTopAln();
            tmp.replace(tmp.length()-2, 1, "T");
            nucaln = Nucaln(tmp, aln.getBottomAln());
         }
      }
   }

   if (aln.getAlnlen() == 24) {
      if (aln.getIdentity() > 0.999) {
         first=aln.getTopAln();
         second.clear();
         return 24;
      }
      else {
         nucaln = Nucaln(aln.getTopAln(), aln.getBottomAln());
      }
   }
   else if (aln.getAlnlen() == 23) {
      if (boost::ends_with(aln.getTopAln(), "CTC") && 
            boost::ends_with(aln.getBottomAln(), "CTC")) {
         //cout << "Little bit shorter with good loop:\n";
         //aln.printAlign(cout);
         if (aln.topAlignBeginIndex() == 1 && aln.bottomAlignBeginIndex() == 0) {
            //cout << "missing top one base\n";
            nucaln = Nucaln(aln.getTopSequence().substring(0,1) + aln.getTopAln(), "-" + aln.getBottomAln());
         }
         else if (aln.bottomAlignBeginIndex() == 1 && aln.topAlignBeginIndex() == 0) {
            //cout << "missing bottom one base\n";
            nucaln = Nucaln("-" + aln.getTopAln(), aln.getBottomSequence().substring(0,1) + aln.getBottomAln());
         }
         else if (aln.topAlignBeginIndex() == 0 && aln.bottomAlignBeginIndex() == 0) {
            // no where to find missing base
            return 23;
         }
         else {
            cout << aln.topAlignBeginIndex() << " " << aln.bottomAlignBeginIndex() << endl;
            cout << "Nothing is done!\n";
            exit(1);
         }
      }
   }
   else if (aln.getAlnlen() > 24) {
      nucaln = Nucaln(aln.getTopAln(), aln.getBottomAln());
   }

   if (!nucaln.empty()) {
      pair<string, string> result = nucaln.getConsensus();
      first = result.first;
      second = result.second;
      if (first.length() > 24) {
         if (aln.getSeq1Length() == 24) {
            first = aln.getSequence1AsString();
            second="";
         }
         else if (aln.getSeq2Length() == 24) {
            first = aln.getSequence2AsString();
            second="";
         }
         else {
            first.clear();
            second.clear();
         }
      }
   }
   else {
      first.clear();
      second.clear();
   }
   return first.length();
}

namespace po = boost::program_options;

void usage() {
   cout << "stemloop bamfile\n"
        << "Optons -i input bam file name default is rawlib.basecaller.bam\n"
        << "       -o output result file\n"
        << "       -r reference file name. The default is a path to the file:\n"
        << "           /home/zhouke/work/nextgen/shrna/CP0001_20130801_reference.csv\n"
        << "       --identity identity cutoff default 0.72\n"
        << "       --for-solexanome or\n"
        << "           --so  this will generate input for Solxanome program\n"
        << "           without doing the stemloop algorithm\n"
        << "       --help or ? will generate this help message\n";
   exit(1);
}

/**
 * Breaking up the main program so that we can work on different kind of input
 * files. This will reduce performance, otherwise, we have to do copy and
 * pasting.
 * We will keep the bam version efficient and the fasta version slower.
 * @return the align status for debug.
 */
template<class T>
int processOneSequence(DNA &rawseq, Dynaln<T> &alignl, Dynaln<T> &alignr, Dynaln<T> &aligns,
      ostream &ous, float idencut, map<pair<string,string>, int> &doubleHit,
      map<string, int> &categoryCount, map<string, int> &refcount, map<int, int> &lengthCount,
      map<int, int> &partialCount, map<QCInfo, int> &qc) 
{
   int insertB, insertE;
   int alnstatus = 0;
   QCInfo qci;
   if (alignWith(alignl, rawseq, ous, idencut)) {
      ++alnstatus;
      qci.setLeftIdentity(round2(alignl.getIdentity()));
      insertB = alignl.bottomEndIndex() + 1;
   }
   if (alignWith(alignr, rawseq, ous, idencut)) {
      alnstatus += 2;
      qci.setRightIdentity(round2(alignr.getIdentity()));
      insertE = alignr.bottomBeginIndex() - 1;
   }
   //cerr << "processOneSequence: alnstatus " << alnstatus << endl;

   if (alnstatus == 3) { // align to both left and right primer 
      ++categoryCount["alnboth"];
      ++lengthCount[findStemLoop(rawseq, insertB, insertE, ous, aligns, refcount, qci, doubleHit)];
   }
   else if (alnstatus == 1) { // align to left primer only
      ++categoryCount["alnleft"];
      ++partialCount[findStemLoopRight(rawseq, insertB, ous, aligns, refcount, qci, doubleHit)];
   }
   else if (alnstatus == 2) { // align to right primer only, very rare
      ++categoryCount["alnright"];
      findStemLoopLeft(rawseq, insertE, ous, aligns, refcount, qci, doubleHit);
   }
   else { // junk
      ous << "junk sequence not matchin either flank\n";
      ++categoryCount["junk"];
   }
   ++qc[qci];
   return alnstatus;
}

void analyzeResult(const map<pair<string,string>, int> &doubleHit,
      const map<string, int> &categoryCount,
      map<string, int> &refcount, const map<int, int> &lengthCount,
      const map<int, int> &partialCount, const map<QCInfo, int> &qc);

int tallyStemLoopFasta(const DNA &left, const DNA &right, const string &refile, 
      const string &inputfile, const string &outfile, float identityCut);


int main(int argc, char* argv[]) {
   //string inputfile = "/home/zhouke/work/nextgen/hnrna/shRNA_2_2_8ul_044/rawlib.basecaller.bam";
   string inputfile = "rawlib.basecaller.bam"; // run in the present directory
   string outfile, confile("stemloop.cfg");
   string referenceFile = "/home/zhouke/work/nextgen/shrna/CP0001_20130801_reference.csv";
   //int lengthCut=15;
   //float identityCut=0.72;
   bool forSolexanome = false;

   po::positional_options_description pd;
   po::options_description common("Common");
   po::options_description cmdline("Command options");
   po::options_description config("Configuration");

   po::variables_map vrbm;
   try {
      common.add_options()
         ("infile,i", po::value<string>(&inputfile)->default_value("rawlib.basecaller.bam"),
                   "input bam file containg raw reads")
         ("outfile,o", po::value<string>(&outfile),
                   "output result file containing the exact algorithm decisions for each sequence. Should be eliminated when the algorithm matures.")
         ("identity,y", po::value<float>()->default_value(0.72), 
                   "identity cutoff for matching to the flaking region sequences")
         ("left,l", po::value<string>(), "flanking sequence 5' to the hairpin")
         ("right,r", po::value<string>(), "flanking sequence 3' to the hairpin")
         ("format,m", po::value<string>()->default_value("bam"), "input sequence format, so far bam, and fasta. May add fastq in the future.")
         ("reference,f", po::value<string>(&referenceFile)->default_value("/home/zhouke/work/nextgen/shrna/CP0001_20130801_reference.csv"), 
                   "reference file for known hairpins. Used for final tabulation.");
      pd.add("infile", 1).add("outfile", 1);
      cmdline.add_options()
         ("help", "produce help message")
         ("conf,c", po::value<string>(&confile)->default_value("stemloop.cfg"), "Stemloop configuration file");
         ("solexanome,s", po::bool_switch(&forSolexanome), "generate input for Soxanome program without running the stemloop algorithm");
      cmdline.add(common);
      config.add(common);

      po::store(po::command_line_parser(argc, argv).options(cmdline).positional(pd).run(), vrbm);
      po::store(po::parse_config_file<char>(confile.c_str(), config), vrbm);
      po::notify(vrbm);
   }
   catch (exception &ex) {
      cerr << cmdline << endl << ex.what() << endl
         << "Failed the command line parsing\n";
      return 1;
   }
/*
   int i = 1;
   while (i < argc) {
      if (string(argv[i]) == "-i") inputfile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      //else if (string(argv[i]) == "-l") lengthCut = atoi(argv[++i]);
      else if (string(argv[i]) == "-r") referenceFile = argv[++i];
      else if (string(argv[i]) == "--identity") identityCut = atof(argv[++i]);
      else if (string(argv[i]) == "--help" 
            || argv[i][0] == '?') usage();
      else if (string(argv[i]) == "--for-solexanome"
            || string(argv[i]) == "--so") forSolexanome = true;
      else {
         inputfile = argv[i];
      }
      ++i;
   }
   //cout << "output file name: " << outfile << endl;
*/
   if (vrbm.count("help") || inputfile.empty()) {
      cerr << cmdline << endl;
      return 1;
   }
   if (outfile.empty()) {
      string::size_type s = inputfile.find('.');
      if (s != string::npos) {
         outfile = inputfile.substr(0,s) + ".slp";
      }
      else outfile = inputfile + ".slp";
   }


   //DNA left("left", "CCATCTCATCCCTGCGTGTCTCCGACTCAGCCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGAC");
   // for performance, using shorter oligos
   //DNA left("left", "CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATCGATGTAAAAGGAC");
   // without the Barcode + 17 nt spacer
   //                                             *
   //DNA left("le   "CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATC"); // 41 nt
   //                1   5    10   15   20   25   30   35   40
   //DNA left("left", "CCGTAACTTGAAAGTATTTCGATTTCTTGGCTTTATATATC"); // 41 nt
   //DNA right("right", "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGATGCTCTTCCGATCTTGTGGATGAATACTGCCATTTGTCTC");
   //DNA right("right", "GAGACAAATGGCAGTATTCATCCACAAGATCGGAAGAGC");
   //                  GCTCTTCCGATCTTGTG|GATGAATACTGCCATTTGTCTC GAGGTCGAGAATTCAAAAA 
   //                        0    6                22
   //                        42   47      55       64
         //                  |BC  | 17 nt          |
   //const string leftFlank = "GATGTAAAAGGACGAAACACCGG";
   // running solxanome barcode start 42, stemloop start 
   //const string rightFlank =   "TTTTTGAATTCTCGACCTC";


   DNA leftSide("leftside", vrbm["left"].as<string>());
   //leftSide.setName("leftside");
   //DNA rightSide = DNA(rightFlank) + right;
   DNA rightSide("rightside", vrbm["right"].as<string>());
   //rightSide.setName("rightSide");
   if (forSolexanome) {
      generateTrimmedFastaFile(leftSide, rightSide, inputfile, 
            leftSide.length(), rightSide.length(), vrbm["identity"].as<float>());
     return 0; // for testing.
   }
   if (vrbm["format"].as<string>() == "bam") {
      if (tallyStemLoop(leftSide, rightSide, referenceFile, inputfile, outfile,
            vrbm["identity"].as<float>()) != 0) {
         cerr << cmdline << endl;
         return 1;
      }
   }
   else if (vrbm["format"].as<string>() == "fasta") {
      cout << "using fasta file: " << inputfile << " as input\n";
      if (tallyStemLoopFasta(leftSide, rightSide, referenceFile, inputfile, outfile,
            vrbm["identity"].as<float>()) != 0) {
         cerr << cmdline << endl;
         return 1;
      }
   }

   return 0;
}

/**
 * Only use left and right context sequence to pick out stem loop region.
 * No more two stage process; this simplifies things. All the work is done at
 * the level of manipulating the context sequence.
 * @param inputfile is the bam input file containing all hairpin reads with
 *        flanking left and right sequences. The raw reads could be junk,
 *        primer dimer that are only the left and right flanks, or containing
 *        either left or righ flank only.
 * @return 0 for success and any other value for failure.
 */
int tallyStemLoop(const DNA &left, const DNA &right, const string &refile, 
      const string &inputfile, const string &outfile, float identityCut) 
{
   map<string, int> refcount;
   readReference(refile, refcount);
   
   BamTools::BamReader reader;
   if (!reader.Open(inputfile)) {
      cerr << "Failed to open " << inputfile << endl;
      return 1;
   }
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
      return 1;
   }

   SimpleScoreMethod scoreMethod(10, -9, -11, -3);
   // use three aligners
   Dynaln<SimpleScoreMethod> alignl(scoreMethod);
   alignl.setSeq1(left);
   Dynaln<SimpleScoreMethod> alignr(scoreMethod);
   alignr.setSeq1(right);
   Dynaln<SimpleScoreMethod> aligns(scoreMethod);

   BamTools::BamAlignment aln;
   int i = 0;
   int every=1000;
   map<int, int> lengthCount;
   map<int, int> partialCount;
   map<QCInfo, int> qc;
   map<string, int> categoryCount;
   map<pair<string, string>, int> doubleHit;
   int insertB, insertE, alnstatus;

   while (reader.GetNextAlignment(aln)) {
      ++i;
      if (i % every == 0) {
         cerr << "working on " << i << "th sequence\n";
         every *= 2;
      }
      if (aln.Length > left.length()) {
         ouf << string(54, '=') << endl << aln.QueryBases << endl;
         // 0 for no match, 1 left, 2 right only, 3 for both
         alnstatus = 0;
         DNA rawseq(aln.Name, aln.QueryBases);
         QCInfo qci;
         try {
            if (alignWith(alignl, rawseq, ouf, identityCut)) {
               ++alnstatus;
               qci.setLeftIdentity(round2(alignl.getIdentity()));
               insertB = alignl.bottomEndIndex() + 1;
            }
            if (alignWith(alignr, rawseq, ouf, identityCut)) {
               alnstatus += 2;
               qci.setRightIdentity(round2(alignr.getIdentity()));
               insertE = alignr.bottomBeginIndex() - 1;
            }

            if (alnstatus == 3) { // align to both left and right primer 
               ++categoryCount["alnboth"];
               ++lengthCount[findStemLoop(rawseq, insertB, insertE, ouf, aligns, refcount, qci, doubleHit)];
            }
            else if (alnstatus == 1) { // align to left primer only
               ++categoryCount["alnleft"];
               ++partialCount[findStemLoopRight(rawseq, insertB, ouf, aligns, refcount, qci, doubleHit)];
            }
            else if (alnstatus == 2) { // align to right primer only, very rare
               ++categoryCount["alnright"];
               findStemLoopLeft(rawseq, insertE, ouf, aligns, refcount, qci, doubleHit);
            }
            else { // junk
               ouf << "junk sequence not matchin either flank\n";
               ++categoryCount["junk"];
            }
         }
         catch (exception &ex) {
            cerr << "failed tabulation at " << i 
               << " " << rawseq << ", stopping the program " 
               << " alnstatus " << alnstatus << " " << ex.what() << endl
               << "production version not stopping, just skipping this sequence\n";
            //return 1;
         }
         ++qc[qci];
      }
      else {
         //ouf << "too short, discarded\n";
         ++categoryCount["tooshort"];
      }
   }
   ouf.close();
   cout << "Detailed algorithm action for programmer written to " << outfile << endl;
   analyzeResult(doubleHit, categoryCount, refcount, lengthCount, partialCount, qc);
   cerr << "Done.\n";

   return 0;
}


void analyzeResult(const map<pair<string,string>, int> &doubleHit,
      const map<string, int> &categoryCount,
      map<string, int> &refcount, const map<int, int> &lengthCount,
      const map<int, int> &partialCount, const map<QCInfo, int> &qc) {
   //////// result output section ////////////
   // output double match, top/bottom of the stem differ and both match to
   // library, this should be small
   ofstream ouf;
   ouf.open("refDoubleHitcount.tab");
   ouf << "topStemSeq\tbottomStemSeq\tcount\n";
   map<pair<string, string>, int>::const_iterator dit = doubleHit.begin();
   int match2Sum = 0;
   while (dit != doubleHit.end()) {
      ouf << dit->first.first << '\t' << dit->first.second << '\t' << dit->second << endl;
      match2Sum += dit->second;
      ++dit;
   }
   ouf.close();

   // summary file
   ouf.open("stemloop.summary");
   printMap(categoryCount, ouf);
   int total = 0;
   map<string, int>::const_iterator mit = categoryCount.begin();
   while (mit != categoryCount.end()) {
      total += mit->second;
      ++mit;
   }

   ouf  << "total sequences\t" << total << endl
      << "match2\t" << doubleHit.size() << endl
      << "match2Sum\t" << match2Sum << endl;

   // matched library seq, positive number, 
   // and novel input (not in library) count (negative number)
   MismatchMapping mm(refcount);
   mm.printForwardToFile("refcount.tab");
   //mm.summarize();
   mm.printSummary(ouf);

   // Considering those sequences that missed the library by one nucleotide as
   // true hit; thus the statistics will look a little bit better.
   mm.run();
   mm.printMappingToFile("refcount_miss1.tab");
   mm.printFinalSummary(ouf);
   ouf.close();

   // number of hits => count, negative number is for non-library sequences.
   ouf.open("numhits_distribution.tab");
   ouf << "numhits\tcount\n";
   mm.printDistribution(ouf);
   ouf.close();

   // count from stemloop with both flanking regions
   ouf.open("insertlencomplete.tab");
   //ouf << "Statistics of insert length\nInsert_length\tCount\n";
   ouf << "length\tcount\n";
   printMap(lengthCount, ouf);
   ouf.close();

   ouf.open("insertlenpartial.tab");
   ouf << "length\tcount\n";
   printMap(partialCount, ouf);
   ouf.close();

   // identity length multiple dimention count
   ouf.open("QCInfo.tab");
   ouf << "leftIden\trightIden\tstemIden\tstemLength\tcount\n";
   map<QCInfo, int>::const_iterator qit = qc.begin();
   map<QCInfo, int>::const_iterator qq;
   ofstream ouf3d("QCInfo_3d.tab");

   qq = qit;
   int sum = 0;
   while (qit != qc.end()) {
      ouf << qit->first << '\t' << qit->second << endl;
      if (qq->first.sameIdentity(qit->first)) {
         sum += qit->second;
      }
      else {
         qq->first.printIdentity(ouf3d) << '\t' << sum << endl;
         sum = qit->second;
         qq = qit;
      }
      ++qit;
   }
   qq->first.printIdentity(ouf3d) << '\t' << sum << endl;
   ouf.close();
   ouf3d.close();

   cerr << "Analysis Done.\n";
}

/**
 * Use fasta instead of bam as input file.
 */
int tallyStemLoopFasta(const DNA &left, const DNA &right, const string &refile, 
      const string &inputfile, const string &outfile, float identityCut) 
{
   map<string, int> refcount;
   readReference(refile, refcount);
   
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
      usage();
      return 1;
   }

   SimpleScoreMethod scoreMethod(10, -9, -11, -3);
   // use three aligners
   Dynaln<SimpleScoreMethod> alignl(scoreMethod);
   alignl.setSeq1(left);
   Dynaln<SimpleScoreMethod> alignr(scoreMethod);
   alignr.setSeq1(right);
   Dynaln<SimpleScoreMethod> aligns(scoreMethod);

   int i = 0;
   //int every=1000;
   int every=1;
   map<int, int> lengthCount;
   map<int, int> partialCount;
   map<QCInfo, int> qc;
   map<string, int> categoryCount;
   map<pair<string, string>, int> doubleHit;

   ifstream inf(inputfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open fasta input file: " << inputfile << endl;
      return 1;
   }

   DNA rawseq;
   string header;
   int alnstatus;

   while (rawseq.read(inf, header)) {
      ++i;
      if (i % every == 0) {
         cerr << "working on " << i << "th sequence " << rawseq << "\n";
         every *= 2;
      }
      if (rawseq.length() > left.length()) {
         ouf << string(54, '=') << endl << rawseq << endl;
         try {
            alnstatus = processOneSequence(rawseq, alignl, alignr, aligns, ouf, identityCut, 
                  doubleHit, categoryCount,refcount, lengthCount, partialCount, qc);
         }
         catch (exception &ex) {
            cerr << "failed tabulation at " << i 
               << " " << rawseq << " alnstatus: " << alnstatus 
               << ", stopping the program " << ex.what() << endl;
            return 1;
         }
      }
      else {
         ++categoryCount["tooshort"];
      }
   }
   inf.close();
   ouf.close();
   cout << "Detailed algorithm action for programmer written to " << outfile << endl;
   analyzeResult(doubleHit, categoryCount, refcount, lengthCount, partialCount, qc);
   cerr << "Done.\n";

   return 0;
}


/**
 * @return number of mapping in the reference library. -1 for empty input.
 *    0 for no hit, 1 for one hit, and 2 for two hits.
 * @param rfc is the REFSEQ -> count map to record the number of hits by a
 *        particular reference sequence.
 * @first most likely candidate. Should be 24 or longer
 */
int updateRefcount(map<string, int> &rfc, const string &first, const string &second,
      ostream &ouf) 
{
   int foundCount=-1;
   bool match1=false;
   bool match2=false;
   if (!first.empty()) {
      foundCount = 0;
      map<string, int>::iterator mit1, mit2;
      ouf << "stem sequence:\n" << first << endl;
      string stem = first.substr(0,21);
      string stem2;
      if ((mit1=rfc.find(stem)) != rfc.end()) {
         if (mit1->second >= 0) { // library sequence
            match1=true;
         }
         else {
            --(mit1->second);
         }
      }
      else {
         rfc[stem] = -1;
      }
      if (!second.empty()) {
         stem2 = second.substr(0,21);
         if (stem != stem2) {
            ouf << " another candidate:\n" << second << endl;
            if ((mit2=rfc.find(stem2)) != rfc.end()) {
               if (mit2->second >= 0) {
                  match2=true;
               }
               else {
                  --(mit2->second);
               }
            }
            else {
               rfc[stem2] = -1;
            }
         }

      }
      if (match1 && match2 ) {
         foundCount = 2;
         ouf << "found both candidates, not using either one\n";
         if (!boost::ends_with(first, "CTC")) {
            ++(mit2->second);
            foundCount = 1;
            ouf << "discarding " << first << " not ending with CTC\n";
         }
         else if (!boost::ends_with(second, "CTC")) {
            ++(mit1->second);
            foundCount = 1;
            ouf << "discarding " << second << " not ending with CTC\n";
         }
      }
      else if (match1) {
         foundCount = 1;
         ++(mit1->second);
         ouf << "found exact hit\n";
      }
      else if (match2) {
         foundCount = 1;
         ++(mit2->second);
         ouf << "found exact hit\n";
      }
      else {
         foundCount = 0;
         ouf << "Did not find in reflib!\n";
      }
   }
   return foundCount;
   // else both are empty, first got priority
}

/**
 * @param lengthcut is the length of left flank + right flank without insert
 */
void generateTrimmedFastaFile(const DNA &left, const DNA &right, const string &bamfile,
      const int leftFlankLength, const int rightFlankLength, const float identityCut) 
{
   BamTools::BamReader reader;
   if (!reader.Open(bamfile)) {
      cerr << "Failed to open " << bamfile << endl;
      return;
   }
   static const int lencut = 45;
   /*
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << endl;
      return;
   }
   */
   ofstream junkout("junk.fas");
   ofstream shortout("short.fas");
   ofstream partialout("partial.fas");
   ofstream dimerout("primerdimer.fas");
   ofstream completeout("complete.fas");
   ofstream complfasq("complete.fastq");
   ofstream partifasq("partial.fastq");

   //right.revcomp();
   //Matrix scmatrix;
   //scmatrix.setPath("/home/zhouke/src/proj/seqaln/matrix");
   //scmatrix.setMatrix("NUC.4.4");
   
   NucleicScoreMethod scoreMethod;
   Dynaln<NucleicScoreMethod> alignl(scoreMethod);
   //alignl.setMatrix(scmatrix);
   alignl.setSeq1(left);
   Dynaln<NucleicScoreMethod> alignr(scoreMethod);
   //alignr.setMatrix(scmatrix);
   alignr.setSeq1(right);

   BamTools::BamAlignment aln;
   int i = 0;
   int every=1000;
   int alnboth = 0;
   int alnleft = 0;
   int alnright = 0;
   int tooshort = 0;
   int junk = 0;

   map<int, int> lengthCount;
   map<int, int> partialCount;
   //map<float, int> leftIdentity;
   //map<float, int> rightIdentity;

   map<pair<float, float>, int> primerMatch;

   //int leftqual[3] = { 0, 0, 0 };
   //int middlequal[2] = { 0, 0 };
   //int rightqual[3] = { 0, 0, 0 }; // good, bad, diff by < 2
   //int middlegood = 0;
   //int middlebad = 0;
   int start, len;
   while (reader.GetNextAlignment(aln)) {
      ++i;
      if (i % every == 0) {
         cerr << "working on " << i << "th sequence "
            << alnboth << " have both ends match\n";
         every *= 2;
      }
      DNA rawseq(aln.Name, aln.QueryBases);
      if (aln.Length > left.length()-2) { // match left
         //ouf << string(54, '=') << endl 
         //   << "Original sequence\n" << endl;
         // 0 for no match, 1 left, 2 right only, 3 for both
         int alnstatus = 0;
         pair<float,float> twoiden = make_pair(0.0f,0.0f);
         //if (alignWith(alignl, rawseq, ouf, identityCut)) {
         if (alignWithNoOutput(alignl, rawseq, identityCut)) {
            ++alnstatus;
            twoiden.first=round2(alignl.getIdentity());
         }
         //if (alignWith(alignr, rawseq, ouf, identityCut)) {
         if (alignWithNoOutput(alignr, rawseq, identityCut)) {
            alnstatus += 2;
            twoiden.second=round2(alignr.getIdentity());
         }
         ++primerMatch[twoiden];

         if (alnstatus == 3) { // align to both left and right primer 
            ++alnboth; 
            len = alignr.bottomBeginIndex() - alignl.bottomEndIndex() - 1;
            ++lengthCount[len];
            try {
               if (len < 1) { //cerr << "primer dimer, no insert\n";
                  dimerout << rawseq << endl;
                  continue;
               }
               if (len < leftFlankLength + rightFlankLength + lencut) {
                  shortout << rawseq << endl;
               }
               else { // worth to be processed
                  start = alignl.bottomEndIndex()+1;
                  DNA insert = rawseq.subsequenceWithName(start, len);
                  insert.appendTitle("length=" + itos(len));
                  insert.appendTitle("left_identity=" + toString(round2(alignl.getIdentity()))
                        + " right_identity=" + toString(round2(alignr.getIdentity())));
                  completeout << insert; // sequence will output end of line
                  complfasq << "@" << insert.getName() << " " << insert.getTitle() << endl
                     << insert.toString() << endl << "+\n"
                     << aln.Qualities.substr(start, len) << endl;

               }
            }
            catch (bioseqexception &e) {
               cerr << e.what() << endl;
               cerr << "bad operation on bioseq: " << rawseq << endl;
               return;
            }
         }
         else if (alnstatus == 1) { // align to left primer only
            ++alnleft; 
            len = rawseq.length() - alignl.bottomEndIndex() - 1;
            //ouf << "Align to left primer only\n";
            if (len > leftFlankLength + 5) {
               start = alignl.bottomEndIndex()+1;
               DNA tail = rawseq.subsequenceWithName(start);
               tail.appendTitle("Only 5'-primer match. length=" + itos(len));
               tail.appendTitle("left_identity=" + toString(round2(alignl.getIdentity())));
               partialout << tail;
               partifasq << "@" << tail.getName() << " " << tail.getTitle() << endl
                  << tail.toString() << endl << "+\n"
                  << aln.Qualities.substr(start) << endl;
            }
            else {
               shortout << rawseq;
            }
            ++partialCount[len];
         }
         else if (alnstatus == 2) { // align to right primer only, very rare
            ++alnright; 
            len = alignr.bottomBeginIndex();
            DNA head = rawseq.subsequenceWithName(0, len);
            head.appendTitle("Align to right primer only. length=" + itos(len));
            head.appendTitle("right_identity=" + toString(round2(alignr.getIdentity())));
            partialout << head;
            partifasq << "@" << head.getName() << " " << head.getTitle() << endl
               << head.toString() << endl << "+\n"
               << aln.Qualities.substr(0, len) << endl;
         }
         else { // does not alight with left or right above cutoff identity, junk
            junkout << rawseq;
            ++junk;
         }
      }
      else {
         //ouf << "too short, discarded\n";
         ++tooshort;
         junkout << rawseq;
      }
   }
   //ouf.close();
   junkout.close();
   shortout.close();
   partialout.close();
   dimerout.close();
   completeout.close();
   complfasq.close();
   partifasq.close();

   ofstream ouf;
   ouf.open("primermatch.summary");
   ouf << alnboth << "\tmatches both ends\n" 
      << alnleft << "\tmatches left ends only\n" 
      << alnright << "\tmatches right ends only\n" 
      << tooshort << "\ttooshort\n" 
      << junk << "\tjunk\n" 
      << i << "\ttotal sequences\n";
   ouf.close();

   ouf.open("insertlencomplete.tab");
   //ouf << "Statistics of insert length\nInsert_length\tCount\n";
   ouf << "length\tcount\n";
   printMap(lengthCount, ouf);
   ouf.close();

   ouf.open("insertlenpartial.tab");
   ouf << "length\tcount\n";
   printMap(partialCount, ouf);
   ouf.close();

   ouf.open("identityLeftRight.tab");
   ouf << "left_identity\tright_identity\tcount\n"; // header
   map<pair<float,float>, int>::const_iterator it = primerMatch.begin();
   while (it != primerMatch.end()) {
      ouf << it->first.first << '\t' << it->first.second << '\t' << it->second << endl;
      ++it;
   }
   ouf.close();
   cout << "Bam file processed, fasta file produced\n";
}

template<class T1, class T2> 
void printMap(const map<T1, T2> &m, ostream &ous) {
   typename map<T1, T2>::const_iterator it = m.begin();
   while (it != m.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}



void readReference(const string &file, map<string, int> &refcount) {
   ifstream inf(file.c_str());
   if (!inf) {
      cerr << "Failed to open input file: " << file << endl;
      exit(1);
   }
   string line;
   getline(inf, line);
   if (line[0] != 'A' || line[0] != 'C' || line[0] != 'G' || line[0] != 'T') {
      // read potentiall header
      getline(inf, line);
   }
   string seq;
   while (!inf.eof()) {
      seq = line.substr(0, line.find(','));
      //cout << seq << endl;
      if (seq.length() != 21) {
         cout << "reference incorrect: " << seq << endl;
         exit(1);
      }
      refcount[seq] = 0;
      getline(inf, line);
   }
   cout << refcount.size() << " references read from file: " << file << endl;
}      

int numdiff(const string &longstr, const string &shortstr) {
   int diff = 0;
   for (string::size_type i = 0; i<shortstr.length(); ++i) {
      if (longstr[i] != shortstr[i]) {
         ++diff;
      }
   }
   return diff;
}

int numdiffTail(const string &longstr, const string &shortstr) {
   int diff = 0;
   string::const_reverse_iterator it = longstr.rbegin();
   for (string::size_type i = shortstr.length() - 1; i > -1; --i) {
      if (*it != shortstr[i]) {
         ++diff;
      }
      ++it;
   }
   return diff;
}



