#include <iostream>
#include <cstring>
#include <fstream>
#include <boost/program_options.hpp>
#include <algorithm>
#include <ctime>

#include <bioseq.h>
#include <fastq.h>
#include <dynalnt.h>

using namespace std;
using namespace orpara;

ostream& outputCurrentTime(ostream &ous);
class NonrepeatBaseException : public exception {
   private:
      string message;

   public:
      NonrepeatBaseException(const string &msg) throw() 
         : message(msg) { }
      ~NonrepeatBaseException() throw() { }
      const char* what() const throw() {
         return message.c_str();
      }
};

/**
 * Simple representation.
 * Index, number of repeats, Base repeated.
 */
class RepeatParam {
   public:
      RepeatParam()
         : index(0), rep(0), base(' ') { }

      RepeatParam(int i, int r, char b)
         : index(i), rep(r), base(b) { }

      string getRepeatSequence() const {
         return string(rep, base);
      }

      int index;
      int rep;
      char base;
};

//using namespace boost::program_options;
namespace po = boost::program_options;
void printHelpHeader(ostream &ous);

bool matchSite(const Dynaln<SimpleScoreMethod> &aligner, const po::variables_map &vm);
int alignSite(Dynaln<SimpleScoreMethod> &aligner1, Dynaln<SimpleScoreMethod> &aligner2, Fastq &fsq, po::variables_map &parm, ostream &ous1, ostream &ous2);
/**
 * Will remove path prefix from the output
 * So the files will be in the current directory.
 */
string nameOutputFile(const string &infile);
int countRepeat(const Dynaln<SimpleScoreMethod> &aligner, const RepeatParam &reparam, ostream &ous, string &variant);
int numberOfRepeat(const string &alignedseq, const char B) throw(NonrepeatBaseException);
bool countingRepeat(const po::variables_map &vrbm);
bool validateRepeatParam(const po::variables_map &vrbm);
void transformParameters(const po::variables_map &vrbm, RepeatParam &fwd, RepeatParam &bwd);

/**
 * extract file name with extensions
 */
string extractFileName(const string &infile) {
   string::size_type i = infile.rfind('/');
   string fileName = infile.substr(i+1);
   return fileName;
}

string replaceExtensionWith(const string &infile, const string &tail) {
   string outfile = infile;
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      outfile = infile.substr(0,i);
   }
   outfile += tail;
   return outfile;
}

void countsub(const string &fqfile, const string &bait, const string &range,
      float identitycut);

/**
 * Looking for fastq sequences that match a given site: short sequence.
 */
int main(int argc, char* argv[]) {
   if (argc == 1 || (argc == 2 && !strcmp(argv[1], "--help"))) 
      printHelpHeader(cerr);
   //int i = 1;
   //string infile("bar4.fastq"), outfile("bar4cut.fastq");
   string infile, outfile1, outfile2, outfile, baitstr;
   //infile="R01sid.fastq"; // for debug
   //string baitstr="AGAAGCTTGCTTCTTTGCTGACGAGTGGCGG";
   bool doCountRepeat = false;
   po::options_description desc("Allowed options");
   po::positional_options_description pd;
   po::variables_map vrbm;
   //string cutsite="cccgggttgcgccttttccaaggcagccctgggtttgcgcagggacgcggct";
   //Above is for testing
   try {
      //cerr << "processing command line ...\n";
      desc.add_options()
         ("help", "produce help message")
         ("infile,i", po::value<string>(&infile)->default_value("raw.fastq"), "input fastq file")
         ("outfile1,f", po::value<string>(&outfile1), "file for alignments in the forward direction")
         ("outfile2,b", po::value<string>(&outfile2), "file for alignments in the backward direction")
         ("outfile,o", po::value<string>(&outfile), "file for storing matched fastq sequences")
         ("site,s", po::value<string>(&baitstr), "The ambiguous site sequence. A length of 35-45 works the best. A good example: CGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGAC")
         ("identity-cutoff,y", po::value<float>()->default_value(0.88), "The sequence identity cutoff bewtween site and target sequence above which target will be processed. A fraction (0-1].")
         ("alnlen-cutoff,l", po::value<int>()->default_value(16), "The length cutoff of the alignment between the cut-site and target sequence above which target will be processed.")
         ("cov-cutoff,c", po::value<float>()->default_value(0.7), "cutoff for the coverage of the site.")
         ("trim,t", "trim unmatched (to bait) part of target. Default not trimming.")
         ("repeat-start,p", po::value<int>(), "Repeat start in position in the site sequence, 1-based index")
         ("repeat-count,n", po::value<int>(), "Number of repeats")
         ("repeat-base,r", po::value<char>(), "Base in the tandem repeat,such as A, T")
         ("countsub,k", po::value<string>(), "Count the subsequece of the bait given by b-e")
         ;

      //po::positional_options_description pd;
      pd.add("infile", 1).add("outfile", 1);
      //po::store(po::parse_command_line(argc, argv, desc), vrbm);
      po::store(po::command_line_parser(argc, argv).options(desc).positional(pd).run(), vrbm);
      po::notify(vrbm);
      //cerr << "command line processed\n";
   }
   catch (exception &ex) {
      cerr << desc << endl;
      cerr << ex.what() << " Failed the command line processing with boost\n";
      return 1;
   }

   if (vrbm.count("help")) {
      cout << desc << endl;
      return 1;
   }
   if (baitstr.empty()) {
      cerr << desc << endl << "Positional parameters: infile outfile\n";
      cerr << "You must enter a cut site sequence with the -site option\n";
      return 1;
   }
   if (vrbm.count("countsub")) {
      //string posrange="9-18";
      //countsub(infile, baitstr, posrange, 0.85);
      countsub(infile, vrbm["site"].as<string>(), vrbm["countsub"].as<string>(),
            vrbm["identity-cutoff"].as<float>());
      return 0;
   }

   // repeate counting operation validation
   RepeatParam repeatParam, reverseRepeatParam;
   if (countingRepeat(vrbm)) {
      doCountRepeat = true;
      if (validateRepeatParam(vrbm)) {
         transformParameters(vrbm, repeatParam, reverseRepeatParam);
      }
      else {
         cerr << "Invalid repeat parameters given!\n";
      }
   }

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "failed to open " << infile << " for reading\n";
      return 1;
   }
   // these file naming should be more specific
   if (!infile.empty() && outfile1.empty()) {
      string fn = extractFileName(infile);
      outfile1 = replaceExtensionWith(fn, ".match");
      outfile2 = replaceExtensionWith(fn, "rc.match");
   }
   if (outfile.empty()) {
      outfile = nameOutputFile(infile);
   }

   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writting. Check for permission.\n";
      return 1;
   }
   ofstream ouf1(outfile1.c_str());
   if (ouf1.fail()) {
      cerr << "failed to open " << outfile1 << " for writting\n";
      return 1;
   }
   ofstream ouf2(outfile2.c_str());
   if (ouf2.fail()) {
      cerr << "failed to open " << outfile2 << " for writting\n";
      return 1;
   }
   ouf1 << "Sequence matching the site in forward direction\n\n";
   ouf2 << "Sequence matching the site in backward direction\n\n";

   Fastq fasq;
   //DNA cutseq("cutsite", boost::algorithm::to_upper_copy(vrbm["site"].as<string>()));
   DNA baitseq("bait", vrbm["site"].as<string>());
   baitseq.toUpperCase();
   DNA baitseqrc = baitseq.revcompCopy();

   //MatrixScoreMethod::setDefaultPath("/home/zhouke/src/proj/seqaln/matrix");
   //Matrix scoreMatrix("NUC.4.4", true); // nucleoic acid 4x4
   // give less penalty to gap because this is the type of error of iron
   // torrent that is more likely.
   SimpleScoreMethod scoreMatrix(7, -7, -5, -3); // identity matrix
   //scoreMatrix.show();

   Dynaln<SimpleScoreMethod> alignerForward(scoreMatrix);
   alignerForward.setSeq1(baitseq);

   Dynaln<SimpleScoreMethod> alignerBackward(scoreMatrix);
   alignerBackward.setSeq1(baitseqrc);

   int count[3]={0,0,0};
   int countbackward[2]={0,0};
   map<int, int> repeatStat;
   map<string, int> repeatVariant;

   try {
      int nr;
      int numread = 0;
      int every = 10;
      string rvar;
      while (fasq.read(inf)) {
         ++numread;
         if (numread % every == 0) {
            cerr << "working on " << numread << " read " << fasq.getName()
            << "  ";
            outputCurrentTime(cerr) << endl;
            every *= 2;
         }
         int alnType=alignSite(alignerForward, alignerBackward, fasq, vrbm, ouf1, ouf2);
         if (alnType == 2) { // no align
            ++count[2];
         }
         else if ((alnType & 4) == 0) { // forward align
            if (vrbm.count("trim")) {
               ouf << fasq.sub(alignerForward.bottomBeginIndex(), alignerForward.bottomEndIndex());
            }
            else ouf << fasq;

            ++count[alnType & 3];
            //cout << "working on forward strand ...\n";
            if (doCountRepeat) {
               nr = countRepeat(alignerForward, repeatParam, ouf1, rvar);
               if (!rvar.empty()) ++repeatVariant[rvar];
               ++repeatStat[nr];
            }
         }
         else { // backward align
            if (vrbm.count("trim")) {
               ouf << fasq.sub(alignerBackward.bottomBeginIndex(), alignerBackward.bottomEndIndex());
            }
            else ouf << fasq;

            ++countbackward[alnType & 3];
            if (doCountRepeat) {
               nr = countRepeat(alignerBackward, reverseRepeatParam, ouf2, rvar);
               if (!rvar.empty()) ++repeatVariant[reverseComplement(rvar)];
               ++repeatStat[nr];
            }
         }
      }
      cerr << numread << " sequences analyzed loop end at: ";
      outputCurrentTime(cerr) << endl;
   }
   catch (exception &e) {
      cerr << e.what() << endl << " failed to do fishing\n";
      return 1;
   }
   inf.close();
   ouf.close();
   ouf1.close();
   ouf2.close();

   cerr << "Final min max score:\n";
   cerr << fasq.getMinScore() << "  " << fasq.getMaxScore() << endl;
   cout << "result written to " << outfile1 << " and " << outfile2 << endl
      << "   " << count[2] << " no match \n"
      << "   Align forward: " << count[0] << " perfect match to site\n"
      << "                  " << count[1] << " match site with defects\n"
      << "   Align backward: " << countbackward[0] << " perfect match to site\n"
      << "                   " << countbackward[1] << " match site with defects\n";
   
   if (doCountRepeat) {
      int sum = 0;
      int n = 0;
      map<int, int>::const_iterator mit = repeatStat.begin();
      cout << "repeat\tcount\n";
      while (mit != repeatStat.end()) {
         cout << mit->first << '\t' << mit->second << endl;
         if (mit->first > 0) {
            n += mit->second;
            sum += (mit->first * mit->second);
         }
         ++mit;
      }
      cout << "average: " << (double)sum/n << " n: " 
         << n << " " << vrbm["repeat-base"].as<char>() << endl;
      map<string, int>::const_iterator it = repeatVariant.begin();
      cout << "Variant of the repeat:\n";
      while (it != repeatVariant.end()) {
         cout << it->first << '\t' << it->second << endl;
         ++it;
      }
   }
   //outputCurrentTime(cerr) << endl;
   cerr << "fishing done\n";

   return 0;
}

void displayProgress(int &n, int &e, const string &msg) {
   if (++n % e == 0) {
      cerr << "working on " << n << " " << msg << endl;
      e *= 2;
   }
}


/**
 * Given a fastq file, it counts the subsequence of the bait
 * @param range 1-based start and end position inside bait
 * @param bait is the short sequene to use
 */
void countsub(const string &fqfile, const string &bait, const string &range,
      float identitycut) {
   cout << "bait: " << bait << " " << bait.length() << endl << " range: " << range << endl;
   vector<int> bound = getAllInt(range);
   if (bound.size() != 2) {
      throw invalid_argument("sub range must be specified with 123-456 or 123,456");
   }
   int reverseB=bait.length()-bound[1]+1;
   int reverseE=bait.length()-bound[0]+1;
   cout << "reverse range: " << reverseB << "-" << reverseE << endl;
   DNA baitseq("bait", bait);
   baitseq.toUpperCase();
   DNA baitseqrc = baitseq.revcompCopy();
   SimpleScoreMethod scoreMatrix(5, -3, -13, -1);
   Dynaln<SimpleScoreMethod> alignerF(scoreMatrix);
   alignerF.setSeq1(baitseq);
   Dynaln<SimpleScoreMethod> alignerB(scoreMatrix);
   alignerB.setSeq1(baitseqrc);
   map<string, int> subcnt; // subsequence count
   map<string, int>::iterator mit;
   Fastq fasq;
   ifstream inf(fqfile);
   if (inf.fail()) {
      throw runtime_error("Failed to open fastq file: " + fqfile);
   }
   int numread=0;
   int every=5;
   outputCurrentTime(cerr) << endl;
   while (fasq.read(inf)) {
      DNA seq(fasq.getName(), fasq.getSequence());
      displayProgress(numread, every, fasq.getName());
      alignerF.setSeq2(seq); 
      alignerF.runlocal(); 
      string bottomStr;
      if (alignerF.getIdentity() > identitycut && alignerF.getCov1()> 0.99) {
         bottomStr = alignerF.getBottomAlnByTopPosition(bound[0], bound[1]);
      }
      else {
         alignerB.setSeq2(seq);
         alignerB.runlocal();
         if (alignerB.getIdentity()>identitycut && alignerB.getCov1()>0.99) {
            bottomStr = alignerB.getBottomAlnByTopPosition(reverseB, reverseE);
            reverseComplementInPlace(bottomStr);
         }
      }
      if (!bottomStr.empty()) {
         if ((mit=subcnt.find(bottomStr)) == subcnt.end()) {
            subcnt.insert(make_pair(bottomStr, 1));
         }
         else {
            ++(mit->second);
         }
      }
   }
   cerr << numread << " sequences analyzed. ";
   outputCurrentTime(cerr) << endl;
   cout << "fished target count\n";
   vector<pair<string,int> > tmp;
   tmp.reserve(subcnt.size());
   for (auto i=subcnt.begin(); i != subcnt.end(); ++i) {
      //cout << i->first << " " << i->second << endl;
      tmp.push_back(make_pair(i->first, i->second));
   }
   sort(tmp.begin(), tmp.end(), [](const pair<string,int> &p1, const pair<string,int> &p2)->bool{ return p1.second > p2.second; });
   for (size_t i=0; i<tmp.size(); ++i) {
      cout << tmp[i].first << " " << tmp[i].second << endl;
   }
}

ostream& outputCurrentTime(ostream &ous) {
   time_t now = time(0);
   tm* ltm = localtime(&now);
   ous << "time: " 
      << ltm->tm_hour << ":"
      << ltm->tm_min << ":"
      << ltm->tm_sec;
   return ous;
}

string nameOutputFile(const string &infile) {
   string outfile = extractFileName(infile);
   string::size_type i = outfile.rfind('.');
   if (i != string::npos) {
      outfile = outfile.substr(0,i);
   }
   outfile += "_fish.fastq";
   return outfile;
}


bool matchSite(const Dynaln<SimpleScoreMethod> &aligner, const po::variables_map &vm) {
   return aligner.getIdentity() > vm["identity-cutoff"].as<float>()  
         && aligner.getAlnlen() > vm["alnlen-cutoff"].as<int>()
         && aligner.getCov1() > vm["cov-cutoff"].as<float>();
}

/*
 * right 2 bits for align type 0 for perfect, 1 for non-perfect 2 for no align
 * the rigth 3rd bits is for direction, 0 for forward, 1 for backward
 * |D|T|T|
 *  0 0 0 perfect
 *  1 0 1 non-perfect
 *    1 0 no align
 */
int alignSite(Dynaln<SimpleScoreMethod> &aligner1, Dynaln<SimpleScoreMethod> &aligner2, Fastq &fsq, 
      po::variables_map &parm, ostream &ous1, ostream &ous2) 
{
   DNA rawseq(fsq.getName(), fsq.getSequence());
   aligner1.setSeq2(rawseq); 
   aligner2.setSeq2(rawseq); 
   aligner1.runlocal();
   aligner2.runlocal();

   // 0 for perfect, 1 for non-perfect
   int alignType=0;
   try {
      if (aligner1.isPerfectSeq1Align()) { // forward perfect 000
         fsq.write(ous1);
         aligner1.printAlign(ous1);
         ous1 << "Matched perfectly to site in forward direction.\n";
      }
      else if (aligner2.isPerfectSeq1Align()) { // backward perfect 100
         alignType |= 4; // backward direction using right 3rd bit
         fsq.write(ous2);
         aligner2.printAlign(ous2);
         ous2 << "Matched perfectly to site in backward direction\n";
      }
      else if (matchSite(aligner1, parm) || matchSite(aligner2, parm)) {
         // 001
         alignType = 1;
         if (aligner1.getScore() > aligner2.getScore()) {
            fsq.write(ous1);
            aligner1.printAlign(ous1);
         }
         else { // 101
            alignType |= 4;
            fsq.write(ous2);
            aligner2.printAlign(ous2);
         }
      }
      else { // 010 no align regardless of direction, the third (from right)
         // bit is not needed
         alignType = 2;
      } 
   }
   catch (exception &e) {
      cerr << e.what() << " failed to process the forard cut-site\n";
   }

   return alignType;
}

int numberOfRepeat(const string &alignedseq, const char B) throw(NonrepeatBaseException) {
   int count = 0;
   for (unsigned int i = 0; i<alignedseq.length(); ++i) {
      if (alignedseq[i] == B) ++count;
      else if (alignedseq[i] == '-') { }
      else {
         throw NonrepeatBaseException(alignedseq + " contains bases other than repeat base: " + B);
      }
   }
   return count;
}

/**
 * the bioseq in aligner is out of scope so many functions cannot be used
 * @param tops index of the repeat start in the top unaligned sequence.
 */
int countRepeat(const Dynaln<SimpleScoreMethod> &aligner, const RepeatParam &reparam, ostream &log, 
      string &variant) {
   string Tn = reparam.getRepeatSequence();
   int beginIndex = reparam.index - aligner.topBeginIndex();
   string top = aligner.getTopAln();
   string bottom = aligner.getBottomAln();
   //log << "top length " << top.length() << " " << beginIndex << endl;
   int i,j; // j count ungapped site sequence, i count aligned seq
   i = 0; j = 0;
   while (i < (int)top.size() && j < beginIndex) {
      if (top[i] == '-') { ++i; }
      else { ++i; ++j; }
   }
   int x, y; // outer boundary
   x = i;
   while (i < (int)top.size() && top[i] == '-') ++i; // walk over gap
   log << i - j << " gaps on the top before repeat\n";
   // i at the frist Repeat Base, such as A or T
   j = i + 1;
   // walk over the repeat that my have gap in the alighment
   // no more than rep bases of T or A
   int bb = 1; // count repeat base
   while (j < (int)top.size() && (top[j] == reparam.base || top[j] == '-')
         && bb <= reparam.rep) {
      ++j;
      if (top[j] == reparam.base) ++bb;
   }
   // now j is at the first non-repeat base or gapt char -
   --j;
   // walk back in case walked over the boundary
   while (j > i && top[j] == '-' && bottom[j] != reparam.base) --j;
   y = j; // j passed the repeat range

   // extend potential extension of repeat in target
   while (y < (int)top.size() && top[y] == '-' && bottom[y] == reparam.base) {
      ++y; 
   }
   int numrep = 0;
   /* for debug slows down performance
   if (top.substr(i, reparam.rep) != Tn) {
      log << "WARN: site index " << beginIndex
         << " top align starting at " << i << " "
         << top.substr(i, reparam.rep) << " is not " << Tn << endl;
      //exit(1);
   }
   */
   //log << "bottom align starting at " << i << " "
   //   << bottom.substr(i, reparam.rep) << endl;
   if (y-x + 1 > reparam.rep) {
      log << "expanded target match: " << top.substr(x, y-x+1) 
         << " | " << bottom.substr(x, y-x+1) << endl;
   }
   //else {
   //   log << "target region: " << bottom.substr(x, y-x+1) << endl;
   //}

   try {
      numrep = numberOfRepeat(bottom.substr(x, y-x+1), reparam.base);
      variant.clear();
   }
   catch (NonrepeatBaseException &ne) {
      log << "WARN: " << ne.what() << endl; // just return 0
      variant = bottom.substr(x, y-x+1);
      remove(variant.begin(), variant.end(), '-');
      //trim_if(variant, is_any_of('-'));
   }
   //log << numrep << " repeats " << bottom.substr(x, y-x+1) << endl;
   return numrep;
}


void printHelpHeader(ostream &ous) {
   ous << "Fishing - is a program to pick sequences from a given input file\n"
      << "guided by a bait sequence (-s option)\n"
      << "There will be two match files that contains the details of the \n"
      << "alignment of bait to the target sequence.  The main output file\n"
      << "contains all the sequences pull out from the input file. \n"
      << "The output file (opton -o) has extension *.fastq\n\n";
}

void transformParameters(const po::variables_map &vrbm, RepeatParam &fwd, RepeatParam &bwd) {
   int idx = vrbm["repeat-start"].as<int>()-1;
   int n = vrbm["repeat-count"].as<int>();
   char B = vrbm["repeat-base"].as<char>();

   fwd.index = idx;
   fwd.rep = n; 
   fwd.base = B;

   bwd.index = (int)vrbm["site"].as<string>().length() - idx - n;
   bwd.rep = n;
   bwd.base = DNA::toRevcompBase(B);
}

bool countingRepeat(const po::variables_map &vrbm) {
   if (vrbm.count("repeat-start") && vrbm.count("repeat-count")
         && vrbm.count("repeat-base")) {
      cout << "I will count repeat\n";
      return true;
   }
   return false;
}

bool validateRepeatParam(const po::variables_map &vrbm) {
   string bait = vrbm["site"].as<string>();
   string repInBait = bait.substr(vrbm["repeat-start"].as<int>()-1, vrbm["repeat-count"].as<int>());
   string repMade = string(vrbm["repeat-count"].as<int>(), vrbm["repeat-base"].as<char>());
   if (repInBait != repMade) {
      cout << "repeat in the bait: " << repInBait << " not the same as stated: "
         << repMade << endl;
      return false;
   }
   return true;
}
