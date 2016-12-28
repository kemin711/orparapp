#include <iostream>
#include <cstring>
#include <fstream>
#include <boost/program_options.hpp>

#include <bioseq.h>
#include <dynalnt.h>
#include <fastq.h>

//#include <boost/program_options/options_descripton.hpp>
//#include <boost/program_options/parsers.hpp>
//#include <boost/program_options/variables_map.hpp>
//#include <boost/algorithm/string.hpp>
//#define DEBUG

using namespace std;
using namespace orpara;
//using namespace boost::program_options;
namespace po = boost::program_options;

string getFileStem(const string &fileName);
//bool matchSite(Dynaln &aligner, Fastq &fsq, const Progparam &parm, int direction);
bool isBadAlign(const Dynaln<SimpleScoreMethod> &aligner, po::variables_map &vm);
bool matchSite(const Dynaln<SimpleScoreMethod> &aligner, const po::variables_map &vm);

int bestCut(Dynaln<SimpleScoreMethod> &alnForward1, Dynaln<SimpleScoreMethod> &alnForward2, 
      Dynaln<SimpleScoreMethod> &alnBackward1, Dynaln<SimpleScoreMethod> &alnBackward2,
      Fastq &fsq, Fastq &rightq, const po::variables_map &parm, ostream &ous);
int cutAtTail(Dynaln<SimpleScoreMethod> *aligner, Fastq &fsq, Fastq &rightq, 
      const po::variables_map &parm, ostream &ous);
int cutAtHead(Dynaln<SimpleScoreMethod> *aligner, Fastq &fsq, Fastq &rightq, 
      const po::variables_map &parm, ostream &ous);

/**
 * Cut raw reads sequence from next generation sequences according to
 * a site: a DNA sequence from 30-50 nucleotides.
 * Cut sites: head, tail
 * head is the first n-nt and tail is the n-nt of the tail end of the
 * linearized plasmid.
 */
int main(int argc, char* argv[]) {
   //int i = 1;
   //string infile("bar4.fastq"), outfile("bar4cut.fastq");
   string infile, outfile, confile("linearize.cfg");

   po::options_description common("Common");
   po::options_description cmdline("Command options");
   po::options_description config("Configuration");

   po::variables_map vrbm;
   //string cutsite="cccgggttgcgccttttccaaggcagccctgggtttgcgcagggacgcggct";
   //Above is for testing
   try {
      common.add_options()
         ("head,h", po::value<string>(), "The CUT site sequence from beginning end. Will cut on the left.  It is best that the sequence be devoid of single nucleotide repeats. A length of 35-45 works the best. Too long slows down and discards too many raw reads. Too short will be non-specific. A good example: CGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGAC")
         ("tail,t", po::value<string>(), "The Cut site sequence from the tail end. Will cut on the right.  Same requirements from the site1 option.  5'=head-site=>--circular plasmid----=tail-site=>3'")
         ("identity-cutoff,y", po::value<float>()->default_value(0.84), "The sequence identity bewtween site and target sequence above which target will be processed.")
         ("alnlen-cutoff,l", po::value<int>()->default_value(17), "The length of the alignment between the cut-site and target sequence above which target will be processed.")
         ("iden-discard,d", po::value<float>()->default_value(0.92), "Identity between site and target sequence below which the sequence from the left margin to the end of target sequence will be discarded.")
         ("gap-discard,g", po::value<int>()->default_value(2), "Number of gaps in site or target above which the sequence from the left margin to the end of target sequence will be discarded.")
         ("shortest,x", po::value<unsigned int>()->default_value(18), "The shortest fastq sequence to keep after the cut operation.");

      //cerr << "processing command line ...\n";
      cmdline.add_options()
         ("help", "produce help message")
         ("infile,i", po::value<string>(&infile), "input fastq file")
         ("outfile,o", po::value<string>(&outfile), "output cut fastq file")
         ("conf,c", po::value<string>(&confile)->default_value("linearize.cfg"), "The configuration file.");
      cmdline.add(common);

      // should use a sub option to eliminate the duplication
      config.add(common);

      po::positional_options_description pd;
      pd.add("infile", 1).add("outfile", 1);

      //po::store(po::parse_command_line(argc, argv, desc), vrbm);
      // for positional parameters, you need a parser class
      po::store(po::command_line_parser(argc, argv).options(cmdline).positional(pd).run(), vrbm);
      cout << "processing config file: " << confile << "\n";
      po::store(po::parse_config_file<char>(confile.c_str(), config), vrbm);
      po::notify(vrbm);
      //cerr << "command line processed\n";
   }
   catch (exception &ex) {
      cerr << ex.what() << " Failed the command line processing with boost\n";
      cerr << cmdline << endl;
      return 1;
   }

   if (vrbm.count("help")) {
      cout << cmdline << endl;
      return 1;
   }
   if (!vrbm.count("head") || !vrbm.count("tail")) {
      cout << cmdline << endl;
      cerr << "You must enter both head and tail cut site sequences with the -h and -t option\n";
      return 1;
   }

   if (infile.empty()) {
      cerr << cmdline << endl << "You must provide input file\n";
      return 1;
   }
   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "failed to open " << infile << " for reading\n";
      return 1;
   }
   if (outfile.empty()) {
      string::size_type i = infile.find('.');
      if (i == string::npos) outfile = infile;
      else outfile = infile.substr(0,i);
      outfile += "cut.fastq";
   }
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "failed to open " << outfile << " for writting\n";
      return 1;
   }

   Fastq fasq;
   //DNA cutseq("cutsite", boost::algorithm::to_upper_copy(vrbm["site"].as<string>()));
   cout << "using head cut site: " << vrbm["head"].as<string>() << endl
      << "using tail cut site: " << vrbm["tail"].as<string>() << endl;
   DNA headseq("headcutsite", vrbm["head"].as<string>());
   headseq.toUpperCase();
   DNA headseqrc = headseq.revcompCopy();
   cout << "reverse complement of head cut site:\n" << headseqrc << endl;
   // tail site is the opposite
   // reverse complement is the topcut site
   DNA tailseq("tailcutsite", vrbm["tail"].as<string>());
   tailseq.toUpperCase();
   DNA tailseqrc = tailseq.revcompCopy();
   cout << "tail " << tailseq << endl << "tailrc " << tailseqrc << endl;
   // --tail->|--head->
   // <-t-rc--|<-h-rc-- in terms of cutting, head && t-rc are in the same
   // configuration.

   // penalize more on gaps
   SimpleScoreMethod scmethod(5, -5, -8, -4);

   // for head 
   Dynaln<SimpleScoreMethod> alignerForwardh(scmethod);
   alignerForwardh.setSeq1(headseq);
   // for tail reverse complement, tail-rc
   Dynaln<SimpleScoreMethod> alignerForwardt(scmethod);
   alignerForwardt.setSeq1(tailseqrc);

   // head reverse complement head-rc
   Dynaln<SimpleScoreMethod> alignerBackwardh(scmethod);
   alignerBackwardh.setSeq1(headseqrc);

   // tail
   Dynaln<SimpleScoreMethod> alignerBackwardt(scmethod);
   alignerBackwardt.setSeq1(tailseq);

   int count_cut = 0;
   //int count_discard_partial = 0;
   int count_discard = 0;
   int count_nocut = 0;
   string logfile = getFileStem(infile) + "_linearize.log";
   ofstream LOG(logfile.c_str());
   int every = 5;

   try {
      int i = 0; 
      while (fasq.read(inf)) {
         ++i;
         if (i % every == 0) {
            cout << "processing " << i << "th fastq\n";
            every *= 2;
         }
         // try dynamic string match algorithm
         int cutType;
         Fastq rightq; // this must start from an empty sequence.
         // or cleared on every round if decared on a higher level scope.
         cutType = bestCut(alignerForwardh, alignerForwardt, 
                           alignerBackwardh, alignerBackwardt, fasq, rightq, vrbm, LOG);
         if (cutType == 0) ++count_cut;
         else if (cutType == 1) ++count_discard;
         else {
            ++count_nocut;
         }

         if (fasq.length() >= vrbm["shortest"].as<unsigned int>()) {
            fasq.write(ouf);
         }
         if (!rightq.empty() && rightq.length() >= vrbm["shortest"].as<unsigned int>()) {
            rightq.write(ouf);
         }
      }
      cerr << i << " sequences processed\n";
   }
   catch (exception &e) {
      cerr << e.what() << endl << " failed to cut sequence\n";
      return 1;
   }
   inf.close();
   ouf.close();
   cout << "Final min max score:\n";
   cout << fasq.getMinScore() << "  " << fasq.getMaxScore() << endl;
   cout << "result written to " << outfile << endl
      << "   " << count_nocut << " no cut\n"
      << "   " << count_cut << " cut\n"
      << "   " << count_discard << " discarded\n";
   float ratioDC = count_discard/(float)count_cut;
   if (ratioDC > 0.2) {
      cout << "WARNING: Cut site may be wrong!\n";
   }
   cout << " discard/cut ratio = " << ratioDC << endl;
   LOG.close();
   cerr << "log file: " << logfile << " written for debug\n";

   return 0;
}

bool isBadAlign(const Dynaln<SimpleScoreMethod> &aligner, const po::variables_map &vm) {
   pair<int,int> gaps = aligner.numgaps();
   return aligner.getIdentity() < vm["iden-discard"].as<float>() 
         || gaps.first > vm["gap-discard"].as<int>() 
         || gaps.second > vm["gap-discard"].as<int>();
}

bool matchSite(const Dynaln<SimpleScoreMethod> &aligner, const po::variables_map &vm) {
   /* too many non-specific matches
   return (aligner.getIdentity() > vm["identity-cutoff"].as<float>()  
         && aligner.getAlnlen() > vm["alnlen-cutoff"].as<int>())
      || (aligner.getIdentity() > 0.999999 && aligner.getAlnlen() > 8);
      */
   return aligner.getIdentity() > vm["identity-cutoff"].as<float>()  
         && aligner.getAlnlen() > vm["alnlen-cutoff"].as<int>();
}



/*
 * Objects cannot be const, because other function will alter 
 * the object.
 * This is just a helper function.
 */
Dynaln<SimpleScoreMethod> *betterAligner(Dynaln<SimpleScoreMethod> &a1, 
      Dynaln<SimpleScoreMethod> &a2) 
{
   Dynaln<SimpleScoreMethod> *p;
   if (a1.getScore() > a2.getScore()) {
      p = &a1;
   }
   else {
      p = &a2;
   }
   return p;
}

int bestCut(Dynaln<SimpleScoreMethod> &alnForward1, Dynaln<SimpleScoreMethod> &alnForward2, 
      Dynaln<SimpleScoreMethod> &alnBackward1, Dynaln<SimpleScoreMethod> &alnBackward2,
      Fastq &fsq, Fastq &rightq, const po::variables_map &parm, ostream &ous) 
{
   DNA rawseq(fsq.getName(), fsq.getSequence());
   alnForward1.setSeq2(rawseq); 
   alnForward2.setSeq2(rawseq); 
   alnBackward1.setSeq2(rawseq);
   alnBackward2.setSeq2(rawseq);

   // run the local algorithm without generating the actual alignment info
   alnForward1.local(); 
   alnForward2.local(); 
   alnBackward1.local();
   alnBackward2.local();

   Dynaln<SimpleScoreMethod> *betterf, *betterb;

   betterf = betterAligner(alnForward1, alnForward2);
   betterb = betterAligner(alnBackward1, alnBackward2);

   if (betterf->getScore() > betterb->getScore()) {
      betterf->buildResult();
      return cutAtHead(betterf, fsq, rightq, parm, ous);
   }
   else {
      betterb->buildResult();
      return cutAtTail(betterb, fsq, rightq, parm, ous);
   }
}

int cutAtHead(Dynaln<SimpleScoreMethod> *aligner, Fastq &fsq, Fastq &rightq, const po::variables_map &parm, ostream &ous) {
   int actionType=0;
   try {
      if (matchSite(*aligner, parm)) {
         // top is site  =====>
         // bottom seq -----------
         //aligner.printAlign(ous);
#ifdef DEBUG
         ous << endl;
         aligner->printAlign(ous);
#endif
         if (isBadAlign(*aligner, parm)) {
#ifdef DEBUG
            ous << "forward poor align discarding:\n";
            ous << endl << fsq << endl;
#endif
            // discard left margin to right end
            rightq.clear();
            fsq.clear();
            actionType = 1;
         }
         else { // cut at left margin
            int cutidx=aligner->bottomBeginIndex();
            if (aligner->topBeginIndex() > 0) {
               cutidx -= aligner->topBeginIndex();
            }
            if (cutidx <= 0) {
               // ======
               // ------------ no need to cut
#ifdef DEBUG
               ous << "forward site match sequence beginning, no need to cut:\n"
                  << fsq << endl;
#endif
               rightq.clear();
            }
            else {
               fsq.cutAt(cutidx, rightq); 
#ifdef DEBUG
               ous << "forward broken into two:\n"
                  << fsq << endl << rightq << endl;
#endif
            }
         }
      }
      else {
         actionType = -1;
      }
   }
   catch (exception &e) {
      cerr << e.what() << " failed to process the forard cut-site\n";
   }

   return actionType;
}

int cutAtTail(Dynaln<SimpleScoreMethod> *aligner, Fastq &fsq, Fastq &rightq, const po::variables_map &parm, ostream &ous) {
   int actionType = 0;
   try {
      int cutIdx;
      if (matchSite(*aligner, parm)) {
         //   <=====
         // -----------
#ifdef DEBUG
         ous << endl;
         aligner->printAlign(ous); // only log matched 
#endif
         if (isBadAlign(*aligner, parm)) { // disard all
#ifdef DEBUG
            ous << "backward align is bad:\n";
            ous << "\ndiscarding fastq " << fsq << endl;
#endif
            rightq.clear();
            fsq.clear();
            actionType=1;
         }
         else { // break apart at tail (right-hand side) end
            //   <======
            //----------x-------
            cutIdx = aligner->bottomEndIndex() + 1;
            if (aligner->topEndIndex() < aligner->getSeq1Length()-1) {
               //   <===//
               // ------x----
               //   ->    x forward x assuming unmatched part matched
               cutIdx += (aligner->getSeq1Length() -1 - aligner->topEndIndex());
            }
            if (cutIdx > aligner->getSeq2Length() - 1) {
               rightq.clear();
               // no need to cut
            }
            else {
               fsq.cutAt(cutIdx, rightq);
            }
#ifdef DEBUG
            ous << "backward of original sequence broken into two:\n"
               << fsq << endl << rightq << endl;
#endif
         }
      }
      else {
         actionType = -1;
      }
   }
   catch (exception &e) {
      cerr << "Failed to process the backward cut\n" << e.what() << endl;
      throw;
   }
   return actionType;
}

string getFileStem(const string &fileName) {
   string::size_type i = fileName.rfind('.');
   string stem = fileName;
   if (i != string::npos) {
      stem = fileName.substr(0, i);
   }
   return stem;
}

