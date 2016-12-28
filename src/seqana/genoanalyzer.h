#ifndef GENOANALYZER_H
#define GENOANALYZER_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <cmath>
#include <iostream>
#include <fstream>

#include <dynalnt.h>
#include <alninfo.h>
#include <bioseq.h>
#include <fastq.h>
#include <stddev.h>

#define DEBUG

using namespace std;
using namespace orpara;


float roundPercent(float frac);
template<class T>
class SortPairByInt {
   public:
      bool operator()(const pair<T, int> &p1, const pair<T, int> &p2) {
         return p1.second > p2.second;
      }
};

class BestHaplotype {
   private:
      DNA haplotype;
      string gbacc;
      int score;
      double identity; 

   public:
      BestHaplotype(const DNA &haplo, const string &name, int sc, double iden) 
         : haplotype(haplo), gbacc(name), score(sc), identity(iden) { }
      BestHaplotype(const BestHaplotype &bh) 
         : haplotype(bh.haplotype), 
           gbacc(bh.gbacc), score(bh.score), 
           identity(bh.identity)
      { }
      BestHaplotype(BestHaplotype &&bh) 
         : haplotype(std::move(bh.haplotype)), 
           gbacc(std::move(bh.gbacc)), 
           score(bh.score), 
           identity(bh.identity)
      { }
      friend ostream& operator<<(ostream &ous, const BestHaplotype &bh) {
         ous << "Best Haplotype\nGBACC\tscore\tidentity\n";
         ous << bh.gbacc << "\t" << bh.score << "\t" << bh.identity << endl
            << bh.haplotype;
         return ous;
      }
};

// helper function

bool gapInPoly(const string &top, const string &bottom, size_t i);
void crossBottom(const string &top, string &bottom, size_t idx);

/**
 * This is an accumulator type of analyzer
 * The first version was designed for amplicaon sequence where only one
 * direction of the reads are sequenced.
 * In my second version I am trying to add shotgun approach sequence sequence
 * to the algorithm.  This will make the code hard to maintain.  
 * Only when I am having time, I may need to redesign the algorithm so that
 * I can easily extend and maintain the code. This will require at least 
 * two weeks of time so far I have not the time yet.
 */
class GenotypeAnalyzer {

   protected:
      /**
       * Reference used for analysis.
       * At this point, we did not use the actual sequence.
       * This should be a random one from the library if not given
       * by the caller. Should have another version of the constructor.
       */
      DNA refseq;
      /**
       * Using fastq so that I can compute scores.
       */
      vector<DNAQual> reads;
      /**
       * This is the position that the user want to focus on.
       * If the reads are not long enough, this is meaningless.
       * This should only be defined for amplicon type experiments.
       * For shotgun type experiment, this should be set to empty 
       * because you will have all possible combination of two or more
       * residues. This not make the user unpleasant.
       */
      vector<int> genotypeIndex; //{51, 57, 60, 63, 129, 141, 153, 243, 246};

      static SimpleScoreMethod scoreMethod; //(10, -9, -20, -3);
      /**
       * Identity cutoff for identity of the alignment.
       * Default is 0.85
       */
      static double identityCut;
      /**
       * cutoff value for alignment on the reference sequence
       * below this value, algiments will not be counted.
       * Default is 200 for Iontorrent platforms where 400 nt reads
       * are generated. For Illumina
       * this value should be set by the user to about 1/3 of the
       * read length.
       */
      static int alnlengthCut;
      /**
       * The unique refseq amplicon file containing all ref sequences.
       */
      static string reflibFile;
      /**
       * The library of reference sequence indexed by the GenBank 
       * Accession 
       * This is the result from reading the reflibFile.
       */
      static map<string, DNA> reflib;

      /**
       * this is the machine used in this class
       */
      Dynaln<SimpleScoreMethod> aligner;
      /** 
       * whether to fill the gap when calcularting the genotype
       */
      bool fillgap;

      /**
       * For quality control.
       */
      map<AlignInfo, int> alnsummary;

      /**
       * For debug.
       * Output for the actual alignment.
       */
      ofstream oaln;
      //string alnoutFile;

      /**
       * Single base count
       */
      vector<map<char, int> > bases;
      /**
       * Insertion of bases one after the base.
       */
      vector<map<string, int> > inserts;
      /**
       * Codon results of the accumulator
       * Each position is a table of codon => count
       */
      vector<map<string, int> > codons;
      /**
       * Fastq quality score average, std
       */
      vector<stddev> quality;

      /**
       * This is determined by the size of the targetCodonIndex.
       * The string element of the combined codon of the positions
       * specified by targetCodonIndex
       * The genotype is sorted by the codon combined string.
       * This only make sense for amplicon type of experiments.
       */
      map<string, int> genotype;
      /**
       * Derived result from genotype
       */
      map<string, int> genotypeAA;

      /**
       * Codons are sorted by their count, highest count appears first.
       * Another format of the result. Stored for efficiency.
       */
      vector<vector<pair<string, int> > > codonsByCount;

      /**
       * Genotype sorted by the highest to lowest count
       * This is by codon. The by amino acids should be smaler.
       */
      vector<pair<string, int> > genotypeByCount;
      /**
       * For amino acid level from highest to lowest frequency.
       */
      vector<pair<string, int> > genotypeAAByCount;

      /**
       * final consensus sequence of the most frequent one
       * Calculated result. Stored for efficiency.
       */
      string consensus;
      
      /**
       * Generate genotype sorted by count descending.
       * the genotype field will be updated.
       */
      void sortGenotype();

      /**
       * After calling this function, the codonsByCount and consensus
       * field will be updated.
       */
      void sortCodons();
      /**
       * Helper function to build the genotype on Amino Acids.
       * This will collpase the codon to have fewer categories
       * The genotypeAAByCount field will be modified as a result
       * of calling this function.
       */
      void buildAAGenotype();

      // helper functions
      void updateAlninfo();
      bool saveFailedRead(const DNAQual &read, ostream &oufail) const;
      /**
       * 0 for forward, 1 for backward, 2 for both forward and backward
       */
      int readDirection;
      /**
       * the reading frame of the reference sequence.
       * This if for codon recognition
       * Default 0.
       */
      int frame;
      /**
       * process one read in the forward direction.
       * @return true if passed quality cutoff, else return false.
       */
      bool consumeOneSequence(DNAQual &read); 
      /**
       * Consume one sequence in the reverse complement direction.
       */
      bool consumeOneSequenceRC(DNAQual &read);
      /**
       * try to align sequence to reference in both directions of the read.
       * @return + for aligned to the plus, - for minus, ? for none.
       */
      char consumeOneBothDirections(DNAQual &read);
      string inputFile1;

      /**
       * Helper only used by this class
       */
      static void shrinkGap(string &top, string &bottom);
      /**
       * Treat gap differently based on the length
       * if gap is corrected then it is in lowercase so this
       * information will be used by the downstream procedure.
       * @param top is always the reference.
       * @param bottom is the read.
       */
      static void correctGap(string &top, string &bottom);

   public:

      /**
       * default constructor.
       * Will set the reference sequence from the default container.
       */
      GenotypeAnalyzer()  
         : refseq(reflib.begin()->second), aligner(scoreMethod), fillgap(false), 
            alnsummary(), 
            bases(refseq.length()), inserts(refseq.length()),
            codons(ceil((float)refseq.length()/3)), 
            quality(refseq.length()), 
            genotype(), genotypeAA(),
            genotypeIndex(), 
            codonsByCount(ceil((float)refseq.length()/3)),
            genotypeByCount(), genotypeAAByCount(), 
            consensus(), readDirection(0), inputFile1(), frame(0)  { 
            aligner.setSeq1(refseq); 
      }
      /**
       * @param ref reference sequence to use.
       * @param targetIdx target codon position of linked interest.
       */
      GenotypeAnalyzer(const DNA &ref, const vector<int> &targetIdx) 
         : refseq(ref), aligner(scoreMethod), fillgap(false), alnsummary(),  
           bases(ref.length()), inserts(ref.length()),
           codons(ceil((float)ref.length()/3)), 
           quality(refseq.length()), 
           genotype(), genotypeAA(),
           genotypeIndex(targetIdx),
           codonsByCount(ceil((float)ref.length()/3)),
           genotypeByCount(), genotypeAAByCount(), consensus(),
           readDirection(0), inputFile1(), frame(0)
      {
           aligner.setSeq1(refseq); 
      }
      /**
      * This version uses the first sequence from static interal reflib.
      * @param targetIdx target codon position of linked interest.
      */
      GenotypeAnalyzer(const vector<int> &targetIdx) 
         : refseq(reflib.begin()->second), aligner(scoreMethod), fillgap(false), 
           alnsummary(),  
           bases(refseq.length()), inserts(refseq.length()),
           codons(ceil((float)refseq.length()/3)), 
           quality(refseq.length()), 
           genotype(), genotypeAA(),
           genotypeIndex(targetIdx),
           codonsByCount(ceil((float)refseq.length()/3)),
           genotypeByCount(), genotypeAAByCount(), consensus(),
           readDirection(0), inputFile1(), frame(0)
      {
           aligner.setSeq1(refseq); 
      }
      
      ~GenotypeAnalyzer() { }

      /**
       * Accumulate input from alignment string.
       */
      //void accumulate(const string &topaln, const string &bottomaln, size_t tb=0);
      /**
       * Version for counting quality scores
       * @param topaln alignment of the top sequence.
       * @param bottomaln alignment string from the bottom sequence.
       * @tb top sequene begin index 0-based
       * @bb bottom sequence begin index 0-based.
       */
      void accumulate(const int* bqual);

      /**
       * Will align all sequences from reads to the reference sequence
       * @param extReads external read buffer.
       * This way multiple read sources can be used.
       */
      //void consume(const vector<DNA> &extReads);
      /**
       * Use internal read buffer.
       * @return the number of sequences passed cutoffs and number of sequence
       * > 250nt and failed the cutoff.
       */
      pair<int,int> consume();
      pair<int,int> consumeForward();

      /**
       * @return a copy of the refseq. If the getBestHaplotype has been called
       *   The returned refseq will the best match from the consensus to the 
       *   haplotype.
       */
      DNA getRefseq() const {
         return refseq;
      }

      /**
       * Set a new reference sequence
       * All reference sequences should be the same length as the ones
       * in the library.
       */
      void setRefseq(const DNA &ref);
      /**
       * Set the reading frame of the reference sequence
       */
      void setFrame(int f) { frame=f; }
      int getFrame() const { return frame; }

      /**
       * Reformat the result of the analysis for user consumption.
         void reorderResult() {
            sortCodons();
            sortGenotype();
         }
       */
      void reorderResult();
      BestHaplotype computeBestHaplotype();
      void buildConsensus();
      string getConsensus() { 
         if (consensus.empty()) buildConsensus();
         return consensus; 
      }
      void setGapfill(bool gf) { fillgap=gf; }
      void fillgapYes() { fillgap=true; }
      void fillgapNo() { fillgap=false; }
      /**
       * Clear all fields for reuse.
       */
      void clear();

      //void run();

      /**
       * @return the length of the reference sequence.
       */
      size_t getRefseqLength() const { return refseq.length(); }
      /** help for making up fail file name
       */
      string getRefseqShortName() const;
      /**
       * Output finctions.
       */
      void printBases(ostream &ous) const;
      void printCodons(ostream &ous) const;
      void printGenotype(ostream &ous) const;
      void printConsensus(ostream &ous) const;
      void printGenotypeAA(ostream &ous) const;
      void printQuality(ostream &ous) const;

      /**
       * print the result of this algorithm
       *  This include: genotype, codons frequency, and consensus sequence.
       *  If not calculating linked mutations frequencies, then
       *  genofile will be empty.
       *  @param codonfiles output file name for codon
       */
      void printResult(const string &codonfiles, const string &genofile, 
            const string &consensusfile, const string &genoaafile, 
            const string &basefile, const string &qualfile) const;

      /**
      * Global variable use by tabulateCodon for efficiency
      * Start 0-based index of the codon of interest in the reference seq.
      * This could be changed before each analysis
      */

      /**
       * Helper method to convert map to sorted vector by count
       * This method could be used by the public.
       */
      //static vector<pair<string, int> >  map2vectorByCount(const map<string, int>& src); 
      void setAlignOutputFile(const string &fname);
      void printAlninfo(const string &alninfoFile) const;

      void setGenotypeIndex(const vector<int>& li) {
         genotypeIndex=li;
      }
      vector<int> getGenotypeIndex() const {
         return genotypeIndex;
      }

      /**
       * Set read input source
       * Make a copy, expensive. try not do this.
       */
      //void setReads(const vector<DNAQual> &reads) { this->reads = reads; }
      /**
       * Read from a fastq file.
       * @param file name of a fastq file.
       */
      string suckupReads(const string &file);

      /**
       * tell the algorithm the direction of the reads
       */
      void setReadDirection(int direction) { readDirection=direction; }
      void setBothDirections() { readDirection=2; }
      int getReadDirection() const { return readDirection; }
      bool isReadBothDirections() const { return readDirection == 2; }

      string constructFailFile() const;
      void trimAlignmentTail(string &top, string &bottom);

      /**
       * select a particular reference sequence library file.
       */
      static void setReflibFile(const string &file) { 
         reflibFile = file; loadReflib(); 
      }
      /**
       * @return the reflib file name. This is a full path usually.
       */
      static string getReflibFile() { return reflibFile; }
      /**
       * For loading the reference sequence library.
       * This methods must be called before usage of the algorithm.
       */
      static void loadReflib();

      static void setAlnlengthCut(int alct) { alnlengthCut=alct; }
      /**
       * Length cutoff for alignment.
       */
      static int getAlnlengthCutoff() { return alnlengthCut; }
      static int getAlnlenCutoff() { return alnlengthCut; }
      static int getAlnlengthCut() { return alnlengthCut; }
      /**
       * 1.5x alnlengthcut
       */
      static int getExpandedAlnlengthCutoff() { return ceil(1.5*alnlengthCut); }

      static void setIdentityCut(double icut) { identityCut = icut; }
      static double getIdentityCut() { return identityCut; }
      static double getIdentityCutoff() { return identityCut; }
};

/**
 * The input is paired end read.
 */
class GenotypeAnalyzerPair : public GenotypeAnalyzer {
   private:

      /** for the second pair */
      map<string, DNAQual> reads2;

   public:
      GenotypeAnalyzerPair() : GenotypeAnalyzer() { }
      /**
       * Reads in the second pair 
       */
      void suckupReads2(const string &file);
      void suckupReadPair(const string &forwardFile, const string &backwardFile) {
         suckupReads(forwardFile);
         suckupReads2(backwardFile);
      }
      /*
      void setReads(const vector<Fastq> &reads1, const vector<Fastq> &reads2) {
         GenotypeAnalyzer::setReads(reads1); 
         this->reads2=&reads2;
      }
      */

      /**
       * The paired version of consume.
       * This method can be very picky such that if only one of the
       * two ends match to the reference, it will give you a warning.
       * @return match detail as a map. + for the same direction as the 
       *   reference. - for opposite direction as reference. ? for not
       *   mapped to reference.  
       */
      map<string, int> consumePair();
      /**
       * read in reads2 but not in reads
       * @return a vector of pointers to the second entry for reads2.
       *   basically, you could modify these objects stored in the map.
       */
      vector<DNAQual*> getSingles();

      /** 
       * For testing, when done should be moved to the parent class.
       * @return true if alignment to reference successful.
      bool consumeOneSequence(const Fastq &read, ofstream &oufail);
      bool consumeOneSequenceRC(const Fastq &read, ofstream &oufail);
       */
};

#endif
