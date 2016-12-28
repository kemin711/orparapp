#ifndef MISMATCHMAPPING_H
#define MISMATCHMAPPING_H

#include <string>
#include <map>

using namespace std;

/*
template<class T1, class T2> 
void printMap(const map<T1, T2> &m, ostream &ous) {
   typename map<T1, T2>::const_iterator it = m.begin();
   while (it != m.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}
*/

// thinking replacing this one with the above one
//void printMissedOne(const map<pair<string, string>, int> &miss, ostream &ous);

/**
 * for sorting string from right to left
 */
class ReverseSort {
   public:
      ReverseSort() { }
      bool operator()(const string &s1, const string &s2) const {
         for (int i = s1.size()-1; i > -1; --i) {
            if (s1[i] < s2[i]) { 
               return true; 
            }
            else if (s1[i] > s2[i]) {
               return false;
            }
         }
         return false;
      }
};

/**
 * Two string are different. 
 * @return true if they differ by one base.
 */
bool diffByOne(const string &s1, const string &s2);
/**
 * @return true if 2 or more differences.
 */
bool diffByTwoOrMore(const string &s1, const string &s2);

///// reverse direction /////////////
bool diffByOneReverse(const string &s1, const string &s2);
bool diffByTwoOrMoreReverse(const string &s1, const string &s2);

/*
template<class T1, class T2, class cmp> 
void printMapComp(const map<T1, T2, cmp> &m, ostream &ous) {
   typename map<T1, T2, cmp>::const_iterator it = m.begin();
   while (it != m.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}
*/


/**
 * Algorithm to find mismatch by sorting.
 * The algorithm simply retabulate the result by considering input sequences
 * that are different from the library sequence by one as a true hit.
 * The the numbers are a little better.
 */
class MismatchMapping {
   public:
      /**
       * constructor 
       * @param rfc input containing the count of for library sequences with
       *        negative number representing non-library sequences. This is the
       *        reuslt from the stemloop algorithm.
       *        This object will make a reference to this data. The algorithm
       *        may also modify the input data!
       */
      MismatchMapping(map<string, int> &rfc) : refcnt(rfc), refcntrev(0), missOne(), 
            summary(), summaryrev(), numhitDist() { }
      ~MismatchMapping() { delete refcntrev; }

      /**
       * summarize the mapping result stored in refcnt.
       */
      void summarize() const ;
      /**
       * print the summary result for refcnt into ous.
       */
      void printSummary(ostream &ous) const;
      /**
       * print the refcnt result to file.
       */
      void printForwardToFile(const string &ouf) const;
      /**
       * summarize the mapping result stored in refcntrev.
       */
      void summarizeReverse() const;

      /**
       * The following methods are public usage.
       */
      void run();
      /**
       * Output the final result after the reverse algorithm.
       * The field label will be suffixed with _miss1.
       */
      void printFinalSummary(ostream &ous) const { printReverseSummary(ous); }

      /**
       * if summaryrev is empty, this function will do the summary for you.
       */
      void printFinalSummaryToFile(const string &ouf) { printReverseSummaryToFile(ouf); }
      void printMappingToFile(const string &ouf) const { printBackwardToFile(ouf); }
      void printDistribution(ostream &ous) const;
      //////////// above are for public usage /////////////////////////////////

      /// the following are for testing and debugging
      void printReverseSummary(ostream &ous) const;
      void printBackwardToFile(const string &ouf) const;
      void printSummaryToFile(const string &ouf) const;
      void printReverseSummaryToFile(const string &ouf) const;
      //void MismatchMapping::addSuffixToFinalField(const string &sfx) const;

   private:
      /// output functions
      void printForward(ostream &ous) const;
      void printBackward(ostream &ous) const;
      /**
       * will run the forward algorithm, then update the summary.
       */
      void forward();
      /**
       * Run the backward algorithm, then update the reverse summary.
       */
      void backward();

      //MismatchMapping(const MismatchMapping &o) {}
      //MismatchMapping() { }
      /**
       * Run one round of forward algorithm
       */
      int forwardRun();
      /**
       * Runn one round of backward algorithm
       */
      int backwardRun();

      /**
       * A reference to reference count map. So not making a copy of the
       * original data. Potentially change the original data too.
       */
      map<string, int> &refcnt;
      map<string, int, ReverseSort> *refcntrev;
      map<pair<string, string>, int> missOne;
      // following two serves as buffer
      mutable map<string, int> summary; // summary of refcnt
      mutable map<string, int> summaryrev; // summary of refcnt
      mutable map<int, int> numhitDist;  // number of hits distribution
};

#endif
