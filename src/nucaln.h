#ifndef NUCALN_H
#define NUCALN_H

#include <queue>
#include <string>
#include <stdexcept>

using namespace std;

class CandidateRemoval {
   public:
      CandidateRemoval(bool top, int pos, int numBase) : isTop(top), idx(pos), count(numBase) { }
      bool less(CandidateRemoval r1, CandidateRemoval r2) const { return r1.count < r2.count; }
      bool operator<(const CandidateRemoval &r) const { return count < r.count; }
      int getIndex() const { return idx; }
      void lowerIndex() { --idx; }

   private:
      bool isTop;
      int idx;
      int count;
};

class Nucaln {
   public:
      Nucaln() : top(), bottom() { }
      Nucaln(const Nucaln &nca);
      Nucaln(string aln1, string aln2) : top(aln1), bottom(aln2),
            candidate() { 
         if (top.length() != bottom.length()) throw runtime_error("top bottom not same length"); }

      /**
       * @return difference between length and expected length (24).
       */
      int getDiff() const;
      size_t length() const { return top.length(); }
      int getNumCandidate() const { return candidate.size(); }
      bool fixLength();
      /**
       * return pair of top,bottom sequences if they differ.
       * If the top bottom sequences are identical, the only the first element
       * will be set.
       * If there is no way for this function to shorten the longer alignments,
       * then the original sequences with gap char will be returned.
       */
      pair<string, string> getConsensus();
      void updateCandidate(const size_t fx);
      friend ostream& operator<<(ostream &ous, const Nucaln &aln);
      Nucaln& operator=(const Nucaln &nca);
      bool empty() const { return top.empty(); }
      bool hasCandidate() const { return candidate.size()>0; }

      string getUngappedTop() const;
      string getUngappedBottom() const;

   private:
      void findCandidate();
      static const char gapchar;
      static const int len = 24;
      /**
       * top alighment sequence with gapchar */
      string top; // top and bottom should be the same length
      string bottom;
      priority_queue<CandidateRemoval> candidate;
};

#endif
