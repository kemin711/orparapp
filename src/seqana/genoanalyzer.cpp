#include "genoanalyzer.h"
#include "stagergap.h"
#include <algorithm>
#include <stack>
#include <codon.h>
#include <bioseq.h>
#include <iomanip>
#include <strformat.h>
#include <cassert>

SimpleScoreMethod GenotypeAnalyzer::scoreMethod(10, -9, -29, -3);
double GenotypeAnalyzer::identityCut=0.85;
int GenotypeAnalyzer::alnlengthCut=200;
string GenotypeAnalyzer::reflibFile="/ngs/ngs16/zhouke/virology/hcvseq/amplicon_ulib";
//string GenotypeAnalyzer::reflibFile="/ng14/zhouke/virology/hcvseq/amplicon_ulib";
map<string,DNA> GenotypeAnalyzer::reflib = map<string,DNA>();

// helper
float roundPercent(float frac) {
      return roundf(frac*100)/100;
}

// helper
//vector<pair<string, int> > GenotypeAnalyzer::map2vectorByCount(const map<string, int>& src) { 
//vector<pair<string, int> > map2vectorByCount(const map<string, int>& src) { 
template<class T>
vector<pair<T, int> > map2vectorByCount(const map<T, int>& src) { 
   /*
    * There is still a bug in the compiler
    * the following code works but cause a double free
    * segmentation fault. I have to use the traditional
    * old way that works fine.
   vector<pair<string,int> > elems(src.size());
   transform(src.begin(), src.end(), elems.begin(), [](const pair<string,int> &p) { return p; });
   sort(elems.begin(), elems.end(), [](const pair<string, int>& l, const pair<string,int>& r) {
            return l.second > r.second; });
   return elems;
   */
   vector<pair<T,int> > arr;
   for (auto itr=src.begin(); itr != src.end(); ++itr) {
      arr.push_back(pair<T, int>(itr->first, itr->second));
   }
   //SortPairByInt sorter;
   //sort(arr.begin(), arr.end(), sorter);
   sort(arr.begin(), arr.end(), SortPairByInt<T>());
   return arr;
}

/**
 * the top is reference and the bottom is read.
 * I will also fix some gaps that are caused by Ion Torrent
 * platform.
 * After adjustment of alignment parameters, this case should
 * no happed, so this method is no longer used.
 */
void GenotypeAnalyzer::shrinkGap(string &top, string &bottom) {
//void shrinkGap(string &top, string &bottom) {
   static const unsigned int w = 6;
   string::size_type i, j;
   stack<pair<unsigned int, unsigned int> > eraseIndex;

   for (i=0; i<top.size(); ++i) {
      if (top[i] == '-') {
         j=i-1;
         while (j > 0 && j > i-w && bottom[j] != '-') --j;
         if (bottom[j] == '-') { // case for shrinking
            //GCACGCCCT--CCCCGGC
            //GCAC-CCCTTCCCCCAGC
           eraseIndex.push(make_pair(i,j));
         }
         else { // try the right size
            j=i+1;
            while (j<top.size() && j < i+w && bottom[j] != '-') ++j;
            if (bottom[j] == '-') {
               eraseIndex.push(make_pair(i,j));
            }
         }
      }
      else if (bottom[i] == '-') {
         unsigned int b, e;
         if (i-w > 0) b = i-w;
         else b = 0;
         e = i + w;
         if (e > top.length()-1) e = top.length()-1;
         bool topHasGap = false;
         for (j=b; j<e; ++j) {
            if (top[j] == '-') {
               topHasGap = true;
               break;
            }
         }
         // correction algorithm
         // ACCACCTGCCCATGTGGAGCAC
         // ACCACCTGCC-ATGTGGAGCGC
         if (!topHasGap && (
            (
               ((i-2 > 0 && bottom[i-1] == bottom[i-2])
                     || (i+2 < bottom.length() && bottom[i+1] == bottom[i+2]))
               && ((i>0 & top[i] == bottom[i-1]) 
                  || (i<bottom.length()-1 && top[i] == bottom[i+1]))
            )
            || (i-2>0 && bottom[i-2] == 'C' && bottom[i-1] == 'A'
               && top[i-2] == 'C' && top[i-1] == 'A' && top[i] == 'A')
            )
         ) 
         {
            bottom[i] = top[i];
         }
      }
   }
   // do the erasing
   while (!eraseIndex.empty()) {
      pair<unsigned int, unsigned int> ij = eraseIndex.top();
      top.erase(ij.first, 1);
      bottom.erase(ij.second, 1);
      eraseIndex.pop();
   }
}

bool topHasGap(const string &top, size_t i, size_t window) {
   size_t j = i-1;
   while (j>0 && j > i-window) {
      if (top[j] == '-') return true;
      --j;
   }
   j = i+1;
   while (j<top.size() && j< i+window) {
      if (top[j] == '-') return true;
      ++j;
   }
   return false;
}
bool topNoGap(const string &top, size_t i, size_t window) {
   return !topHasGap(top, i, window);
}

/*
 * Example cases
 * GATTTCAAA    GATTTCAA  TCATAAAAAACGACT TCATAAAAAAAACGACT
 * GATT-CAAA    GA-TTCAA  TCATAAAA--CGACT TCAT--AAAAAACGACT
 * ----
 *  AAACCA
 *  AAAC-A
 *  bottom has gap, top no gap. 
 *  Call this method with top bottom swaped if checking the
 *  other direction.
 */
bool gapInPoly(const string &top, const string &bottom, size_t i) {
   // gap not at start or end
   if (int(i)-1 < 0 || i+1 >= top.size()) return false; 
   if (bottom[i+1] != '-') { // single gap
      //cerr << top[i-1] << ' ' << top[i] << ' ' << bottom[i-1] << endl;
      if (top[i] == top[i-1] && top[i] == bottom[i-1]) {
         //cerr << "gap at right side end of repeat\n";
         return true;
      }
      //cerr << top[i+1] << ' ' << top[i] << ' ' << bottom[i+1] << endl;
      if (top[i] == top[i+1] && top[i] == bottom[i+1]) {
         //cerr << "gap at left side end of repeat\n";
         return true;
      }
   }
   // multiple gap
   size_t x = i+1;
   while (x<top.size() && bottom[x]=='-') {
      if (top[x] != top[i]) return false;
      ++x;
   }
   // now all gap in the poly
   if (x-1 > i) { // more than one gap
      if (top[i] == top[i-1] && top[i] == bottom[i-1]) {
         //cerr << "multiple gaps at right side end of repeat\n";
         return true;
      }
      else {
         i=x-1;
         if (i+1 >= top.size()) return false;
         if (top[i] == top[i+1] && top[i] == bottom[i+1]) {
            //cerr << "multiple gaps at left end of repeat\n";
            return true;
         }
      }
   }
   return false;
}

void displayGap(const string &top, const string &bottom, size_t i, size_t window) {
   size_t b, e;
   int bb;
   bb = i-window;
   if (bb < 0) bb = 0;
   b = bb;
   e = i+window;
   if (e > top.size()) e = top.size();
   cout << "gap region:\n";
   cout << top.substr(b, e-b) << endl;
   cout << bottom.substr(b, e-b) << endl << endl;
}

int countGapsAt(const string &str, size_t idx) {
   int total=0;
   while (idx < str.size() && str[idx] == '-') {
      ++total;
      ++idx;
   }
   return total;
}

int gapsInBound(const string &str, int idx, int w) {
   int b = idx-w;
   if (b < 0) b=0;
   int e = idx+w;
   if (e > str.size()) e=str.size();
   int total=0;
   for (size_t i=(size_t)b; i<(size_t)e; ++i) {
      if (str[i] == '-') {
         ++total;
      }
   }
   return total;
}

void fillBottomGap(const string &top, string &bottom, size_t idx) {
  while (idx < bottom.size() && bottom[idx] == '-') {
     bottom[idx]=tolower(top[idx]);
     ++idx;
  }
}

// the bottom is extra insertion replace Base with X
void crossBottom(const string &top, string &bottom, size_t idx) {
   while (idx < bottom.size() && top[idx] == '-') {
      bottom[idx] = 'X';
      ++idx;
   }
}

void GenotypeAnalyzer::trimAlignmentTail(string &top, string &bottom) {
   static const int window=12;
   if (top.length() < 3*window) return;
   int identical=0;
   unsigned int i=top.length() - window - 1;
   while (i < top.length()) {
      if (top[i] == bottom[i]) ++identical;
      ++i;
   }
   while ((float)identical/window < 0.9*identityCut
         || (float)identical/window < 0.9*aligner.getIdentity()) {
      i=top.length()-window;
      while (top[i] == bottom[i]) ++i;
      //cout << "before trim tail:\n" << top << endl << bottom << endl;
      top.erase(i);
      bottom.erase(i);
      //cout << "after trim tail:\n" << top << endl << bottom << endl;
      // more rounds
      if (top.length() < 3*window) return;
      identical=0;
      i=top.length() - window - 1;
      while (i < top.length()) {
         if (top[i] == bottom[i]) ++identical;
         ++i;
      }
   }
}

void GenotypeAnalyzer::correctGap(string &top, string &bottom) {
   static const unsigned int gapw = 6;
   string::size_type i=0, j, x;
   while (i<top.size()) {
      if (top[i] == '-') {
         int gapshere = countGapsAt(top, i);
         if (gapshere % 3 == 0) { // multiple of 3
            //cerr << "reads has insertion of 3n\n";
            i += gapshere;
            continue;
         }
         else if (gapInPoly(bottom, top, i)) {
            // cross bottom
            crossBottom(top, bottom, i);
            i += gapshere;
            continue;
         }
      }
      if (bottom[i] == '-') {
         int gaplen=countGapsAt(bottom, i);
         if (gaplen % 3 == 0) {
            //cerr << "read has 3n deletion!\n";
            i += gaplen;
            continue;
         }
         int topgaps=gapsInBound(top, i, gapw);
         if (topgaps == 0) {
            if (gapInPoly(top, bottom, i)) {
               fillBottomGap(top, bottom, i);
               i += gaplen;
               continue;
            }
            else if (gaplen > 2) { // too long not filling
               //cerr << "not filling gaps longer than 2\n";
               i += gaplen;
               continue;
            }
            else {
               fillBottomGap(top, bottom, i);
               i += gaplen;
            }
         }
         else { // top also got gaps, make sure top gap is not facing
            Stagergap stagergap(top, bottom, i);
            if (stagergap.findAndFixLeft()) {
               ++i;
            }
            else if (stagergap.findAndFixRight()) {
               ++i;
            }
            else {
               // bottom X-ed read on the left
               // top has gaps with window
               //       S        Top same as bottom left or right 2 bases
               // NNNNNNNNNNNNNN
               // NNNNNN-NNNNNNN
               //     SS
               //        SS
               // NNNNNCAANNNNNN
               // NNNNNCA-NNNNNN
               if ((
                     ((int(i)-2 > 0 && bottom[i-1] == bottom[i-2])
                           || (i+2 < bottom.length() && bottom[i+1] == bottom[i+2]))
                     && ((i>0 & top[i] == bottom[i-1]) 
                        || (i < bottom.length()-1 && top[i] == bottom[i+1]))
                  )
                  || (int(i)-2 > 0 && bottom[i-2] == 'C' && bottom[i-1] == 'A'
                     && top[i-2] == 'C' && top[i-1] == 'A' && top[i] == 'A')
                  ) {
                  bottom[i] = tolower(top[i]);
               }
               //else {
               //   displayGap(top, bottom, i, gapw);
               //}
            }
         }
      }
      ++i;
   }
}
/*
void GenotypeAnalyzer::accumulate(const string &topaln, const string &bottomaln, size_t tb) 
{
   string top(topaln);
   string bottom(bottomaln);
   if (fillgap) {
      //GenotypeAnalyzer::shrinkGap(top, bottom);
      cerr << "before correct gap\n"
         << topaln << endl << bottomaln << endl;
      GenotypeAnalyzer::correctGap(top, bottom);
      //assert(top.length() == bottom.length());
      cerr << "after correct gap\n"
         << top << endl << bottom << endl;
   }

   // i for index in the alignment
   // refi is for index in the refseq
   string::size_type i, refi, readi;
   i=0;
   refi=0; // reference base index 
   if (tb > 0) {
      refi = tb;
      while (refi % 3 != 0) {
         if (top[i] == '-') {
            ++inserts[refi][string(1, bottom[i])];
         }
         else {
            ++refi;
            ++bases[refi][bottom[i]];
         }
         ++i;
      }
   }

   string codon(3, 'N');
   // selected target codon combined
   string combined; // all target codons
   unsigned int t = 0;
   bool insideInsert;
   string aInsert;
   while (i < top.length()) {
      if (top[i] != '-') {
         if (insideInsert) { // save the insert if came from an insert state
            ++inserts[refi][aInsert];
            insideInsert = false;
         }
         ++bases[refi][bottom[i]];
         codon[refi%3] = bottom[i];
         if (refi%3 == 2) {
            // if target codon reached
            if (!genotypeIndex.empty() && refi-2 == genotypeIndex[t]) {
               combined += codon;
               ++t;
            }
            ++(codons[refi/3][codon]);
         }
         ++refi;
      }
      else { 
         if (!insideInsert) {
            aInsert = string(1, bottom[i]); // clear the old one by new value
            insideInsert = true;
         }
         else {
            aInsert += bottom[i];
         }
      }
      ++i;
   }
   if (refi == getRefseqLength() && refi%3 > 0) { // dealing with full length
      //cerr << refi << " " << refi/3+1 << " " << codons.size() << endl;
      // input, we are going to take the partial codon.
      ++(codons[refi/3][codon.substr(0, refi%3)]);
   }
   // must be 100% coverage of all interested codons
   if (!genotypeIndex.empty() && combined.size() == genotypeIndex.size()*3) {
      //cerr << combined << endl;
      ++genotype[combined];
   }
   *
    * If the reads are not long enough you will go here,
    * not a problem.
#ifdef DEBUG
   else {
      cerr << "Genostring: " << combined
         << " should be " << genotypeIndex.size()*3 << " long "
         << " but you have " << combined.length() << "! I will discard it.\n";
   }
#endif
*
}
*/

// version for quality score
void GenotypeAnalyzer::accumulate(const int* bqual) 
{
   string top(aligner.getTopAln());
   string bottom(aligner.getBottomAln());
   size_t tb = (unsigned int)aligner.topBeginIndex();
   size_t bb = (unsigned int)aligner.bottomBeginIndex();
   //const int* bqual=static_cast<DNAQual>(aligner.getBottomSequence()).getQuality();
   // cannot cast bioinfo to DNAQual type.  

   if (fillgap) {
      correctGap(top, bottom);
      trimAlignmentTail(top,bottom);
      // for efficiency only doing these fine-tuned operations in the second
      // step
      // fix edge effect
      if (tb < 4 && bb < 4 && tb == bb) {
         top.insert(0, aligner.getTopSequence().substring(0,tb));
         bottom.insert(0, aligner.getBottomSequence().substring(0,bb));
         tb=0;
         bb=0;
      }
      int topUnalignedTail = aligner.getSeq1Length() - aligner.topEndIndex() - 1;
      int bottomUnalignedTail = aligner.getSeq2Length() - aligner.bottomEndIndex() - 1;
      if (topUnalignedTail > 0 && topUnalignedTail < 4 && bottomUnalignedTail >= topUnalignedTail) {
         top += aligner.getTopSequence().substring(aligner.getSeq1Length()-topUnalignedTail);
         bottom += aligner.getBottomSequence().substring(aligner.bottomEndIndex()+1, topUnalignedTail);
      }
   }

   // i for index in the alignment
   // refi is for index in the refseq
   string::size_type i, refi, readi;
   i=0;
   refi = tb;
   readi= bb;
   //cerr << "readi at start: " << readi << endl;
   // advance 1 or 2 bases, frame can be 0, 1, or 2
   while (refi % 3 != frame) {
      if (top[i] == '-') {
         if (bottom[i] == 'X') { // do nothing
         }
         else {
            ++inserts[refi][string(1, bottom[i])];
         }
         ++readi;
      }
      else { // bottom cannot be gap
         ++bases[refi][bottom[i]];
         if (bottom[i] == '-') {
            //quality[refi](0); // zero quality
            quality[refi](floor(0.8*bqual[readi])); // del quality
         }
         else {
            quality[refi](bqual[readi]);
            // use quality of preceeding base for gaps filled
            if (isupper(bottom[i])) ++readi;
         }
         ++refi;
      }
      ++i;
   }

   string codon(3, 'N');
   // selected target codon combined
   string combined; // all target codons
   unsigned int t = 0;
   bool insideInsert=false;
   string aInsert;
   while (i < top.length()) {
      if (top[i] != '-') { // out of insert state
         if (insideInsert) { // save the insert if came from an insert state
            ++inserts[refi][aInsert];
            insideInsert = false;
         }
         ++bases[refi][bottom[i]];
         if (bottom[i] == '-') {
            //quality[refi](0);
            quality[refi](floor(0.8*bqual[readi])); // del quality
         }
         else {
            quality[refi](bqual[readi]);
            if (isupper(bottom[i])) ++readi;
         }
         codon[(refi-frame)%3] = bottom[i];
         if ((refi-frame)%3 == 2) { // read end of a full codon
            if (!genotypeIndex.empty()) {
               if (t<genotypeIndex.size()) {
                  // if target codon reached, and asked for linked genotype
                  if ((refi-2) == genotypeIndex[t]) {
                     combined += codon;
                     ++t; // advance to next target codon
                  }
                  //else {
                  //   cerr << "t " << t << " wanted index: " << genotypeIndex.size() 
                  //      << " " << refi-2 << " | " << genotypeIndex[t] << endl;
                  //   cerr << "not collecting\n";
                  //}
               }
            }
            ++(codons[(refi-frame)/3][codon]);
         }
         ++refi;
      }
      else { // top is gap char, bottom cannot be gap
         if (bottom[i] == 'X') { // ignore
         }
         else {
            if (!insideInsert) {
               aInsert = string(1, bottom[i]); // clear the old one by new value
               insideInsert = true;
            }
            else {
               aInsert += bottom[i];
            }
         }
         ++readi;
      }
      ++i;
   }
   // for debug
   //cerr << "readi at end: " << readi << " seq2len: " << aligner.getSeq2Length() << endl;
   //assert(refi <= refseq.length());
   //assert(readi <= aligner.getSeq2Length());
   //assert(i <= bottom.length());
   //cout << "at end i and length of align: " << i << " " << aligner.getSeq2Length() << endl;

   if (refi == getRefseqLength() && (refi-frame)%3 > 0) { // dealing with full length
      //cerr << refi << " " << refi/3+1 << " " << codons.size() << endl;
      // input, we are going to take the partial codon.
      ++(codons[(refi-frame)/3][codon.substr(0, (refi-frame)%3)]);
   }
   // must be 100% coverage of all interested codons
   if (!genotypeIndex.empty() && combined.size() == genotypeIndex.size()*3) {
      ++genotype[combined];
   }
   //else {
   //   cerr << "combined.size() " << combined.size() << endl;
   //   cerr << "failed to collect genotype index asked " << genotypeIndex.size() << endl;
   //}
}

/*
void GenotypeAnalyzer::consume(const vector<DNA> &extReads) {
   clear();
   int every = 10;
   size_t i;
   for (i=0; i<extReads.size(); ++i) {
      if (extReads[i].length() < alnlengthCut) continue;
      aligner.setSeq2(extReads[i]);
      aligner.runlocal();
      if (aligner.getIdentity() > identityCut && aligner.getSeq1AlignedLength() > alnlengthCut) {
         accumulate(aligner.getTopAln(), aligner.getBottomAln(),
               (unsigned int)(aligner.topBeginIndex()));
      }
      ++alnsummary[AlignInfo(aligner.getScore(), roundPercent(aligner.getIdentity()),
            roundPercent(aligner.getCov1()), extReads[i].length())];
#ifdef DEBUG
      aligner.printAlign(oaln, 80);
#endif
      if (i%every == 0) {
         cerr << "working on reads " << i << " identity=" << roundPercent(aligner.getIdentity())
            << " refcov=" << roundPercent(aligner.getCov1()) << endl;
         every *= 2;
      }
   }
   cerr << i << " reads processed\n";
}
*/

bool GenotypeAnalyzer::consumeOneSequence(DNAQual &read) {
   aligner.setSeq2(read);
   aligner.runlocal();
   if (aligner.getIdentity() > identityCut 
         && aligner.getSeq1AlignedLength() > alnlengthCut) 
   {
      accumulate(read.getQuality());
            //aligner.getTopAln(), aligner.getBottomAln(),
         //read.getQuality(),
         //(unsigned int)(aligner.topBeginIndex()),
         //(unsigned int)(aligner.bottomBeginIndex())
         //);
      return true;
   }
   return false;
}

bool GenotypeAnalyzer::consumeOneSequenceRC(DNAQual &read) {
   //cerr << "using reverse complement\n";
   return consumeOneSequence(*read.getRevcomp());
}

char GenotypeAnalyzer::consumeOneBothDirections(DNAQual &read) {
   bool aligned=consumeOneSequence(read);
   char direc='?';
   if (aligned) direc='+';
   else {
      aligned=consumeOneSequenceRC(read);
      if (aligned) direc='-';
   }
   return direc;
}

bool GenotypeAnalyzer::saveFailedRead(const DNAQual &read, ostream &oufail) const {
   if (read.length()> getExpandedAlnlengthCutoff()) {
      if (fillgap) {
         oufail << read;
      }
      return true;
   }
   return false;
}

void GenotypeAnalyzer::updateAlninfo() {
   ++alnsummary[AlignInfo(aligner.getScore(), roundPercent(aligner.getIdentity()),
         roundPercent(aligner.getCov1()), aligner.getSeq2Length())];
#ifdef DEBUG
   if (fillgap) aligner.printAlign(oaln, 80);
#endif
}

string GenotypeAnalyzer::constructFailFile() const {
   string failFile = refseq.getName();
   if (refseq.getName().find('|') != string::npos) {
      vector<string> fields=split(failFile, '|');
      if (fields.size() > 3) {
         failFile=fields[3];
      }
      else {
         cerr << "failed to extract shorter name for " << failFile << endl;
      }
   }
   string inf=inputFile1.substr(0, inputFile1.rfind('.'));
   if (inf[inf.length()-2] == 'R' && isdigit(inf[inf.length()-1])) {
      inf=inf.substr(0, inf.length()-2);
   }
   if (inf.find('/') != string::npos) {
      inf=inf.substr(inf.find('/'));
   }
   //failFile += (itos(time(0)) + ".fail.fas");
   failFile += "_" + inf + ".fail.fas";
   return failFile;
}

// for the single direction experiment where all reads 
// are from the same direction
pair<int,int> GenotypeAnalyzer::consumeForward() {
   if (!genotypeIndex.empty())
      cerr << genotypeIndex.size() << " position to look after\n";
   clear();
   int every = 100;
   size_t i;
   int goodcnt = 0;
   int failcnt = 0;

   string failFile = constructFailFile();
   ofstream oufail(failFile.c_str());
   for (i=0; i < reads.size(); ++i) {
      bool aligned=consumeOneSequence(reads[i]);
      if (aligned) ++goodcnt;
      else {
         if (saveFailedRead(reads[i], oufail)) ++failcnt;
      }
      updateAlninfo();
      if ((i+1)%every == 0) {
         cerr << "working on reads " << i+1 << " identity=" << roundPercent(aligner.getIdentity())
            << " refcov=" << roundPercent(aligner.getCov1()) << endl;
         every *= 2;
      }
   }
   cerr << i << " reads longer than " << alnlengthCut << " processed\n"
        << goodcnt << " are good\n";
   cerr << failcnt << " longer than " << getExpandedAlnlengthCutoff()
      << " did not match to refseq well\n";
   return make_pair(goodcnt, failcnt);
}

pair<int,int> GenotypeAnalyzer::consume() {
   if (!genotypeIndex.empty())
      cerr << genotypeIndex.size() << " position to look after\n";
   clear();
   int every = 100;
   size_t i;
   int goodcnt = 0;
   int failcnt = 0;

   string failFile = constructFailFile();
   ofstream oufail(failFile.c_str());
   for (i=0; i < reads.size(); ++i) {
      char alndirection=consumeOneBothDirections(reads[i]);
      if (alndirection == '?') {
         if (saveFailedRead(reads[i], oufail)) ++failcnt;
      }
      else { // good alignment
         ++goodcnt;
      }
      updateAlninfo();
      if ((i+1)%every == 0) {
         cerr << "working on reads " << i+1 << " identity=" << roundPercent(aligner.getIdentity())
            << " refcov=" << roundPercent(aligner.getCov1()) << endl;
         every *= 2;
      }
   }
   cerr << i << " reads longer than " << alnlengthCut << " processed\n";
   cerr << failcnt << " longer than " << getExpandedAlnlengthCutoff()
      << " did not match to refseq well\n";
   return make_pair(goodcnt, failcnt);
}

// this function should be broken up into
// pure calculation, and output
BestHaplotype  GenotypeAnalyzer::computeBestHaplotype() {
   DNA cons("consensus", getConsensus());
   cout << "consensus " << endl << cons << endl;
   int score = 0;
   double identity=0;
   string best;
   aligner.setSeq1(cons);
   for (auto itr=reflib.begin(); itr != reflib.end(); ++itr) {
      //cout << "comparing consensus to " << itr->first << endl;
      aligner.setSeq2(itr->second);
      int tmp = aligner.runlocal();
      //cout << "score " << tmp << endl;
      if (tmp > score) {
         score = tmp;
         best = itr->first;
         identity = aligner.getIdentity();
      }
   }
   cout << "best ref: " << best << " score: " << score
      << " identity: " << identity << endl;
   auto bitr = reflib.find(best);
   //if (refseq.length() == bestHaplotype.length()) {
   if (refseq.length() == bitr->second.length()) {
      refseq = bitr->second;
      aligner.setSeq1(refseq); // now use the new refseq
   }
   else {
      cerr << "new haplotype has a different length from the original refseq!\n"
         << refseq << bitr->second << endl;
      exit(1);
   }
   //return refseq;
   return BestHaplotype(bitr->second, best, score, identity);
}

template<class T>
void clearVectorMap(vector<map<T, int> > &vm) {
   if (!vm.empty()) {
      for (size_t i=0; i<vm.size(); ++i) {
         vm[i].clear();
      }
   }
}

void GenotypeAnalyzer::clear() {
   clearVectorMap(bases);
   clearVectorMap(inserts);
   clearVectorMap(codons);
   if (!genotype.empty()) {
      genotype.clear();
   }
}

//void GenotypeAnalyzer::run() {
//   nogapfill();
//}

// produce a vector
void GenotypeAnalyzer::sortGenotype() {
   genotypeByCount = map2vectorByCount(genotype);
   /*
   std::transform(genotype.begin(), genotype.end(), genotypeByCount.begin(), 
         [](const pair<string,int> &p) { return p; }
   );
   sort(genotypeByCount.begin(), genotypeByCount.end(), 
         [](const pair<string, int>& l, const pair<string,int>& r) {
            return l.second > r.second;
         }
   );
   */
}

void GenotypeAnalyzer::sortCodons() {
   //if (!consensus.empty()) consensus.clear();
   for (size_t i=0; i < codons.size(); ++i) {
      vector<pair<string, int> > cod = map2vectorByCount(codons[i]);
      //consensus += cod.front().first;
      codonsByCount[i] = cod;
   }
}

void GenotypeAnalyzer::reorderResult() {
   sortCodons();
   if (!genotypeIndex.empty()) {
      cerr << "inside reorderResult\n";
      sortGenotype();
      buildAAGenotype();
   }
}

/*
 * count all nongap characters.
 */
int getTotalValue(const map<char, int> &bc) {
   int sum = 0;
   for (auto itr = bc.begin(); itr != bc.end(); ++itr) {
      sum += itr->second;
   }
   return sum;
}

void GenotypeAnalyzer::buildConsensus() {
   // user may give the wrong reference!
   consensus.clear();
   for (size_t i=0; i<bases.size(); ++i) {
      if (bases[i].empty()) {
         cerr << " no bases collected from read for index " << i << endl;
         cerr << " size of bases: " << bases.size() << endl;
         //exit(1);
         consensus.append(1, '-');
      }
      else {
         vector<pair<char, int> > bcount = map2vectorByCount(bases[i]);
         consensus.append(1, bcount[0].first);
      }
   }
}

void GenotypeAnalyzer::printBases(ostream &ous) const {
   ous << "position\tbase\tcount(frequency)\n";
   for (size_t i=0; i<bases.size(); ++i) {
      int totalBase = getTotalValue(bases[i]);
      ous << i+1;
      vector<pair<char, int> > bcount = map2vectorByCount(bases[i]);
      for (size_t b = 0; b<bcount.size(); ++b) {
         ous << '\t' << bcount[b].first << '\t' << bcount[b].second << '(' 
            << setprecision(8) << double(bcount[b].second)/totalBase << ')';
      }
      // output potentional insertion
      auto jtr = inserts[i].begin();
      if (jtr != inserts[i].end()) {
         bool markerOut = false;
         if (jtr->second > 100 && double(jtr->second)/totalBase > 0.1) {
            markerOut = true;
            ous << "\tinsert: " << jtr->first << '\t' << jtr->second 
               << '(' << setprecision(7) << double(jtr->second)/totalBase << ')';
         }
         ++jtr;
         while (jtr != inserts[i].end()) {
            if (jtr->second > 100 && double(jtr->second)/totalBase > 0.1) {
               if (markerOut) ous << '\t';
               else {
                  ous << "\tinsert: ";
                  markerOut = true;
               }
               ous << jtr->first << '\t' << jtr->second 
                  << '(' << double(jtr->second)/totalBase << ')';
            }
            ++jtr;
         }
      }
      ous << endl;
   }
}

void GenotypeAnalyzer::printCodons(ostream &ous) const {
   for (unsigned int i=0; i<codonsByCount.size(); ++i) {
      ous << i+1 << "\t";
      int totalCodons = 0;
      unsigned int j;
      for (j=0; j < codonsByCount[i].size(); ++j) {
         totalCodons += codonsByCount[i][j].second;
      }
      float frac;
      for (j=0; j < codonsByCount[i].size(); ++j) {
         ous << codonsByCount[i][j].first 
            << ":" << codonsByCount[i][j].second 
            << " (" << setprecision(8) 
            << (double)codonsByCount[i][j].second/totalCodons << ")\t";
      }
      ous << endl;
   }
}

// helper for printGenotype
// formats one line of genotype.
string formatGenotype(const string &geno) {
   // convert codon to translation
   codon cd;
   string ret;
   for (unsigned int i=0; i<geno.length(); i += 3) {
      ret += (geno.substr(i, 3) + ":" + cd[geno.substr(i,3)]);
      if (i < geno.length()-3) {
         ret += " | ";
      }
   }
   return ret;
}

void GenotypeAnalyzer::buildAAGenotype() {
   // make sure we don't combine multiple run results.
   if (!consensus.empty()) {
      genotypeAA.clear();
   }
   if (!genotypeAA.empty()) {
      genotypeAA.clear();
   }

   for (auto itr=genotype.begin(); itr != genotype.end(); ++itr) {
      string aa;
      translate(aa, itr->first, 1, 0); // translate all 
      if (genotypeAA.find(aa) != genotypeAA.end()) {
         genotypeAA[aa] += itr->second;
      }
      else {
         genotypeAA[aa] = itr->second;
      }
   }
   genotypeAAByCount = map2vectorByCount(genotypeAA);
}

void GenotypeAnalyzer::printGenotype(ostream &ous) const {
   if (genotypeIndex.empty()) {
      cerr << "empty genotype request, nothing will be done\n";
      return;
   }
   //ous << "Codons\tCount\tFrequency\n";
   size_t i;

   cerr << "outputing genotype ...\n";
   // output positions user specified for header
   for (i=0; i<genotypeIndex.size(); ++i) {
      ous << setw(8) << left << genotypeIndex[i];
      cerr << "base location " << genotypeIndex[i] << " ";
   }
   cerr << endl;

   // the fixed part of the header
   ous << "\tCount\tFrequency\n";

   int totalGenotype = 0;
   for (i=0; i<genotypeByCount.size(); ++i) {
      totalGenotype += genotypeByCount[i].second;
   }

   for (i=0; i<genotypeByCount.size(); ++i) {
      ous << formatGenotype(genotypeByCount[i].first) << "\t" 
         << genotypeByCount[i].second 
         << "\t" << setprecision(8) 
         << (double)genotypeByCount[i].second/totalGenotype << endl;
   }
}

void GenotypeAnalyzer::printGenotypeAA(ostream &ous) const {
   ous << "AA\tcount\n";
   for (size_t i=0; i<genotypeAAByCount.size(); ++i) {
      ous << genotypeAAByCount[i].first << "\t" 
         << genotypeAAByCount[i].second << endl;
   }
}

void GenotypeAnalyzer::printConsensus(ostream &ous) const {
   ous << consensus << endl;
}

void GenotypeAnalyzer::printQuality(ostream &ous) const {
   ous << "position\tqavg\tqstddev\tcount\n";
   for (auto i=0; i<quality.size(); ++i) {
      ous << i+1 << '\t';
      quality[i].print(ous) << endl;
   }
}

void GenotypeAnalyzer::printResult(const string &codonfile, 
      const string &genofile, const string &consensusfile, 
      const string &genoaafile, const string &basefile,
      const string &qualfile) const {
   ofstream oucod(codonfile.c_str());
   ofstream ouc(consensusfile.c_str());
   ofstream ougen, ouga;
   if (!genotypeIndex.empty()) {
      ougen.open(genofile.c_str());
      ouga.open(genoaafile.c_str());
   }
   ofstream oub(basefile.c_str());

   // output result
   //analyzer.reorderResult();
   printCodons(oucod);
   printConsensus(ouc);
   // the following only make sense if the genotypeIndex is not empty!
   if (!genotypeIndex.empty()) {
      printGenotype(ougen);
      printGenotypeAA(ouga);
   }
   printBases(oub);

   oucod.close();
   ouc.close();
   if (!genotypeIndex.empty()) {
      ougen.close();
      ouga.close();
   }
   oub.close();
   ofstream ouqual(qualfile.c_str());
   if (ouqual.fail()) {
      cerr << "Failed to open " << qualfile << endl;
      exit(1);
   }
   printQuality(ouqual);
   ouqual.close();
}

void GenotypeAnalyzer::setAlignOutputFile(const string &fname) {
   oaln.open(fname.c_str());
   if (oaln.fail()) {
      cerr << "Failed to open " << fname << endl;
      exit(1);
   }
}

void GenotypeAnalyzer::setRefseq(const DNA &ref) {
   if (refseq.empty()) { refseq = ref; }
   else {
      if (ref.length() != reflib.begin()->second.length()) {
         cerr << "you cannot switch to a reference sequence of different length!\n";
         exit(1);
      }
      refseq = ref;
   }
}

void GenotypeAnalyzer::printAlninfo(const string &alninfoFile) const {
   ofstream ous(alninfoFile.c_str());
   if (ous.fail()) {
      cerr << "Failed to open alninfo file: " << alninfoFile << endl;
      exit(1);
   }

   ous << "score\tidentity\tcoverage\treadlen\tcount\n";
   for (auto itr = alnsummary.cbegin(); itr != alnsummary.cend(); ++itr) {
      ous << itr->first << "\t" << itr->second << endl;
   }
   ous.close();
}

// static function to use a particular reference library
void GenotypeAnalyzer::loadReflib() {
   DNA dna;
   string header;
   ifstream ifs(reflibFile.c_str());
   if (ifs.fail()) {
      cerr << "failed to open " << reflibFile << endl;
      exit(1);
   }
   while (dna.read(ifs, header)) {
      if (dna.getName().find('|') != string::npos) {
         vector<string> nameparts = split(dna.getName(), '|');
         reflib.insert(make_pair(nameparts[3], dna));
      } 
      else { // handmade reference seq
         reflib.insert(make_pair(dna.getName(), dna));
      }
   }
   cout << reflib.size() << " references read into memory\n";
   //for (auto itr=reflib.begin(); itr != reflib.end(); ++itr) {
   //   cout << itr->first << "\n" << itr->second << endl;
   //}
}
string GenotypeAnalyzer::getRefseqShortName() const {
   string shortname;
   if (refseq.getName().find('|') != string::npos) {
      vector<string> parts = split(refseq.getName(), '|');
      if (parts.size() > 3) {
         shortname=parts[3];
      }
      else {
         cerr << "seqname " << refseq.getName() << " lack third | break\n";
         exit(1);
      }
   }
   else shortname = refseq.getName();
   size_t x;
   if ((x=shortname.find('.')) != string::npos) {
      shortname = shortname.substr(0, x);
   }
   return shortname;
}

string GenotypeAnalyzer::suckupReads(const string &file) {
   inputFile1=file;
   cerr << "Reading fastq into memory ...\n";
   ifstream inf(file.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << file << endl;
      exit(1);
   }
   if (inf.fail()) {
      cerr << "Failed to open input file: " << file << endl;
      return "";
   }
   Fastq fasq;
   int cnt=0;
   int total=0;
   while (fasq.read(inf)) {
      ++total;
      if (fasq.length() > alnlengthCut) {
         ++cnt;
         // remove the @ from fastq names
         //reads.push_back(DNAQual(fasq.getName().substr(1), fasq.getSequence(), fasq.getQuality()));
         reads.push_back(DNAQual(fasq.getName().substr(1), fasq.getSequence(), fasq.getQscore()));
         if (fasq.hasDescription()) {
            reads.back().setTitle(fasq.getDescription());
         }
      }
   }
   cerr << reads.size() << " Fastq sequences longer than " << alnlengthCut
      << " pulled into memory out of " << total << " from file: " << file << "\n";
   return getRefseqShortName();
}

/////////////////// paired version //////////////////
//
// helper function
/*
bool GenotypeAnalyzerPair::consumeOneSequence(const Fastq &read, ofstream &oufail) {
   DNA dna(read.getName(), read.getSequence());
   aligner.setSeq2(dna);
   aligner.runlocal();
   ++alnsummary[AlignInfo(aligner.getScore(), roundPercent(aligner.getIdentity()),
         roundPercent(aligner.getCov1()), dna.length())];
   if (aligner.getIdentity() > identityCut 
         && aligner.getSeq1AlignedLength() > alnlengthCut) 
   {
      accumulate(aligner.getTopAln(), aligner.getBottomAln(),
         (unsigned int)(aligner.topBeginIndex()));
      return true;
   }
   else if (read.length()> getExpandedAlnlengthCutoff()) {
      if (fillgap) {
         oufail << read;
      }
   }
#ifdef DEBUG
   if (fillgap) aligner.printAlign(oaln, 80);
#endif
   return false;
}

bool GenotypeAnalyzerPair::consumeOneSequenceRC(const Fastq &read, ofstream &oufail) {
   Fastq readrc(read);
   readrc.revcomp();
   return consumeOneSequence(readrc, oufail);
}
*/

map<string,int> GenotypeAnalyzerPair::consumePair() {
   if (!genotypeIndex.empty())
      cerr << genotypeIndex.size() << " position to look after\n";
   clear();
   int every = 100;
   size_t i;
   int goodcnt = 0, failcnt=0;
   string failFile = constructFailFile();
   ofstream oufail(failFile.c_str());
   map<string, int> matchDetail;
   char fdirec, bdirec;

   for (i=0; i < reads.size(); ++i) {
      // align forward direction
      fdirec=consumeOneBothDirections(reads[i]);
      if (fdirec == '?') {
         if (saveFailedRead(reads[i], oufail)) ++failcnt;
      }
      else {
         ++goodcnt;
      }
      updateAlninfo();
      // you may not find the second of the pair
      // align other end backward
      auto ptr=reads2.find(reads[i].getName());
      if (ptr == reads2.end()) {
         ++matchDetail[string(1,fdirec).append(1,'?')];
         continue; // not finding second for pair
      }
      bdirec='?';
      if (fdirec == '+') {
         if (consumeOneSequenceRC(ptr->second)) {
            bdirec='-';
         }
         else {
            if (consumeOneSequence(ptr->second)) {
               bdirec='+';
            }
         }
      }
      else if (fdirec == '-') {
         if (consumeOneSequence(ptr->second)) {
               bdirec='+';
         }
         else {
            if (consumeOneSequenceRC(ptr->second)) {
               bdirec='-';
            }
         }
      }
      else { // have to try both directions
         bdirec = consumeOneBothDirections(ptr->second);
      }
      if (bdirec == '?') {
         if (saveFailedRead(ptr->second, oufail)) ++failcnt;
      }
      else {
         ++goodcnt;
      }
      updateAlninfo();
      ++matchDetail[string(1,fdirec).append(1,bdirec)];
      if ((i+1)%every == 0) {
         cerr << "working on reads " << i+1 << " identity=" << roundPercent(aligner.getIdentity())
            << " refcov=" << roundPercent(aligner.getCov1()) << endl;
         every *= 2;
      }
   }

   // work on singles from the second of the pair
   cerr << "working on reads that don't have the first of the pair ...\n";
   vector<DNAQual*> single2=getSingles();
   for (i=0; i<single2.size(); ++i) {
      bdirec=consumeOneBothDirections(*single2[i]);
      if (bdirec == '?') {
         if (saveFailedRead(*(single2[i]), oufail)) ++failcnt;
      }
      else {
         ++goodcnt;
      }
      updateAlninfo();
      ++matchDetail[string("?").append(1,bdirec)];
   }
   cerr << reads.size() + single2.size() << " read pairs processed\n";
   cerr << failcnt << " longer than " << getExpandedAlnlengthCutoff() 
         << " and did not match to refseq well\n";
   cout << "matching detail\n";
   for (auto ptr=matchDetail.begin(); ptr != matchDetail.end(); ++ptr) {
      cout << ptr->first << " " << ptr->second << endl;
   }
   matchDetail["goodcnt"]=goodcnt;
   matchDetail["failcnt"]=failcnt;
   return matchDetail;
}

void GenotypeAnalyzerPair::suckupReads2(const string &file) {
   cerr << "suckupReads2() Reading fastq into map ...\n";
   ifstream inf(file.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << file << endl;
      exit(1);
   }
   if (inf.fail()) {
      cerr << "Failed to open input file: " << file << endl;
      return;
   }
   Fastq fasq;
   int cnt=0;
   int total=0;
   while (fasq.read(inf)) {
      ++total;
      if (fasq.length() > alnlengthCut) {
         ++cnt;
         // remove the @ from fastq names
         DNAQual tmp(fasq.getName().substr(1), fasq.getSequence(), fasq.getQuality());
         if (fasq.hasDescription()) {
            tmp.setTitle(fasq.getDescription());
         }
         reads2[tmp.getName()] = tmp;
      }
   }
   cerr << reads2.size() << " Fastq sequences longer than " << alnlengthCut
      << " pulled into memory out of " << total << " from file: " << file << "\n";
}

vector<DNAQual*> GenotypeAnalyzerPair::getSingles() {
   set<string> reads1name;
   for (auto i=0; i<reads.size(); ++i) {
      reads1name.insert(reads[i].getName());
   }
   vector<DNAQual*> tmp;
   for (auto ptr=reads2.begin(); ptr != reads2.end(); ++ptr) {
      if (reads1name.find(ptr->first) == reads1name.end())
         tmp.push_back(&(ptr->second));
   }
   return tmp;
}

