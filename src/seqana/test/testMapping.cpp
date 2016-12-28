#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include "MismatchMapping.h"

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

/*
bool diffByOne(const string &s1, const string &s2) {
   int diff = 0;
   for (size_t i=0; i<s1.size(); ++i) {
      if (s1[i] != s2[i]) {
         ++diff;
         if (diff > 1) return false;
      }
   }
   //cout << s1 << endl << s2 << " differ by one\n";
   return true;
}

bool diffByTwoOrMore(const string &s1, const string &s2) {
   int diff = 0;
   for (size_t i=0; i<s1.size(); ++i) {
      if (s1[i] != s2[i]) {
         ++diff;
         if (diff > 1) {
            //cout << s1 << endl << s2 << " differ by 2 or more\n";
            return true;
         }
      }
   }
   return false; // differe exactly by one
}

void printMissedOne(const map<pair<string, string>, int> &miss, ostream &ous) {
   map<pair<string, string>, int>::const_iterator it = miss.begin();
   while (it != miss.end()) {
      ous << it->first.first << '\t' << it->first.second << '\t' << it->second << endl;
      ++it;
   }
}


int mismatch1(map<string, int> &refcount, map<pair<string, string>, int> &missedOne) {
   map<string, int>::iterator mit, before, after, bb, aa;
   mit = refcount.begin();
   int corrections = 0;

   //map<pair<string, string>, int> missedOne;
   // check the beginning condition
   before = refcount.begin();
   ++mit;
   if (mit != refcount.end()) { // incase there is only two elements, this is impossible
      // but I am checking it for sure
      if (mit->second >= 0 && before->second < 0 && diffByOne(before->first, mit->first)) {
         missedOne[make_pair(before->first, mit->first)] = before->second*-1;
         mit->second -= before->second;
         refcount.erase(before);
         ++mit;
         ++corrections;
      }
   }
   int i = 0;
   while (mit != refcount.end()) {
      if (mit->second >= 0) {
         if (i > 1) { // make sure bb not out of bound
            before = mit;
            --before;
            bb = before;
            --bb;
            if (before->second < 0 && diffByOne(before->first, mit->first)
                  && (bb->second < 0 || diffByTwoOrMore(bb->first, before->first)) ) {
               missedOne[make_pair(before->first, mit->first)] = before->second*-1;
               mit->second -= before->second;
               refcount.erase(before);
               ++corrections;
            }
         }
         after = mit;
         ++after;
         if (after == refcount.end()) break;
         aa = after;
         ++aa;
         if (aa == refcount.end()) break;
         if (after->second < 0 && diffByOne(after->first, mit->first)
               && (aa->second < 0 || diffByTwoOrMore(aa->first, after->first))) {
            missedOne[make_pair(after->first, mit->first)] = after->second*-1;
            mit->second -= after->second;
            refcount.erase(after);
            ++corrections;
         }
      }
      ++mit; ++i;
   }
   // end condition
   mit = refcount.end();
   --mit;
   after = mit;
   --mit;
   if (mit->second >= 0 && after->second < 0 && diffByOne(mit->first, mit->first)) {
      missedOne[make_pair(after->first, mit->first)] = after->second*-1;
      mit->second -= after->second;
      refcount.erase(after);
      ++corrections;
   }
   return corrections;
}

void summarizeMapping(const map<string, int> &refc, map<string, int> &result) {
   map<string, int>::const_iterator mit = refc.begin();
   int nomatch=0, match=0; // match species
   // nomatch is the library sequence got no match
   int nomatchInput = 0;
   int nomatchInputSum = 0;
   int matchSum = 0;
   while (mit != refc.end()) {
      if (mit->second > 0) {
         ++match;
         matchSum += mit->second;
      }
      else if (mit->second == 0) ++nomatch;
      else { // negative number
         ++nomatchInput;
         nomatchInputSum += mit->second;
      }
      ++mit;
   }
   result.clear();
   result["nomatch"] = nomatch;
   result["match"] = match;
   result["nomatchInput"] = nomatchInput;
   result["nomatchInputSum"] = nomatchInputSum;
   result["matchSum"] = matchSum;
}

////////////////// reverse direction //////////////////////

class ReverseSort {
   public:
      ReverseSort() { }
      bool operator()(const string &s1, const string &s2) const {
         //cout << "comparing\n" << s1 << endl << s2 << endl;
         for (int i = s1.size()-1; i > -1; --i) {
            //cout << s1[i] << " x " << s2[i] << endl;
            if (s1[i] < s2[i]) { 
               //cout << s1 << " before " << s2 << endl;
               return true; 
            }
            else if (s1[i] > s2[i]) {
               //cout << s2 << " before " << s1 << endl;
               return false;
            }
         }
         return false;
      }
};

bool diffByOneReverse(const string &s1, const string &s2) {
   int diff = 0;
   for (int i=s1.size()-1; i > -1; --i) {
      if (s1[i] != s2[i]) {
         ++diff;
         if (diff > 1) return false;
      }
   }
   //cout << s1 << endl << s2 << " differ by one\n";
   return true;
}

bool diffByTwoOrMoreReverse(const string &s1, const string &s2) {
   int diff = 0;
   for (int i=s1.size()-1; i>-1; --i) {
      if (s1[i] != s2[i]) {
         ++diff;
         if (diff > 1) {
            //cout << s1 << endl << s2 << " differ by 2 or more\n";
            return true;
         }
      }
   }
   return false; // differe exactly by one
}

int mismatch1Reverse(map<string, int, ReverseSort> &refcount, map<pair<string, string>, int> &missedOne) {
   map<string, int, ReverseSort>::iterator mit, before, after, bb, aa;
   mit = refcount.begin();
   int corrections = 0;

   //map<pair<string, string>, int> missedOne;
   // check the beginning condition
   before = refcount.begin();
   ++mit;
   if (mit != refcount.end()) { // incase there is only two elements, this is impossible
      // but I am checking it for sure
      if (mit->second >= 0 && before->second < 0 && diffByOneReverse(before->first, mit->first)) {
         missedOne[make_pair(before->first, mit->first)] = before->second*-1;
         mit->second -= before->second;
         refcount.erase(before);
         ++mit;
         ++corrections;
      }
   }
   int i = 0;
   while (mit != refcount.end()) {
      if (mit->second >= 0) {
         if (i > 1) { // make sure bb not out of bound
            before = mit;
            --before;
            bb = before;
            --bb;
            if (before->second < 0 && diffByOneReverse(before->first, mit->first)
                  && (bb->second < 0 || diffByTwoOrMoreReverse(bb->first, before->first)) ) {
               missedOne[make_pair(before->first, mit->first)] = before->second*-1;
               mit->second -= before->second;
               refcount.erase(before);
               ++corrections;
            }
         }
         after = mit;
         ++after;
         if (after == refcount.end()) break;
         aa = after;
         ++aa;
         if (aa == refcount.end()) break;
         if (after->second < 0 && diffByOneReverse(after->first, mit->first)
               && (aa->second < 0 || diffByTwoOrMoreReverse(aa->first, after->first))) {
            missedOne[make_pair(after->first, mit->first)] = after->second*-1;
            mit->second -= after->second;
            refcount.erase(after);
            ++corrections;
         }
      }
      ++mit; ++i;
   }
   // end condition
   mit = refcount.end();
   --mit;
   after = mit;
   --mit;
   if (mit->second >= 0 && after->second < 0 && diffByOneReverse(mit->first, mit->first)) {
      missedOne[make_pair(after->first, mit->first)] = after->second*-1;
      mit->second -= after->second;
      refcount.erase(after);
      ++corrections;
   }
   return corrections;
}

void summarizeMappingReverse(const map<string, int, ReverseSort> &refc, map<string, int> &result) {
   map<string, int, ReverseSort>::const_iterator mit = refc.begin();
   int nomatch=0, match=0; // match species
   // nomatch is the library sequence got no match
   int nomatchInput = 0;
   int nomatchInputSum = 0;
   int matchSum = 0;
   while (mit != refc.end()) {
      if (mit->second > 0) {
         ++match;
         matchSum += mit->second;
      }
      else if (mit->second == 0) ++nomatch;
      else { // negative number
         ++nomatchInput;
         nomatchInputSum += mit->second;
      }
      ++mit;
   }
   result.clear();
   result["nomatch"] = nomatch;
   result["match"] = match;
   result["nomatchInput"] = nomatchInput;
   result["nomatchInputSum"] = nomatchInputSum;
   result["matchSum"] = matchSum;
}

template<class T1, class T2, class cmp> 
void printMapComp(const map<T1, T2, cmp> &m, ostream &ous) {
   typename map<T1, T2, cmp>::const_iterator it = m.begin();
   while (it != m.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}
*/

int main(int argc, char* argv[]) {

   /*
   ifstream inf("refcount.tab");
   if (!inf) {
      cerr << "Failed to open refcount.tab" << endl;
      return 1;
   }

// loading map from file
   map<string, int> orig;
   string refs;
   int count;
   inf >> refs >> count;
   while (!inf.eof()) {
      orig[refs] = count;
      inf >> refs >> count;
   }

   MismatchMapping mm(orig);
   mm.run();

   mm.printFinalSummary(cout);
   mm.printMappingToFile("refcount_missed1.tab");
   */

   map<char, int> chcnt;
   ++chcnt['a'];
   ++chcnt['a'];
   ++chcnt['b'];
   ++chcnt['t'];
   ++chcnt['t'];
   for (auto itr=chcnt.begin(); itr != chcnt.end(); ++itr) {
      cout << itr->first << "\t" << itr->second << endl;
   }

   return 0;
}
