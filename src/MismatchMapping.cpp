#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include "MismatchMapping.h"
#include <vector>

using namespace std;

/**
 * differ exactly by one. In the contaxt of unique string.
 */
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
/*
// this can be replaced with generic printMap template function
void printMissedOne(const map<pair<string, string>, int> &miss, ostream &ous) {
   map<pair<string, string>, int>::const_iterator it = miss.begin();
   while (it != miss.end()) {
      ous << it->first.first << '\t' << it->first.second << '\t' << it->second << endl;
      ++it;
   }
}
*/


int MismatchMapping::forwardRun() {
   map<string, int>::iterator mit, before, after, bb, aa;
   mit = refcnt.begin();
   int corrections = 0;

   //map<pair<string, string>, int> missOne;
   // check the beginning condition
   before = refcnt.begin();
   ++mit;
   if (mit != refcnt.end()) { // incase there is only two elements, this is impossible
      // but I am checking it for sure
      if (mit->second >= 0 && before->second < 0 && diffByOne(before->first, mit->first)) {
         missOne[make_pair(before->first, mit->first)] = before->second*-1;
         mit->second -= before->second;
         refcnt.erase(before);
         ++mit;
         ++corrections;
      }
   }
   int i = 0;
   while (mit != refcnt.end()) {
      if (mit->second >= 0) {
         if (i > 1) { // make sure bb not out of bound
            before = mit;
            --before;
            bb = before;
            --bb;
            if (before->second < 0 && diffByOne(before->first, mit->first)
                  && (bb->second < 0 || diffByTwoOrMore(bb->first, before->first)) ) {
               missOne[make_pair(before->first, mit->first)] = before->second*-1;
               mit->second -= before->second;
               refcnt.erase(before);
               ++corrections;
            }
         }
         after = mit;
         ++after;
         if (after == refcnt.end()) break;
         aa = after;
         ++aa;
         if (aa == refcnt.end()) break;
         if (after->second < 0 && diffByOne(after->first, mit->first)
               && (aa->second < 0 || diffByTwoOrMore(aa->first, after->first))) {
            missOne[make_pair(after->first, mit->first)] = after->second*-1;
            mit->second -= after->second;
            refcnt.erase(after);
            ++corrections;
         }
      }
      ++mit; ++i;
   }
   // end condition
   mit = refcnt.end();
   --mit;
   after = mit;
   --mit;
   if (mit->second >= 0 && after->second < 0 && diffByOne(mit->first, mit->first)) {
      missOne[make_pair(after->first, mit->first)] = after->second*-1;
      mit->second -= after->second;
      refcnt.erase(after);
      ++corrections;
   }
   return corrections;
}

void MismatchMapping::summarize() const {
   map<string, int>::const_iterator mit = refcnt.begin();
   int nomatch=0, match=0; // match species
   // nomatch is the library sequence got no match
   int nomatchInput = 0;
   int nomatchInputSum = 0;
   int matchSum = 0;
   while (mit != refcnt.end()) {
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
   //summary.clear();
   summary["nomatch"] = nomatch;
   summary["match"] = match;
   summary["nomatchInput"] = nomatchInput;
   summary["nomatchInputSum"] = nomatchInputSum;
   summary["matchSum"] = matchSum;
}

////////////////// reverse direction //////////////////////


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

int MismatchMapping::backwardRun() {
   map<string, int, ReverseSort>::iterator mit, before, after, bb, aa;
   mit = refcntrev->begin();
   int corrections = 0;

   // check the beginning condition
   before = refcntrev->begin();
   ++mit;
   if (mit != refcntrev->end()) { // incase there is only two elements, this is impossible
      // but I am checking it for sure
      if (mit->second >= 0 && before->second < 0 && diffByOneReverse(before->first, mit->first)) {
         missOne[make_pair(before->first, mit->first)] = before->second*-1;
         mit->second -= before->second;
         refcntrev->erase(before);
         ++mit;
         ++corrections;
      }
   }
   int i = 0;
   while (mit != refcntrev->end()) {
      if (mit->second >= 0) {
         if (i > 1) { // make sure bb not out of bound
            before = mit;
            --before;
            bb = before;
            --bb;
            if (before->second < 0 && diffByOneReverse(before->first, mit->first)
                  && (bb->second < 0 || diffByTwoOrMoreReverse(bb->first, before->first)) ) {
               missOne[make_pair(before->first, mit->first)] = before->second*-1;
               mit->second -= before->second;
               refcntrev->erase(before);
               ++corrections;
            }
         }
         after = mit;
         ++after;
         if (after == refcntrev->end()) break;
         aa = after;
         ++aa;
         if (aa == refcntrev->end()) break;
         if (after->second < 0 && diffByOneReverse(after->first, mit->first)
               && (aa->second < 0 || diffByTwoOrMoreReverse(aa->first, after->first))) {
            missOne[make_pair(after->first, mit->first)] = after->second*-1;
            mit->second -= after->second;
            refcntrev->erase(after);
            ++corrections;
         }
      }
      ++mit; ++i;
   }
   // end condition
   mit = refcntrev->end();
   --mit;
   after = mit;
   --mit;
   if (mit->second >= 0 && after->second < 0 && diffByOneReverse(mit->first, mit->first)) {
      missOne[make_pair(after->first, mit->first)] = after->second*-1;
      mit->second -= after->second;
      refcntrev->erase(after);
      ++corrections;
   }
   return corrections;
}

void MismatchMapping::run() { 
   forward(); 
   refcntrev = new map<string, int, ReverseSort>(refcnt.begin(), refcnt.end());
   backward(); 
}

void MismatchMapping::summarizeReverse() const {
   map<string, int, ReverseSort>::const_iterator mit = refcntrev->begin();
   int nomatch=0, match=0; // match species
   // nomatch is the library sequence got no match
   int nomatchInput = 0;
   int nomatchInputSum = 0;
   int matchSum = 0;
   //map<int,int>::iterator hi;
   while (mit != refcntrev->end()) {
      ++numhitDist[mit->second];
      if (mit->second > 0) {
         ++match;
         matchSum += mit->second;
      }
      else if (mit->second == 0) {
         ++nomatch;
      }
      else { // negative number
         ++nomatchInput;
         nomatchInputSum += mit->second;
      }
      ++mit;
   }
   //summaryrev.clear();
   summaryrev["nomatch_miss1"] = nomatch;
   summaryrev["match_miss1"] = match;
   summaryrev["nomatchInput_miss1"] = nomatchInput;
   summaryrev["nomatchInputSum_miss1"] = nomatchInputSum;
   summaryrev["matchSum_miss1"] = matchSum;
}

/* you don't need this one, hard-coding field name to the summaryReverse
 * function.
void MismatchMapping::addSuffixToFinalField(const string &sfx) const {
   map<string, int> tmp;
   map<string, int>::const_iterator it = summaryrev.begin();
   while (it != summaryrev.end()) {
      tmp[it->first + sfx] = it->second;
      ++it;
   }
   summaryrev = tmp;
}
*/
   
void MismatchMapping::forward() {
   /*
   int round = 1;
   int num;
   while ((num=forwardRun()) > 0) {
      cout << "Forward round " << round << " updated " << num << " mismatches\n";
      ++round;
   }
   */
   while (forwardRun() > 0);
   summarize();
}

void MismatchMapping::backward() {
   /*
   int round = 1;
   int num;
   while ((num=backwardRun()) > 0) {
      cout << "Backward round " << round << " updated " << num << " mismatches\n";
      ++round;
   }
   */
   while (backwardRun() > 0);
   summarizeReverse();
}


void MismatchMapping::printForward(ostream &ous) const {
   map<string, int>::const_iterator it = refcnt.begin();
   while (it != refcnt.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}

void MismatchMapping::printBackward(ostream &ous) const {
   map<string, int, ReverseSort>::const_iterator it = refcntrev->begin();
   while (it != refcntrev->end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}

void MismatchMapping::printSummary(ostream &ous) const {
   if (summary.empty()) {
      summarize();
   }
   map<string, int>::const_iterator it = summary.begin();
   while (it != summary.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
}

void MismatchMapping::printReverseSummary(ostream &ous) const {
   if (summaryrev.empty()) {
      summarizeReverse();
   }
   // for defined order not affected by string sorting order.
   vector<string> fields;
   fields.push_back("match_miss1");
   fields.push_back("matchSum_miss1");
   fields.push_back("nomatch_miss1");
   fields.push_back("nomatchInput_miss1");
   fields.push_back("nomatchInputSum_miss1");
   for (size_t i=0; i<fields.size(); ++i) {
      ous << fields[i] << '\t' << summaryrev[fields[i] ] << endl;
   }

   /*
   map<string, int>::const_iterator it = summaryrev.begin();
   while (it != summaryrev.end()) {
      ous << it->first << '\t' << it->second << endl;
      ++it;
   }
   */
}

void MismatchMapping::printForwardToFile(const string &ouf) const {
   ofstream of(ouf.c_str());
   if (of.fail()) {
      cerr << "Failed to open " << ouf << endl;
      throw runtime_error("failed to open file: " + ouf);
   }
   printForward(of);
   of.close();
}

void MismatchMapping::printBackwardToFile(const string &ouf) const {
   ofstream of(ouf.c_str());
   if (of.fail()) {
      cerr << "Failed to open " << ouf << endl;
      throw runtime_error("failed to open file: " + ouf);
   }
   printBackward(of);
   of.close();
}

void MismatchMapping::printSummaryToFile(const string &ouf) const {
   ofstream of(ouf.c_str());
   if (of.fail()) {
      cerr << "Failed to open " << ouf << endl;
      throw runtime_error("failed to open file: " + ouf);
   }
   printSummary(of);
   of.close();
}

void MismatchMapping::printReverseSummaryToFile(const string &ouf) const {
   ofstream of(ouf.c_str());
   if (of.fail()) {
      cerr << "Failed to open " << ouf << endl;
      throw runtime_error("failed to open file: " + ouf);
   }
   printReverseSummary(of);
   of.close();
}

void MismatchMapping::printDistribution(ostream &ous) const {
   map<int, int>::const_iterator it = numhitDist.begin();
   while (it != numhitDist.end()) {
      ous << it->first << "\t" << it->second << endl;
      ++it;
   }
}
