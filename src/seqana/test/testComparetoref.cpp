#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <stack>
#include <algorithm>
#include <fstream>
#include <strformat.h>
#include <iterator>

using namespace std;

// correction algorithm for poor alignment
// or large block mutation.
void shrinkGap(string &top, string &bottom) {
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
   }
   // to the erasing
   while (!eraseIndex.empty()) {
      pair<unsigned int, unsigned int> ij = eraseIndex.top();
      top.erase(ij.first, 1);
      bottom.erase(ij.second, 1);
      eraseIndex.pop();
   }
}
         
void tabulateCodons(vector<map<string, int> > &result, map<string, int> &target, 
      string &top, string &bottom, unsigned int tb=0) 
{
   shrinkGap(top, bottom);
   string::size_type i, refi;
   i=0;
   refi=0; // reference base index 
   if (tb > 0) {
      refi = tb;
      while (refi % 3 != 0) {
         ++i;
         ++refi;
      }
   }
   cout << top.length() << " " << bottom.length() << endl;
   // target codons in original 28, 30, 31, 32, 54, 58, 62, 92, 93
   // target codons in my refseq 18, 20, 21, 22, 44, 48, 52, 82, 83
   // the index should be pos - 1
   // These are protein index we need nucleotid index
   //static vector<int> targetCodonIndex{17, 19, 20, 21, 43, 47, 51, 81, 82};
   //*3
   static vector<int> targetCodonIndex{51, 57, 60, 63, 129, 141, 153, 243, 246};

   string codon(3, 'x');
   string combined;
   unsigned int t = 0;
   while (i < top.length()) {
      //cout << i << " " << refi << endl;
      if (top[i] != '-') {
         codon[refi%3] = bottom[i];
         if (refi%3 == 2) {
            // if target codon reached
            if (refi-2 == targetCodonIndex[t]) {
               combined += codon;
               ++t;
            }
            
            ++(result[refi/3][codon]);
         }
         ++refi;
      }
      ++i;
   }
   //cout << combined << endl;
   if (combined.size() == targetCodonIndex.size()*3) {
      ++target[combined];
   }
   else {
      cerr << "combined target codon not the expected size: " << combined
         << " should be " << targetCodonIndex.size()*3 << " long!\n";
   }
}

void showMapByCount(const map<string, int> &src) {
   vector<pair<string,int> > elems;
   transform(src.begin(), src.end(), back_inserter(elems), 
         [](const pair<string,int> &p) { 
            return p; 
         }
   );
   for (unsigned int i=0; i<elems.size(); ++i) {
      cout << elems[i].first << " " << elems[i].second << endl;
   }

   sort(elems.begin(), elems.end(), 
         [](const pair<string, int>& l, const pair<string,int>& r) {
            return l.second > r.second;
         }
   );
   cout << "after sorting\n";
   for (unsigned int i=0; i<elems.size(); ++i) {
      cout << elems[i].first << " " << elems[i].second << endl;
   }
}

class CountSorter {
   public:
      bool operator()(const pair<string,int>&p1, const pair<string,int>&p2) {
         return p1.second > p2.second;
      }
};

vector<pair<string, int> > map2vectorByCount(const map<string, int>& src) {
   vector<pair<string,int> > arr;
   for (auto itr=src.begin(); itr != src.end(); ++itr) {
      arr.push_back(pair<string, int>(itr->first, itr->second));
   }   
   cout << "before sort:\n";
   for (size_t i=0; i<arr.size(); ++i) {
         cout << arr[i].first << " -> " << arr[i].second <<  " ";
            }   
   cout << endl;
   CountSorter sorter;
   sort(arr.begin(), arr.end(), sorter);
   cout << "after sort:\n";
   for (size_t i=0; i<arr.size(); ++i) {
                     cout << arr[i].first << " -> " << arr[i].second <<  " ";
   }   
   cout << endl;
   return arr; 
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

bool gapInPoly(const string &top, const string &bottom, size_t i) {
   if (i-1 < 0) return false;
   if (top[i] != bottom[i-1])  return false; 
   if (top[i-1] != bottom[i-1]) return false;
   return true;
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

void correctGap(string &top, string &bottom) {
   static const unsigned int gapw = 6;
   string::size_type i, j;
   stack<pair<unsigned int, unsigned int> > eraseIndex;
   
   i=5;
   cerr << "i-7=" << i-7 << " (int)(i-7)=" << (int)(i-7) << endl;

   for (i=0; i<top.size(); ++i) {
      if (bottom[i] == '-') {
         if (topNoGap(top, i, gapw)) {
            if (gapInPoly(top, bottom, i)) {
               bottom[i] = top[i];
               if (i+1 < top.size() && bottom[i+1] == '-') {
                  bottom[i+1] = top[i+1];
               }
            }
            // AAAAATGG
            // AAAAC-GG
            else if (((int)i-5) > -1 && top.substr(i-5, 5) == "AAAAA"
                  && bottom.substr(i-5,4) == "AAAA" && bottom[i-1] != 'A') {
               bottom[i] = bottom[i-1];
               bottom[i-1] = 'A';
            }
            // TTCTTCTCATGCC
            // TCATTC--ATGCC
            else if ((i-3) >= 0 && i+1 < top.length() 
                  && top.substr(i-3, 5) == "TTCTC" 
                  && bottom[i+1] == '-') {
               bottom[i]=top[i];
               bottom[i+1] = top[i+1];
            }
            else if (((int)i-1) > -1 && (i+1) < top.length() 
                  && bottom[i-1] != '-' && bottom[i+1]  != '-') {
               bottom[i] = top[i];
            }
            else if (top[i] == 'T' && i+1 < top.size() && top[i+1] == 'C'
                  && bottom[i+1] == '-') {
               bottom[i] = top[i];
               bottom[i+1] = top[i+1];
            }
            else if (top[i] == 'G' && i+1 < top.size() && top[i+1] == 'A'
                  && bottom[i+1] == '-') {
               bottom[i] = top[i];
               bottom[i+1] = top[i+1];
            }
            else if (top[i] == 'T' && i+1 < top.size() && top[i+1] == 'A'
                  && bottom[i+1] == '-') {
               bottom[i] = top[i];
               bottom[i+1] = top[i+1];
            }
            //else {
            //   cout << "Gap not fixed\n";
            //   displayGap(top, bottom, i, gapw);
            //}

         }
         else {
            if ((
                  ((i-2 > 0 && bottom[i-1] == bottom[i-2])
                        || (i+2 < bottom.length() && bottom[i+1] == bottom[i+2]))
                  && ((i>0 & top[i] == bottom[i-1]) 
                     || (i<bottom.length()-1 && top[i] == bottom[i+1]))
               )
               || (i-2>0 && bottom[i-2] == 'C' && bottom[i-1] == 'A'
                  && top[i-2] == 'C' && top[i-1] == 'A' && top[i] == 'A')
               ) {
               bottom[i] = top[i];
            }
            //else {
            //   displayGap(top, bottom, i, gapw);
            //}
         }
      }
   }
}
         
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
            cout << stoi(tmp[i]) << " ";
            result.push_back(stoi(tmp[i]));
         }
      }
      getline(inf, ln);
   }
   cout << endl;
   inf.close();
   cout << "numbers stored in file:\n";
   copy(result.begin(), result.end(), ostream_iterator<int>(cout, " "));
   cout << endl;
   return result;
}

int main(int argc, char* argv[]) {
      string top="TGGATATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAGTCCCCTTCTTCTCATGTCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCACCGGACATGTGAAAAACGGTTCCATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCACGCCCT--CCCCGGCGCCAAATTATTCTAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTACGTGG-AGGTTACGCGG-GTGGGGGAT";
      string bottom = "TGGATATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAGCTCTTGCCACGGTTGCCGGGAGTCCCCTTCTTCTCGTGCCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATTATGCAAACCACCTGCCCATGTGGAGCGCAAATCACCGGACATGTCAAAAACGGTTCCATGAGGATCGTCGGGCCTAGGACCTGCAGCAACACGTGGCATGGAACATTCCCCATCAACGCATACACCACGGGCCCCTGCAC-CCCTTCCCCCAGCGCCCAACTATTCCAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTATGTGGGAGGTTACGCGGAGTGGGGGAT";

   vector<map<string, int> > result(top.length()/3);
   map<string, int> target;
   tabulateCodons(result, target, top, bottom);
   top="TGGATATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAGTCCCCTTCTTCTCATGTCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCACCGGACATGTGAAAAACGGTTCCATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCACGCCCT-CCCCGGCGCCAAATTATTCTAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTACGTGGAGGTTACGCGGGTGGGGGATTTCCACTA";
   bottom="TGGATATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAACTCTTGCCACGGTTGTCGGGAGTCCCCTTCCTCTCGTGCCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATTATGCAAACCACCTGCCCATGTGGAGCGCAAATCACCGGCCATGTCAAAAACGGTTCCATGAGAATCGTCGGGCCTAGGACCTGCAGCAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCAC-CCCTTCCCCGGCGCCAAACTATTCCAGGGCGCTGTGGCGGGTAGCTGCTGAGGAGTATGTGGAGGTTACGCGAGTGGGGGATTTCCACTA";
   if (result.size() < top.length()/3)  {
      cout << "expanding result container\n";
      result.resize(top.length()/3);
   }
   tabulateCodons(result, target, top, bottom);
   top="TGGATATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAGTCCCCTTCTTCTCATGTCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCACC-GGACATGTGAAAAACGGTTCCATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCACGCCCT--CCCCGGC-GCCAAATTATTCTAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTACGTGGAGGTT-ACGCGGGTGGGGGATTTCCACTA";
   bottom="TGGATATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAGCTCTTGCCACGGTTGCCGGGAGTCCCCTTCTTCTCGTGCCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATTATGCAAACCACCTGCCCATGTGGAGCGCAAATCACCCGGACATGTCAAAAACGGTTCCATGAGGATCGTCGGGCCTAGGACCTGCAGCAACACGTGGCATGGAACATTCCCCATCAACGCGTACACCACGGGCCCCTGCAC-CCCTTCCCCCAGCCGCCCAACTATTCCAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTATGTGGAGGTTTACGCGTGTGGGGGATTTCCACTA";
   tabulateCodons(result, target, top, bottom);

   top="ATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAGTCCCCTTCTTCTCATGTCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCACCGGACATGTGAAAAACGGTTCCA-TGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCACGCCC-TCCCCGGCGCCAAATTATTCTAGGGCGCTGTGGCGGGTGGCTGCTGAGGAGTACGTGGAGGTTACGCG-GGTGGGGGATTTCCACTA";
   bottom="ATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAACTCTTGCCACGGTTGCCGGGAGTTCCTTTCTTCTCGTGCCAACGCGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAAATCACCGGACATGTCAAAAACGGTTCCAATGAGAATCGTCGGGCCTAGGGCCTGCAGCAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCAC-CCCATCCCCGGCGCCAAACTATTCCAGGGCGCTGTGGCGGGTAGCTGCTGAGGAGTATGTGGAGGTTACGCGAGGTGGGGGATTTCCACTA";
   tabulateCodons(result, target, top, bottom, 5); 

      top="TGGATATGCACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGG-AGTCCCCTTCTTCTCATGTCAACGTGGGTACAAGGGAGTCTGGCGGGGCGACGGCATCATGCAAACCACCTGCCCATGTGGAGCACAGATCACCGGACATGTGAAAAACGGTTCC-ATGAGGATCGTGGGGCCTAGGACCTGTAGTAACACGTGGCATGGAACATTCCCCATTAACGCGTACACCACGGGCCCCTGCACGCCCT--CCCCGGCGCCAAATTATTCTAGGGCGCTGTGGCGGGTGGCTGC-TGAGGAGTACGTGGAGGTTACGCGGGTGGGGGATTTCCACTA";
   bottom="TGGATATGCACGGTGTTGACTGACTTCAAGACCTGGCTCCAGTCCAAACTCTTGCCACGGTTGCCGGGGAGTCCCCTTCTTCTCGTGCCAACGCGGGTATAAGGGAGTCTGGCGGGGCGACGGCATTATGCAAACCACCTGCCCATGTGGAGCGCAAATCACCGGACATGTCAAAAACGGTTCCGATGAGGATCGTCGGGCCTAGGACCTGCAGCAACACGTGGCATGGAACATTCCCCATCAACGCATACACCACGGGCCCCTGCAC-CCCTTCCCCCAGCGCCCAACTATTCTAGGGCGCTGTGGCGGGTGGCTGCCTGAGGAGTATGTGGAGGTTACGCGAGTGGGGGATTTCCACTA";
   if (result.size() < top.length()/3)  {
      cout << "expanding result container\n";
      result.resize(top.length()/3);
   }
   tabulateCodons(result, target, top, bottom); 
      
   // output
   for (auto i = 0; i<result.size(); ++i) {
      cout << i << ": ";
      for (auto mi = result[i].begin(); mi != result[i].end(); ++mi) {
         cout << mi->first << " " << mi->second << "  ";
      }
      cout << endl;
   }
   for (auto i = target.begin(); i != target.end(); ++i) {
      cout << i->first << " " << i->second << endl;
   }
   map<string, int> ttt{{"first", 100}, {"second", 99},
      {"blala", 50}, {"foo", 30}, {"dog", 1}, {"cat", 3}};
   showMapByCount(ttt);
   vector<pair<string,int> > xx = map2vectorByCount(ttt);

   top=       "AGTCCCCTTCTTCTCATGCCAACGTGGGTACAA";
   string bot("AGTTCCTTCCTTC--ATGTCAACGTGGGTACAA"); 
   cout << top << endl << bot << endl;
   correctGap(top, bot);
   cout << "after correction\n" << top << endl << bot << endl;

   top= "GGTGTTGACTGACTTCAAGACC";
   bot= "GGTGTTGACCGA-TTCAAGACC"; 
   cout << top << endl << bot << endl;
   correctGap(top, bot);
   cout << "after correction\n" << top << endl << bot << endl;

   // AAAAATGG
   // AAAAC-GG
   top= "TTGACAAAAATGGTGATTCAAGA";
   bot= "TTGACAAAAC-GGCGATTCAAGA"; 
   cout << top << endl << bot << endl;
   correctGap(top, bot);
   cout << "after correction\n" << top << endl << bot << endl;

   // TTCTTCTCATGCC
   // TCATTC--ATGCC
   top= "GACAAAAATTCTTCTCATGCCGGTGATT";
   bot= "GACAAAACTCATTC--ATGCCGGCGATT"; 
   cout << top << endl << bot << endl;
   correctGap(top, bot);
   cout << "after correction\n" << top << endl << bot << endl;

   vector<int> tidx=readTargetIndex("comparetoIndex.txt");

   return 0;
}


