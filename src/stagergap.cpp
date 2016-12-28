#include "stagergap.h"
#include <iostream>

int Stagergap::wdsz=12;

void Stagergap::display() {
   cout << top << endl << bottom << endl;
   cout << bgi << endl;
}

// ==
// CATGCCAACGTGGGTACAAGGG-AGTCTGGCGGGGAGATGGCGTCATGCAC
// |||||||||| ||||||| ||| ||||||||||| |||||| |||||||||
// CATGCCAACGAGGGTACATGGGGAGTCTGGCGGG-AGATGGTGTCATGCAC
// window width = 12
int Stagergap::left() {
   if (bgi < 3) return -1;
   if (bgi+1 < top.length() && bottom[bgi+1] == '-') 
      return -1;
   int lbound = (int)bgi-wdsz;
   if (lbound < 0) lbound = 0;
   int x = bgi;
   while (x > lbound) { // check the first found gap in top
      if (top[x] == '-') {
         if (bottom[x] != 'X' && isupper(bottom[x]) 
               && x-1 > 0 && top[x-1] != '-') {
            //cerr << "detected left stager at " << x << " and " << bgi << endl;
            //display();
            return x;
         }
         else {
            return -1;
         }
      }
      --x;
   }
   return -1;
}
//
//CATGAGGATCGTTGGGCCTAAAACCTGCAGCGTGGC-ACGGAACATTCCCCAT
//|||||||||||||||||||||| ||||||||||||| ||||||||||||||||
//CATGAGGATCGTTGGGCCTAAA-CCTGCAGCGTGGCCACGGAACATTCCCCAT
//
// return -1 for not found, > 0 for the location found in top
int Stagergap::right() {
   if (bgi+1 < top.length() && bottom[bgi+1] == '-')
      return -1;
   size_t rbound = bgi + wdsz;
   if (rbound >= top.size()) rbound=top.size()-1;
   if (rbound <= bgi+3) return -1;
   size_t x = bgi;
   while (x <= rbound) {
      if (top[x] == '-') {
         if (bottom[x] != 'X' && isupper(bottom[x])
               && x+1 < rbound && top[x+1] != '-') {
            return x;
         }
         else return -1;
      }
      ++x;
   }
   return -1;
}
///  topi
// AC-GTACCGCA  ACGTACCGCA   more different, no change
// ACCGTAC-GCA  ACCGTACGCA
//        bgi
void Stagergap::shiftLeft(size_t topi) {
   //cerr << "left gap at " << topi << endl;
   top.erase(topi, 1);
   bottom.erase(bgi, 1);
}

///  topi
// AC-GTACCGCA  AC-GTACCGCA  more similar, two changes, if error greater
// ACCGTAC-GCA  ACXGTACcGCA  X was used for postprocessing.
//        bgi
void Stagergap::fixLeft(size_t topi) {
   //cerr << "fixing left stager: " << topi << " " << bottom[topi]
   //   << " " << bgi << top[bgi] << endl;
   //display();
   bottom[topi]='X';
   bottom[bgi]=tolower(top[bgi]);
   //display();
}

///        bgir
// ACCAGTAC-CA  ACCAGTACCA   more different
// ACC-GTACGCA  ACCGTACGCA
//    bgi
void Stagergap::shiftRight(size_t topi) {
   top.erase(topi, 1);
   bottom.erase(bgi, 1);
}

///        bgir
// ACCAGTAC-CA  ACCAGTAC-CA   more different
// ACC-GTACGCA  ACCaGTACXCA
//    bgi
void Stagergap::fixRight(size_t topi) {
   bottom[topi]='X';
   bottom[bgi]=tolower(top[bgi]);
}

bool Stagergap::findAndFixLeft() {
   int i = left();
   if (i != -1) {
      fixLeft(i);
      return true;
   }
   return false;
}
bool Stagergap::findAndFixRight() {
   int i = right();
   if (i != -1) {
      fixRight(i);
      return true;
   }
   return false;
}

