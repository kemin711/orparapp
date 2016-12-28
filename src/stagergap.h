#ifndef STAGERGAP_H
#define STAGERGAP_H

#include <string>

using namespace std;

class Stagergap {
   private:
      string& top;
      string& bottom;
      /** gap index in the bottom */
      size_t bgi; 
      /** gap index in the top withing wdsz distance */
      int tgi; // gap index in the top, -1 for none
      /** window size */
      static int wdsz; // default 12

   public:
      Stagergap(string &t, string &b, size_t idxb) : top(t), bottom(b), bgi(idxb) { }
      int left();
      //
      //CATGAGGATCGTTGGGCCTAAAACCTGCAGCGTGGC-ACGGAACATTCCCCAT
      //|||||||||||||||||||||| ||||||||||||| ||||||||||||||||
      //CATGAGGATCGTTGGGCCTAAA-CCTGCAGCGTGGCCACGGAACATTCCCCAT
      //
      // return -1 for not found, > 0 for the location found in top
      int right();
      ///  idxl
      // AC-GTACCGCA  ACGTACCGCA   more different, no change
      // ACCGTAC-GCA  ACCGTACGCA
      //       idx
      void shiftLeft(size_t topi);

      ///  idxl
      // AC-GTACCGCA  AC-GTACCGCA  more similar, two changes, if error greater
      // ACCGTAC-GCA  ACXGTACcGCA  X was used for postprocessing.
      //       idx
      void fixLeft(size_t topi);
      ///        idxr
      // ACCAGTAC-CA  ACCAGTACCA   more different
      // ACC-GTACGCA  ACCGTACGCA
      //    idx
      void shiftRight(size_t topi);

      ///        idxr
      // ACCAGTAC-CA  ACCAGTAC-CA   more different
      // ACC-GTACGCA  ACCaGTACXCA
      //    idx
      void fixRight(size_t topi);
      bool findAndFixLeft();
      bool findAndFixRight();
      void display();
      static void setWindowSize(int sz) { wdsz=sz; }
         
};
#endif

