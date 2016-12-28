#include "genoanalyzer.h"
#include "stagergap.h"

using namespace std;

void testGapInPoly() {
    // GATTTCAAA    GATTTCAA  TCATAAAAAACGACT TCATAAAAAAAACGACT
    // GATT-CAAA    GA-TTCAA  TCATAAAA--CGACT TCAT--AAAAAACGACT
    //
    // TCGTATTTTTTTCGGCATCA
    // TCGTACTTTT---GGCATCA
   vector<string> tops={"GATTTCAAA","GATTTCAA", "TCATAAAAAACGACT", "TCATAAAAAAAACGACT",
      "TCGTATTTTTTTCGGCATCA"};
   vector<string> bottoms={"GATT-CAAAA", "GA-TTCAA", "TCATAAAA--CGACT", "TCAT--AAAAAACGACT",
      "TCGTACTTTT---GGCATCA"};
   vector<size_t> idx = {4, 2, 8, 4};
   for (size_t i=0; i<tops.size(); ++i) {
      cout << "looking at gap at index " << idx[i] << tops[i][idx[i]] << endl;
      if (gapInPoly(tops[i], bottoms[i], idx[i])) {
         cout << "gap in tadem repeat:\n" << tops[i] << endl << bottoms[i] << endl;
         cout << endl;
      }
      else {
         cout << "gap not in tadem repeat:\n" << tops[i] << endl << bottoms[i] << endl;
         cout << endl;
      }
   }
}

void testCrossBottom() {
   //CCACAGGCCCCTGCACACCC-TCCCCGGCGCCAAGCTACTCCAGGGCGCTGTGGCGGGTG
   //|||||||||||||||||||| ||||||||||||| |||||| | |||||| |||||||||
   //CCACAGGCCCCTGCACACCCCTCCCCGGCGCCAAACTACTCTAAGGCGCTATGGCGGGTG
   vector<string> top={"ACCC-TCCCCGG"};
   vector<string> bottom={"ACCCCTCCCCGG"};
   vector<size_t> idx={4};
   for (size_t i=0; i<top.size(); ++i) {
      cout << top[i] << endl << bottom[i] << endl;
      crossBottom(top[i], bottom[i], idx[i]);
      cout << "after corssing\n";
      cout << top[i] << endl << bottom[i] << endl;
      cout << endl;
   }
}

void testStager() {
   // CATGCCAACGTGGGTACAAGGG-AGTCTGGCGGGGAGATGGCGTCATGCACACCACCTGCTCATGTGGAGCACA
   // |||||||||| ||||||| ||| ||||||||||| |||||| ||||||||||||| |||| |||| ||||||||
   // CATGCCAACGAGGGTACATGGGGAGTCTGGCGGG-AGATGGTGTCATGCACACCATCTGCCCATGCGGAGCACA
   //
   // AC-GTACCGCA  ACGTACCGCA   more different, no change
   // ACCGTAC-GCA  ACCGTACGCA
   //
   // AC-GTCACCGCA  AC-GTCACCGCA
   // ACCGTCAC-GCA  ACXGTCACcGCA 
   //
   // ACCAGTAC-CA  ACCAGTACCA 
   // ACC-GTACGCA  ACCGTACGCA
   //
   // ACCAGTACATC-CA  ACCAGTACATC-CA
   // ACC-GTACATCGCA  ACCaGTACATCXCA
   //
   vector<string> top={"CATGCCAACGTGGGTACAAGGG-AGTCTGGCGGGGAGATGGCGTCATGCACACCACCTGCTCATGTGGAGCACA",
      "AC-GTACCGCA", "AC-GTCACCGCA", "ACCAGTAC-CA", "ACCAGTACATC-CA"};
   vector<string> bottom={"CATGCCAACGAGGGTACATGGGGAGTCTGGCGGG-AGATGGTGTCATGCACACCATCTGCCCATGCGGAGCACA", 
      "ACCGTAC-GCA", "ACCGTCAC-GCA", "ACC-GTACGCA", "ACC-GTACATCGCA"};
   vector<size_t> idx={34, 7, 9, 3, 3};
   Stagergap::setWindowSize(15);
   int ti;
   for (size_t i=0; i<top.size(); ++i) {
      string t, b; // make a copy
      t=top[i];
      b=bottom[i];
      cout << "input\n";
      cout << top[i] << endl << bottom[i] << endl;
      Stagergap stg(top[i], bottom[i], idx[i]);
      if ((ti=stg.left()) != -1) {
         cout << "Stager left\n";
         stg.shiftLeft(ti);
         cout << "after shift\n";
         stg.display();
         Stagergap stg2(t,b,idx[i]);
         stg2.fixLeft(ti);
         cout << "after fix\n";
         stg2.display();
         cout << endl;
      }
      else if ((ti=stg.right()) != -1) {
         cout << "stager right\n";
         stg.shiftRight(ti);
         cout << "after shift\n";
         stg.display();
         Stagergap stg2(t,b,idx[i]);
         stg2.fixRight(ti);
         cout << "after fix\n";
         stg2.display();
         cout << endl;
      }
      else {
         cout << "no stager detected\n";
      }
   }
}

void readconfig() {
   string conf="comparetoref.conf";
   cerr << "rewading config file: " << conf << endl;
   ifstream ifs(conf.c_str());
   if (ifs.fail()) {
      cerr << "Failed to open " << conf << endl;
      exit(1);
   }
   map<string, string> confmap;
   string ln;
   getline(ifs, ln);
   while (!ifs.eof()) {
      if (ln[0] != '#' && !ln.empty()) {
         cerr << ln << endl;
         size_t s=ln.find('=');
         confmap[ln.substr(0,s)]=ln.substr(s+1);
      }
      getline(ifs, ln);
   }
   cerr << endl;
   for (auto ptr=confmap.begin(); ptr != confmap.end(); ++ptr) {
      cerr << ptr->first << " => " << ptr->second << endl;
   }
}

int main(int argc, char* argv[]) {
   testGapInPoly();
   testCrossBottom();
   testStager();
   readconfig();
   return 0;
}
