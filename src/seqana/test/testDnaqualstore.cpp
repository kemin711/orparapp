#include "dnaqualstore.h"

using namespace std;

void testStore() {
   //string fastqFile="/remote/DataAnalysis/metag/pacbio/Aug20training/C01_1/Cccssid.fastq";
   string fastqFile="/remote/DataAnalysis/metag/pacbio/Duy0824/D1R1/rebdisclustrun/allsid.qcut.fastq";
   DNAQualstore store;
   store.setInputFile(fastqFile);
   store.readAverageQuality();
   store.save("testdnaqualstore.dqc");
   cout << store.getNumberUnique() << " unique sequences " << endl;
   cout << store.getTotalSequences() << " total sequences\n";
}

void testStoreOpen() {
   DNAQualstore store2;
   store2.open("testdnaqualstore.dqc");
   cout << store2.getNumberUnique() << " unique " 
      << store2.getTotalSequences() << " sequences " << endl;
}

int main(int argc, char* argv[]) {
   testStore();
   testStoreOpen();
   return 0;
}
