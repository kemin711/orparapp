#include <iostream>
#include <thread>
#include <dynalnt.h>
#include <bioseq.h>

using namespace std;

class TestA {
   private:
      int id;

   public:
      TestA() : id(1) { }
      int nextId() { return id++; }
};

void foo() {
   srand(time(0));
   TestA ta;
   for (int i=0; i<5; ++i) {
      cout << i << " " << rand() << endl;
      cout << "next id: " << ta.nextId() << endl;
   }
}

void alignseq(const DNA &d1, const DNA &d2) {
   SimpleScoreMethod sm(10, -9, -30, -1);
   Dynaln<SimpleScoreMethod> aligner(d1, d2, sm);
   aligner.runlocal();
   aligner.printAlign(cout);
}


int main(int argc, char* argv[]) {
   string seq1="TGGATATACACGGTGTTGACTGATTTCAAGACCTGGCTCCAGTACAAGCTCCTTCCGCGCTTGCCGGGAG";
   string seq2="TGGATATGCACGGTTTTGGCTGATCTCAAGACCTGCCTCCAGTCCAAGCTCCTGCCGCGATTGCCGGGAG";

   DNA d1("d1", seq1);
   DNA d2("d2", seq2);

   thread first(foo);
   thread second(foo);
   thread thrid(alignseq, d1, d2);
   first.join();
   second.join();
   thrid.join();
   cout << "threads completed\n";
   return 0;
}
