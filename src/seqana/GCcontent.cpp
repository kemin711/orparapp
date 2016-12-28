#include <cstring>
#include <iostream>
#include <fstream>

#include <bioseq.h>

using namespace std;
using namespace orpara;


void usage() {
   cerr << "Usage: GCcontent <fasta file> default genomic.fas\n"
      << " or GCcontent -i fasfile\n"
      << "Options:\n"
      << "   -h or --help print help message\n"
      << "   -i inputfile this can also be given with the positional parameter\n";
   exit(1);
}
void writeResult(ostream &ous, long int A, long int C, long int G, long int T,
      long int N, long int genomiclen);

/** use fasta files to generate GC contesnt for all the sequences
 * in one file
 */
int main(int argc, char* argv[]) {
   //string fasfile="genomic.fas";
   string fasfile;
   string outfile;
   int i=0;
   while (i < argc) {
      if (!strcmp(argv[i], "-i"))
         fasfile = argv[++i];
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h"))
         usage();
      if (!strcmp(argv[i], "-o"))
         outfile = argv[++i];
      else {
         fasfile=argv[i];
      }
      ++i;
   }
   if (fasfile.empty()) usage();
   ifstream IN(fasfile);
   if (IN.fail()) {
      cerr << "Failed to open fasta file: " << fasfile << endl;
      usage();
      return 1;
   }
   if (outfile.empty()) {
      outfile = fasfile.substr(0, fasfile.rfind('.')) + ".GCcontent.res";
   }
   ofstream OU(outfile);
   DNA genomic;
   int sumgenomiclen=0;
   string header;
   long int A,C,G,T,N;
   A=C=G=T=N=0;
   while (genomic.read(IN, header)) {
      sumgenomiclen += genomic.length();
      OU << genomic.getName() << "\t" << genomic.GCContent(A,C,G,T,N) 
         << "\t" << genomic.length() << endl;
   }
   writeResult(cout, A,C,G,T,N, sumgenomiclen);
   writeResult(OU, A,C,G,T,N, sumgenomiclen);
   cerr << "result written to GCcontent.res\n";

   return 0;
}

void writeResult(ostream &ous, long int A, long int C, long int G, long int T,
      long int N, long int genomiclen) 
{
   long int total = A+C+G+T;
   if (total+N != genomiclen) {
      cerr << "total shuld be the same as sum of genomic len!\n";
   }
   ous << "Overal GC content: " << double(G+C)/total << " N content/(A+C+G+T+N) "
      << double(N)/(total+N) << endl
      << "Total genomic length: " << genomiclen << endl;
}

