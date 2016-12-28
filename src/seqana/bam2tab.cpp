#include <iostream>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/SamHeader.h>
#include <api/SamSequenceDictionary.h>
#include <iterator>
#include <string>

using namespace std;


void showCigar(const vector<BamTools::CigarOp> &cigar) {
   for (int i=0; i<cigar.size(); ++i) {
      cout << cigar[i].Type << " " << cigar[i].Length << endl;
   }
}

int getAlignLength(const vector<BamTools::CigarOp> &cigar) {
   int alnlen=0;
   for (int i=0; i<cigar.size(); ++i) {
      alnlen += cigar[i].Length;
   }
   return alnlen;
}

int getDeletionLength(const vector<BamTools::CigarOp> &cigar) {
   int len=0;
   for (int i=0; i<cigar.size(); ++i) {
      if (cigar[i].Type == 'D') 
         len += cigar[i].Length;
   }
   return len;
}

int getInsertionLength(const vector<BamTools::CigarOp> &cigar) {
   int len=0;
   for (int i=0; i<cigar.size(); ++i) {
      if (cigar[i].Type == 'I') 
         len += cigar[i].Length;
   }
   return len;
}
void extractSummary(const string &bamfile, const string &tabfile);

void usage() {
   cerr << "Usage: bam2tab bamfilename.bam\n"
      << " You must specify a bam file as input\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   string bamfile, tabfile; //="bar1toref.bam";
   int i=1;
   while (i < argc) {
      if (strcmp(argv[i], "-b") == 0) {
         bamfile = string(argv[++i]);
      }
      else if (strcmp(argv[i], "-o") == 0) {
         tabfile = string(argv[++i]);
      }
      else if (argv[i][0] == '-') {
         cerr << "option " << argv[i] << " not recognized by this program\n";
         return 1;
      }
      else {
         bamfile = string(argv[i]);
      }
      ++i;
   }
   if (bamfile.empty()) {
      usage();
   }
   if (tabfile.empty()) 
      tabfile=bamfile.substr(0, bamfile.find('.')) + ".tab";
   extractSummary(bamfile, tabfile);

   return 0;
}

void extractSummary(const string &bamfile, const string &tabfile) {
   BamTools::BamAlignment bam;
   BamTools::BamReader breader;
   if (!breader.Open(bamfile)) {
      cerr << "Failed to open " << bamfile << endl;
      exit(1);
   }
   string tagstrval;
   int8_t tagi8val;
   int tagival;
   char tt;
   unsigned int taguival;

   const int MATCH=10;
   const int MMATCH=-9;
   const int GAPO=-29;
   const int GAPE=-3;

   //string tabfile=bamfile.substr(0, bamfile.find('.')) + ".tab";
   ofstream ouf(tabfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << tabfile << endl;
      exit(1);
   }
   BamTools::SamHeader sheader = breader.GetHeader();
   BamTools::SamSequenceDictionary sq=sheader.Sequences;
   // header info will not be repeated with all subsequent lines.
   ouf << "reference name: " << sq.Begin()->Name << "\t length: " << sq.Begin()->Length << endl;
   //ouf << "I have not gotten time to comput other numbers, it is not straight forward\n";

   ouf << "readname\treadlength\tscore\tidentical\talnlen\tnumgap\tgaplen\trefbegin\trefend\treadbegin\treadend\n";
   while (breader.GetNextAlignment(bam)) {
      //ouf << bam.Name << '\t' << bam.Length << '\t' << bam.Position;
      ouf << bam.Name << '\t' << bam.Length << '\t';
      if (!bam.IsMapped()) {
         ouf << "0\t0\n";
         continue;
      }
      vector<BamTools::CigarOp> cigar = bam.CigarData;
      int alnlen = getAlignLength(cigar);
      int readBegin, readend;
      vector<string> tags = bam.GetTagNames();
      if (cigar[0].Type == 'M') {
         readBegin=1;
      }
      else {
         cerr << "please guess read start position from " << cigar[0].Type << endl;
         readBegin=-1;
      }
      // note end pos is one short of begin pos + length
      int readEnd = readBegin + alnlen - getDeletionLength(cigar) - 1;
      int refEnd = bam.Position + alnlen - getInsertionLength(cigar);

      //vector<int> tagv(tags.size, 0);
      // get tag values
      //cout << "tags:\n";
      int score=0;
      int mismatch, gapo, gape, match;
      for (int i=0; i<tags.size(); ++i) {
         bam.GetTagType(tags[i], tt); // tt is tagtype single char
         //cout << "tag type: " << tt << " " << tags[i] << " => ";
         if (tt == 'Z') {
            bam.GetTag(tags[i], tagstrval);
            //cout << tagstrval << endl;
         }
         else if (tt == 'c') {
            bam.GetTag(tags[i], tagi8val);
            tagival=int(tagi8val); // must do a cast
            //cout << tagival << endl;
         }
         else if (tt == 'C') {
            bam.GetTag(tags[i], taguival);
            //cout << taguival << endl;
         }

         if (tags[i] == "XM") {
            mismatch=taguival;
         }
         else if (tags[i] == "XO") {
            gapo=taguival;
         }
         else if (tags[i] == "XG") {
            gape = taguival;
         }
      }
      match = alnlen - mismatch - gape;
      score = match*MATCH + mismatch*MMATCH + gapo*GAPO + (gape-gapo)*GAPE;
      /*
      cout << "mismatch: " << mismatch << " gap open: " << gapo 
         << " gap extend: " << gape 
         << " match: " << alnlen - mismatch - gape
         << endl;
      cout << "score: " << match*MATCH + mismatch*MMATCH + gapo*GAPO + (gape-gapo)*GAPE
         << endl;
      */
      ouf << score << '\t' << match << '\t' << alnlen << '\t' << gapo << '\t'
         << gape << '\t' << bam.Position+1 << '\t' << refEnd 
         << '\t' << readBegin << '\t' << readEnd << endl;
   }
   cout << "Alignment in tabular format written to " << tabfile << endl;
}
