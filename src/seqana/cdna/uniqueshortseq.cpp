#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <algorithm>

using namespace std;

void weedoutBadSeq(const string &infile, vector<char*> &allseq);
void readSeqFromFile(const string &infile, vector<char*> &seqs);
void eliminateSameAndSimple(vector<char*> &v);
void removeSimple(vector<char*> &seq);
bool isRepeat(const char *str);

int main(int argc, char* argv[]) {
   //string infile="all.fas";
   string infile, cleanfile, uniquefile;

   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "-i")) infile=argv[++i];
      else if (!strcmp(argv[i], "-c")) cleanfile=argv[++i];
      else if (!strcmp(argv[i], "-u")) uniquefile=argv[++i];
      else {
         cerr << "unrecognized option " << argv[i] << endl;
         exit(1);
      }
      ++i;
   }
   vector<char*> allgoodseq;
   if (!infile.empty()) weedoutBadSeq(infile, allgoodseq);
   else if (!cleanfile.empty()) {
      readSeqFromFile(cleanfile, allgoodseq);
   }
   else if (!uniquefile.empty()) {
      readSeqFromFile(uniquefile, allgoodseq);
      removeSimple(allgoodseq);
      return 0;
   }
   else {
      cerr << "must get sequence from -i rawIlluminaEST.fas or -c cleanIlluminaEST.fas\n";
      exit(1);
   }

   eliminateSameAndSimple(allgoodseq);
   return 0;
}

void readSeqFromFile(const string &infile, vector<char*> &seqs) {
   ifstream IN(infile.c_str());
   if (IN.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   char title[200], seq[500];
   IN.getline(title, 200);
   char *p, *pp;
   //vector<char*> allseq;
   while (!IN.eof() && title[0] == '>') {
      IN.getline(seq, 500);
      char *x = new char[strlen(seq)+1];
      strcpy(x, seq);
      seqs.push_back(x);
      IN.getline(title, 200);
   }
   IN.close();
   cerr << seqs.size() << " sequences read from "
      << infile << "\n";
}

void weedoutBadSeq(const string &infile, vector<char*> &allseq) {
   //string infile="all.fas";
   string oufile="allESTclean.fas";

   ifstream IF(infile.c_str());
   if (IF.fail()) {
      cerr << "Failed to open " << infile << endl;
      exit(1);
   }
   ofstream OU(oufile.c_str());
   // reading fasta formated file as stream
   char title[200], seq[500];
   IF.getline(title, 200);
   char *p, *pp;
   //vector<char*> allseq;
   while (!IF.eof() && title[0] == '>') {
      IF.getline(seq, 500);
      p=strchr(seq, 'N');
      bool good=false;

      if (p == NULL) good=true;
      else { // has at least one N
         if (*(p+1) == '\0') {
            *p='\0'; // N at end accept
            good=true;
         }
         else {
            ++p;
            pp=strchr(p, 'N');
            if (pp == NULL) good=true;
            else { // find location of second N
               if (pp-seq > 25) {
                  *pp='\0';
                  good=true;
               }
               else { // poor sequence
                  cout << title << endl;
                  cout << seq << endl;
               }
            }
         }
      }
      if (good) {
         OU << title << endl;
         OU << seq << endl;
         char *x = new char[strlen(seq)+1];
         strcpy(x, seq);
         allseq.push_back(x);
      }
      IF.getline(title, 200);
   }
   IF.close();
   OU.close();
}

bool compareCharPtr(const char *c1, const char *c2) {
   return strcmp(c1,c2) < 0;
}
bool samestr(const char *c1, const char *c2) {
   return strcmp(c1,c2) == 0;
}

void eliminateSameAndSimple(vector<char*> &v) {
   cerr << v.size() << " good sequences\n";
   sort(v.begin(), v.end(), compareCharPtr);
   cerr << "sorting done\n";
   vector<char*>::iterator e = unique(v.begin(), v.end(), samestr);
   vector<char*>::const_iterator it=v.begin();

   int cnt=1;
   int discard=0;
   ofstream OU("uniqueEST.fas");
   while (it != e) {
      if (isRepeat(*it)) {
         cerr << "Discarding simple reads\n";
         ++discard;
      }
      else {
         OU << ">IXkz" << cnt++ << endl << *it << endl;
      }
      delete[] *it;
      ++it;
   }
   cerr << cnt << " Unique sequences written to uniqueEST.fas\n"
      << discard << " simple sequences removed\n";
   while (it != v.end()) {
      delete[] *it;
      ++it;
   }
}

/** the seq shold have been cleared 
 * for testing, I am starting from the 
 * intermediate results. **/
void removeSimple(vector<char*> &seq) {
   int cnt=1;
   int discard=0;
   ofstream OU("uniqueESTNoSimple.fas");
   vector<char*>::iterator it=seq.begin();
   while (it != seq.end()) {
      if (isRepeat(*it)) {
         //cerr << "Discarding simple reads\n";
         ++discard;
      }
      else {
         OU << ">IXkz" << cnt++ << endl << *it << endl;
      }
      delete[] *it;
      ++it;
   }
   cerr << discard << " simple sequences removed\n";
}

/** there is not so many simple reads **/
bool isRepeat(const char *str) {
   const char* p=str;
   int t, a, c, g;
   int maxt=0, maxa=0, maxc=0, maxg=0;
   int length=0;
   while (*p != '\0') {
      t=0; a=0; g=0, c=0;
      while (*p != '\0' && (*p == 'T' ||  *p == 'N')) {
         ++t; ++p;
      }
      length += t;
      if (t > maxt) maxt=t;

      while (*p != '\0' && (*p == 'A' ||  *p == 'N')) {
         ++a; ++p;
      }
      length += a;
      if (a > maxa) maxa=a;

      while (*p != '\0' && (*p == 'C' ||  *p == 'N')) {
         ++c; ++p;
      }
      length += c;
      if (c > maxc) maxc=c;

      while (*p != '\0' && (*p == 'G' ||  *p == 'N')) {
         ++g; ++p;
      }
      length += g;
      if (g > maxg) maxg=g;
      if (*p == 'N') { ++length; ++p; continue; }
   }
   if (maxa/static_cast<double>(length) > 0.89) {
      //cerr << str << " is mostly A\n";
      return true;
   }
   else if (maxc/static_cast<double>(length) > 0.89) {
      //cerr << str << " is mostly C\n";
      return true;
   }
   else if (maxg/static_cast<double>(length) > 0.89) {
      //cerr << str << " is mostly G\n";
      return true;
   }
   else if (maxt/static_cast<double>(length) > 0.89) {
      //cerr << str << " is mostly T\n";
      return true;
   }
   else {
      return false;
   }
}

