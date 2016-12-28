#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <strformat.h>

using namespace std;
using namespace orpara;


string nameOutfile(const string &infile) {
   string outfile;
   size_t i=infile.rfind('.');
   if (i != string::npos) {
      outfile = infile.substr(0, i+1) + "cnt.tab";
   }
   else {
      outfile = infile + "_cnt.tab";
   }
   return outfile;
}

void writeTable(const map<pair<int, int>, int>& pcnt, const string &outfile) {
   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writing\n";
      exit(1);
   }
   for (auto itr=pcnt.begin(); itr != pcnt.end(); ++itr) {
      ouf << itr->first.first << '\t' << itr->first.second << '\t' << itr->second << endl;
   }
   cerr << "count written to " << outfile << endl;
}

pair<int,int> parseIntPair(const char* twoint) {
   // -c 1,2 format
   const char *pch = strchr(twoint, ',');
   char tmp[50];
   strncpy(tmp, twoint, pch-twoint);
   int first = atoi(tmp);
   ++pch;
   strcpy(tmp, pch);
   int second = atoi(tmp);
   return make_pair(first, second);
}

int main(int argc, char* argv[]) {
   pair<int, int> dataColumn(1, 2);
   string infile;
   char sep='\t';

   // command line argument tabulatepair 1 2
   int i=1;
   while (i < argc) {
      if (strcmp(argv[i], "-i") == 0) {
         infile=string(argv[++i]);
      }
      else if (strcmp(argv[i], "-1") == 0) {
         dataColumn.first=atoi(argv[++i]);
      }
      else if (strcmp(argv[i], "-2") == 0) {
         dataColumn.second=atoi(argv[++i]);
      }
      else if (strcmp(argv[i], "-c") == 0) { // -c 1,2
         dataColumn=parseIntPair(argv[++i]);
      }
      else {
         infile=string(argv[i]);
      }
      ++i;
   }

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open file " << infile << endl;
      return 1;
   }
   string line;
   getline(inf, line); // read header
   getline(inf, line); // read first line of data
   map<pair<int, int>, int> paircnt;
   while (!inf.eof()) {
      vector<string> row = split(line, sep);
      ++paircnt[make_pair(stoi(row[dataColumn.first]), stoi(row[dataColumn.second]))];
      getline(inf, line);
   }
   writeTable(paircnt, nameOutfile(infile));

   return 0;
}
