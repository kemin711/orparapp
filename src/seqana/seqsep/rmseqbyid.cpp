#include<iostream>
#include<string>
#include<fstream>
#include<set>
#include<cstring>
#include<cstdlib>

#include<strformat.h>

using namespace std;
using namespace orpara;

set<string> readIds(const string &idfile);
void splitFasta(const string &infile, const string &idfile,
      const string &removed, const string &saved);
/** shorter version of the above */
void splitFasta(const string &infile, const string &idfile);

int main(int argc, char* argv[]) {
   int i=1;
   string infile="Accs2.filter.good.align";
   string removeIdFile="Accs2.filter.good.ref.uchime.accnos";

   while (i<argc) {
      if (!strcmp(argv[i], "--id")) removeIdFile=argv[++i];
      else {
         infile=argv[i];
      }
      ++i;
   }
   splitFasta(infile, removeIdFile);

   return 0;
}

set<string> readIds(const string &idfile) {
   set<string> ids;
   cerr << "reading ids from " << idfile << endl;
   ifstream inf(idfile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open idlist file: " << idfile << endl;
      exit(1);
   }
   string line;
   getline(inf, line);
   size_t i;
   while (!inf.eof()) {
      if ((i=line.find_first_of("\t ")) != string::npos) {
         set<string> row = digest2set(line, " \t");
         ids.insert(row.begin(), row.end());
      }
      else ids.insert(line);
      getline(inf, line);
   }
   inf.close();
   cerr << ids.size() << " ids collected\n";
   return ids;
}

void splitFasta(const string &infile, const string &idfile) {
   string rmfile = infile.substr(0, infile.rfind('.')) + ".removed" + infile.substr(infile.rfind('.'));
   string svfile = infile.substr(0, infile.rfind('.')) + ".saved" + infile.substr(infile.rfind('.'));
   splitFasta(infile, idfile, rmfile, svfile);
}

void splitFasta(const string &infile, const string &idfile,
      const string &removed, const string &saved) {
   cerr << "split fasta file: " << infile  << " according to " << idfile << endl;
   set<string> rmids = readIds(idfile);
   set<string>::const_iterator sit;
   ofstream outrm(removed.c_str());
   ofstream outsv(saved.c_str());
   if (outrm.fail()) {
      cerr << "Failed to open output file " << removed << endl;
      exit(1);
   }
   if (outsv.fail()) {
      cerr << "Failed to open output file " << saved << endl;
      exit(1);
   }

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open input fasta file\n";
      exit(1);
   }
   string line;
   getline(inf, line);
   size_t i;
   string id;
   size_t cnt=0;
   while (!inf.eof() && line[0] == '>') {
      if ((i=line.find(' ')) != string::npos) {
         id = line.substr(1, i);
      }
      else id = line.substr(1);
      if (rmids.find(id) != rmids.end()) {
         outrm << line << endl;
         getline(inf, line);
         while (!inf.eof() && line[0] != '>') {
            outrm << line << endl;
            getline(inf, line);
         }
      }
      else {
         ++cnt;
         outsv << line << endl;
         getline(inf, line);
         while (!inf.eof() && line[0] != '>') {
            outsv << line << endl;
            getline(inf, line);
         }
      }
   }
   cerr << "file " << infile << " split into " << saved << " and " << removed << endl;
   cerr << rmids.size() << " sequence removed " << cnt << " saved\n";
}





