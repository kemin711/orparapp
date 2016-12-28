#include <fstream>
#include <iostream>
#include <cstdlib>
#include <string>
#include <map>

#include <fastq.h>
#include <strformat.h>

using namespace std;
using namespace orpara;

void usage();
/**
 * @return the oputput file name with a pattern infile.stem_trimWC.fastq.
 *         where W and C are two integers of window and cutoff values.
 */
string nameOutputFile(const string &infile);
string nameOutputFileWithExtention(const string &infile, const string &ext);
string getFileStem(const string &infile);
int writeUnique(istream &ins, ostream &ous, map<string, int> &names);
int readOnly(istream &ins, map<string, int> &names);

/**
 * Check the duplicated fastq sequences by their id.
 */
int main(int argc, char* argv[]) {
   string infile, outfile, uniqueFile;
   int i=1;
   while (i < argc) {
      if (string(argv[i]) == "-i") infile = argv[++i];
      else if (string(argv[i]) == "-o") outfile = argv[++i];
      else if (string(argv[i]) == "-u") uniqueFile = argv[++i];
      else if (string(argv[i]) == "--help" || argv[i][0]=='?') {
         usage(); return 1;
      }
      else {
         infile=argv[i];
         if (i+1 < argc && argv[i+1][0] != '-') {
            ++i;
            outfile = argv[i];
         }
      }
      ++i;
   }

   if (infile.empty()) {
      usage();
      return 1;
   }
   if (outfile.empty()) {
      outfile = nameOutputFile(infile);
   }

   ifstream inf(infile.c_str());
   if (inf.fail()) {
      cerr << "Failed to open " << infile << endl;
      return 1;
   }
   int count;
   map<string, int> name;
   if (!uniqueFile.empty()) {
      ofstream unq(uniqueFile.c_str());
      if (unq.fail()) {
         cerr << "Failed to open " << uniqueFile << endl;
         return 1;
      }
      count = writeUnique(inf, unq, name);
      cout << "unique fastq sequence writtent to " << uniqueFile << endl;
   }
   else count = readOnly(inf, name);
   inf.close();

   ofstream ouf(outfile.c_str());
   if (ouf.fail()) {
      cerr << "Failed to open " << outfile << " for writing " << endl;
      return 1;
   }
   map<string, int>::const_iterator it;
   i=0;
   for (it = name.begin(); it != name.end(); ++it) {
      if (it->second > 1) {
         ouf << it->first << "\t" << it->second << endl;
         cout << it->first << "\t" << it->second << endl;
         ++i;
      }
   }
   cout << i << " sequence names duplicated out of " << count << endl;
   ouf.close();

   return 0;
}

int readOnly(istream &ins, map<string, int> &names) {
   Fastq fasq;
   int count=0;
   while (fasq.read(ins)) {
      ++count;
      ++names[fasq.getName()];
   }
   return count;
}

int writeUnique(istream &ins, ostream &ous, map<string, int> &names) {
   Fastq fasq;
   int count=0;
   string name;
   while (fasq.read(ins)) {
      ++count;
      name = fasq.getName();
      ++names[name];
      if (names[name] == 1) {
         fasq.write(ous);
      }
   }
   return count;
}

void usage() {
   cout << "Usage: checkdupfasq <fastq filename> \n"
      << "  --help will print a help messasge\n"
      << "  checkdupfasq ? will also print the help message\n"
      << "  You can also do this: checkdupfasq infile outfile\n"
      << "  Options:\n"
      << "      -i input fastq file\n"
      << "      -o output duplicate fastq name file. If not given will be constructed automatically\n"
      << " All options are optional. The input file \n" 
      << " Can also be given directly without the -i options\n";
}

/**
 * @return the trim output file name for a given input file name
 *    and parameters used.
 */
string nameOutputFile(const string &infile) {
   return getFileStem(infile) + "_duplicated.list";
}

string nameOutputFileWithExtention(const string &infile, const string &ext) {
   return getFileStem(infile) + ext;
}

string getFileStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      stem = infile.substr(0, i);
   }
   return stem;
}
