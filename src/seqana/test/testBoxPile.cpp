#include "boxpile.h"
#include <glob.h>

using namespace std;

/**
 * This is only useful for command line version
 * Test the BoxPile class
 * it produces three bar graphs.
 */
string findDefaultInput();

void usage(const string &arg0) {
   string program=arg0;
   if (program.rfind('/') != string::npos) {
      program = program.substr(program.rfind('/'));
   }
   cerr << program << " -o otubar.svg otu_consensus_silvaplus.goodmerged.tab\n"
      << "Options\n"
      << "   -o is the output file name. The file is an SVG file\n"
      << "Positional arguments is the input file.\n"
      << "If you don't specify any argument this program will find one\n"
      << "You must run this program in the directory where you have\n"
      << "run the 16S pipeline or you can specify a full path to \n"
      << "the input file\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   //string infile="german_D.tab";
   string infile;
   string outfile="outbar.svg";

   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-o")) {
         outfile = argv[++i];
      }
      else if (argv[i][0] == '-') {
         cerr << "wrong option " << argv[i] << endl;
         return 1;
      }
      else if (!strcmp(argv[i], "-h") ||
            !strcmp(argv[i], "--help")) {
         usage(argv[0]);
      }
      else {
         infile = argv[i];
      }
      ++i;
   }
   if (infile.empty()) {
      infile=findDefaultInput();
   }

   BoxPile bargraph;
   bargraph.readOtumap(infile);
   bargraph.draw(outfile);
   bargraph.draw2("otubar2.svg");
   bargraph.draw3("otubar3.svg");
   
   return 0;
}

string findDefaultInput() {
   glob_t glob_res;
   glob("otu_consensus_*.goodmerged.tab", GLOB_TILDE, NULL, &glob_res);
   cerr << glob_res.gl_pathc << " results\n";
   for (size_t i=0; i<glob_res.gl_pathc; ++i) {
      cout << glob_res.gl_pathv[i] << endl;
   }
   if (glob_res.gl_pathc > 1) {
      cerr << "you must specify one file as input\n";
      exit(1);
   }
   else {
      return string(glob_res.gl_pathv[0]);
   }
}
      
