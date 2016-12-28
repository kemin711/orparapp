#include <iostream>
#include <fstream>

#include <strformat.h>
#include <bioseq.h>

using namespace std;
using namespace orpara;

/**
 * the sequence identifier used by Silva
 * contain subsequence information.
 * For a small subset, multiple copies of 16S are
 * found from the same genome annotation. They share
 * the same parent sequence id. So you cannot
 * remove them.
 */
void extractHeader(const string &fasFile, const string &newfile, const string &header, int lencut);
string getStem(const string &infile) {
   string stem = infile;
   string::size_type i = infile.rfind('.');
   if ( i != string::npos) {
      stem = infile.substr(0, i);
   }
   return stem;
}

void usage() {
   cerr << "parsesilvassu -h header -o shortid.fas inputfile\n"
      << "parses small subunit including 16S and 18S ribisomal rRNA file\n";
   exit(1);
}

int main(int argc, char* argv[]) {
   int i=1;
   string inputfile, outputfile, headerfile;
   //string inputfile="Bacteria16S";
   //string outputfile="Bacteria16Ssid.fas";
   //string headerfile="Bacteria16Sheader.tab";
   int lengthcut=700;
   while (i < argc) {
      if (!strcmp(argv[i], "-t")) {
         headerfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-o")) {
         outputfile = argv[++i];
      }
      else if (!strcmp(argv[i], "-l")) {
         lengthcut = atoi(argv[++i]);
      }
      else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
         usage();
      }
      else {
         inputfile = argv[i];
      }
      ++i;
   }
   if (inputfile.empty()) usage();
   if (outputfile.empty()) {
      outputfile = getStem(inputfile) + "cut.fas";
   }
   if (headerfile.empty()) {
      headerfile = getStem(inputfile) + "header.tab";
   }
   extractHeader(inputfile, outputfile, headerfile, lengthcut);
}


string cleanName(const string &input) {
   string name=input;
   if (input[0] == '[') {
      name=input.substr(1);
      name.erase(name.rfind(']'), 1);
   }
   if (input[0] == '\'') {
      name=input.substr(1);
      name.erase(name.rfind('\''));
   }
   return name;
}

bool endWithQuotedSpecies(const string &sp) {
   string::size_type i=sp.find('(');
   string::size_type j=sp.rfind(')');
   if (j == sp.length()-1) {
      if (isScientificSpeciesName(sp.substr(i+1,j-i-1))) {
         return true;
      }
   }
   return false;
}

bool isSymbion(const string &species, const string &strain) {
   if (species.find(" sp. ") != string::npos) {
      string::size_type i=species.find('(');
      string::size_type j=species.rfind(')');
      if (isScientificSpeciesName(species.substr(i+1,j-i-1)) 
            && (isupper(strain) || wc(strain) == 1)) {
         return true;
      }
   }
   else if (endWithQuotedSpecies(species)) {
      return true;
   }
   return false;
}

string extractStrain(string &species) {
   string strain;
   if (wc(species) > 2) {
      string::size_type i = species.rfind(" strain ");
      if (i != string::npos) {
         strain = species.substr(i+8);
         species = species.substr(0, i);
      }
      else if ((i=species.rfind(" str. ")) != string::npos) {
         strain = species.substr(i+6);
         species = species.substr(0,i);
      }
      else {
         vector<string> word = split(species, ' ');
         string newsp = word[0] + ' ' + word[1];
         strain = species.substr(newsp.length());
         species = newsp;
      }
   }
   return strain;
}

bool isAllowedName(const string &name) {
   static vector<string> allowedSpecies = { "haloarchaeon", "cyanobiont", "Candidatus",
      "marine", "actinomycete", "archaeon" };
   string tmp = firstword(name);
   for (unsigned i=0; i<allowedSpecies.size(); ++i) {
      if (tmp == allowedSpecies[i]) {
         return true;
      }
   }
   vector<string> word = split(name, ' ');
   if (word[1] == "bacterium" || word[1] == "proteobacterium") return true;
   return false;
}

bool isGenericBacterium(const string &sp) {
   static vector<string> unnamedBacteria = {"coliform", "soil"};
   if (sp.find("bacterium") != string::npos) {
      for (unsigned i = 0; i< unnamedBacteria.size(); ++i) {
         string tmp = unnamedBacteria[i] + " bacterium";
         if (startwith(sp, tmp)) return true;
      }
   }
   return false;
}

/**
 * return species, strain pair
 */
pair<string, string> breakupHeader(const string &hd) {
   vector<string> path=split(hd, ';');
   string species = path[path.size() - 1];
   species = cleanName(species);
   path.resize(path.size()-1);
   string lowerestTaxon = path[path.size()-1];
   while (lowerestTaxon == "uncultured") {
      path.resize(path.size()-1);
      lowerestTaxon = path[path.size()-1];
   }
   string strain = extractStrain(species);
   if (firstword(species) != lowerestTaxon) { // good annotation
      if (firstword(species) == "uncultured") {
         species = lowerestTaxon + " uncultured";
      }
      else if (species == "unidentified") {
         species = lowerestTaxon + "unidenfified";
      }
      else if (isGenericBacterium(species) || endwith(firstword(species), "bacterium")) {
         species = lowerestTaxon + " " + species;
      }
      else if (!isScientificSpeciesName(species) && !endwith(species, " sp.")) {
         if (isAllowedName(species) || isupper(species[0])
               || species.find("bacterium") != string::npos
               || species.find("endosymbiont") != string::npos) {
         }
         else {
            //cout << endl << hd << endl;
            //cout << species << " is not scientific\n";
         }
      }
   }
   return make_pair(species, strain);
}

void extractHeader(const string &fasFile, const string &newfile, const string &header, int lencut) {
   DNA seq;
   ifstream inf(fasFile);
   if (inf.fail()) {
      throw runtime_error("failed to open file: " + fasFile);
   }
   ofstream ouf(newfile);
   ofstream oufh(header);
   int cnt=0;
   while (seq.read(inf)) {
      //cout << seq.getName() << endl
      //   << seq.getTitle() << endl;
      if (seq.getTitle().find("synthetic") != string::npos
            || seq.length() < unsigned(lencut)) 
         continue;
      //cout << seq.getName() << endl;
      pair<string, string> taxon = breakupHeader(seq.getTitle());
      oufh << seq.getName() << '\t' << taxon.first << '\t' << taxon.second << endl;
      ouf << seq;
      ++cnt;
   }
   cout << cnt << " fast sequences with shorter id written to " << newfile << endl;
}

