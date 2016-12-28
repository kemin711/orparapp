#include <thread>
#include <iterator>
#include <cstdlib>

#include <strformat.h>
#include "constaxmap.h"
#include "dnaqualstore.h"

using namespace orpara;

// the filesystem is not installed in gnu
// 5.2.0
//#include <experimental/filesystem>
//namespace fs = std::experimental::filesystem;


ConsensusParameter ConsensusParameter::onlyOne=ConsensusParameter();

void ConsensusParameter::probeFileSystem() {
   cerr << "probing file system settings\n";
   char* metagData = getenv("METAG_DATA");
   if (metagData) { dataPath = metagData; }
   else {
      cerr << "The directory METAG_DATA environment variable is missing!\n";
      //throw runtime_error("cannot find metag data directory");
   }
   setDefaultBindir();
   ifstream tmp(binaryPath);
   if (tmp.fail()) {
      cerr << "binary path " << binaryPath << " invalid\n";
      binaryPath="/usr/local/bin";
   }
   else {
      cerr << "using binary path: " << binaryPath << endl;
   }
}

void ConsensusParameter::setDefaultBindir() {
   char *metagHome = getenv("METAG_HOME");
   if (metagHome) {
      binaryPath = string(metagHome) + "/bin";
   }
   else {
      cerr << "METAG_HOME environment variable not set\n";
      //throw runtime_error("METAG_HOME environment variable not set");
   }
}

//////////////////////////////////////////////////////////////////
////// Blastrowd meaning default for the last d ////////////////
map<string, pair<string,string> > Blastnrowd::taxdic={};
// Note: the header for this class is different from the syster class
// Blastnrow where qlen, and slen are added
vector<string> Blastnrowd::header={"qseqid", "clustersz", "cldepth", "sseqid", "species", 
   "strain", "pident", "length", "mismatch", "gapopen", "qstart", "qend", 
   "sstart", "send", "bitscore"};

void Blastnrowd::writeHeader(ostream &ouf) {
   auto it = header.begin();
   ouf << *it;
   ++it;
   while (it != header.end()) {
      ouf << "\t" << *it;
      ++it;
   }
}

Blastnrowd::Blastnrowd(const vector<string> &row)
   : qseqid(row[0]), sseqid(row[1]),
     pident(stof(row[2])), length(stoi(row[3])),
     mismatch(stoi(row[4])), gapopen(stoi(row[5])),
     qstart(stoi(row[6])), qend(stoi(row[7])),
     sstart(stoi(row[8])), send(stoi(row[9])),
     bitscore(stoi(row[11])),
     taxon(), queryCount(0), consDepth(0)
{
   setTaxon();
}

void Blastnrowd::setTaxon() {
   map<string, pair<string,string> >::const_iterator mit=taxdic.find(sseqid);
   if (mit == taxdic.end()) {
      cerr << "failed to find " << sseqid << " in tax dictioanry\n";
      throw runtime_error(sseqid + " not in taxdictionary");
   }
   taxon=mit->second;
}

ostream& operator<<(ostream& ous, const Blastnrowd& r) {
   static char sep='\t';
   ous << r.qseqid << sep << r.queryCount << sep << r.getConsensusDepth() << sep
      << r.sseqid << sep << r.taxon.first << sep;
   if (r.taxon.second.empty()) ous << "\\N";
   else ous << r.taxon.second;
   ous << sep
      << r.pident << sep << r.length << sep
      << r.mismatch << sep << r.gapopen << sep
      << r.qstart << sep << r.qend << sep
      << r.sstart << sep << r.send << sep << r.bitscore;
   return ous;
}

bool Blastnrowd::operator<(const Blastnrowd& r) const {
   if (getQidNumber() < r.getQidNumber()) {
      return true;
   }
   if (getQidNumber() > r.getQidNumber()) {
      return false;
   }
   // same id
   if (getBitscore() > r.getBitscore()) {
      return true;
   }
   return false;
}

tuple<string, int, string> Blastnrowd::getSpeciesDepth() const { 
   if (taxon.second.empty()) 
      return make_tuple(getQuery(), getConsensusDepth(), taxon.first);
   return make_tuple(getQuery(), getConsensusDepth(), taxon.first + " " + taxon.second); 
}

vector<string> Blastnrowd::toStringVector() {
   vector<string> res;
   res.push_back(qseqid); res.push_back(to_string(queryCount)); res.push_back(to_string(consDepth));
   res.push_back(sseqid); res.push_back(taxon.first); res.push_back(taxon.second);
   res.push_back(ftos(pident,3));
   //res.push_back(to_string(qlen)); // only exist in the sister Blastnrow class
   //res.push_back(to_string(slen));
   res.push_back(to_string(length));
   res.push_back(to_string(qstart));
   res.push_back(to_string(qend));
   res.push_back(to_string(sstart));
   res.push_back(to_string(send));
   res.push_back(to_string(bitscore));
   return res;
}  


//////////////////////////////////////////////////
//// some helper functions /////////

// simply id tax
map<string,string> slurpSilvaTaxonomy() {
ConsensusParameter& param=ConsensusParameter::getInstance();
//string filePath="/remote/DataAnalysis/metag/refseq/silva123.header";
string filePath=param.getTaxFilePath();
   cerr << "reading silva taxonomy " << filePath << " into memory ...\n";
   map<string,string> id2tax;
   ifstream inf(filePath);
   if (inf.fail()) {
      throw runtime_error(filePath + " invalid or no right to read");
   }
   string line;
   getline(inf, line);
   while (!inf.eof()) {
      string id=line.substr(0, line.find(' '));
      string tax=line.substr(line.find(' ')+1);
      id2tax[id]=tax;
      getline(inf, line);
   }
   cerr << id2tax.size() << " entries read\n";
   return id2tax;
}

// id => {species, strain}
map<string, pair<string, string> > loadTaxdic() {
   ConsensusParameter& param=ConsensusParameter::getInstance();
   //string filePath="/remote/DataAnalysis/metag/refseq/silva123.header";
   string filePath=param.getTaxFilePath();
   cerr << "loading taxonomy " << filePath << " into memory ...\n";
   map<string, pair<string, string> > id2tax;
   ifstream inf(filePath);
   if (inf.fail()) {
      throw runtime_error(filePath + " invalid or no right to read");
   }
   string line;
   getline(inf, line);
   while (!inf.eof()) {
      vector<string> row = split(line, '\t');
      id2tax[row[0]]=make_pair(trimSpace(row[1]), trimSpace(row[2]));
      getline(inf, line);
   }
   cerr << id2tax.size() << " entries read\n";
   return id2tax;
}

string runBlast(const string &fastaFile, bool lazy) {
   cerr << "running tabular blast on " << fastaFile << " ...\n";
   int percentIdentityCut=95;
   int numhits=200;
   ConsensusParameter& param=ConsensusParameter::getInstance();
   string blastOutFile = fastaFile.substr(0, fastaFile.rfind('.')) 
      + "_" + param.getRefdb() + ".bln.tab";
   // should use the maximum number threads available in the hardware
   //vector<string> blastParameter={"-db", "silva123.fas", 
   // now the blast db is controlled by the consensus parameter.
   vector<string> blastParameter={"-db", param.getBlastdb(), 
      "-num_threads", to_string(thread::hardware_concurrency()), 
      "-out", blastOutFile, "-evalue", "0.001", "-word_size", "15", 
      "-perc_identity", to_string(percentIdentityCut), 
      "-query", fastaFile, "-max_target_seqs", to_string(numhits), 
      "-outfmt", "'6 qseqid sseqid pident qlen slen length qstart qend sstart send  bitscore'"};
   if (lazy) {
      ifstream tmp(blastOutFile);
      if (tmp.good()) { // now check the headers match or not
         string line;
         getline(tmp, line);
         vector<string> row = split(line, '\t');
         if (row.size() == 11) {
            cerr << "blast done before, lazy version not redoing it\n";
            return blastOutFile;
         }
      }
   }
   //string binarydir="/remote/RSU/sw-cache/metag/bin";
   //string blastnpath=binarydir + "/blastn";
   string blastnpath=param.getBlastnPath();
   string blastcmd = blastnpath + " "
      + accumulate(blastParameter.begin(), blastParameter.end(), string(),
         [](const string& a, const string& b)->string { return a+(a.length()>0? " " : "") + b; });
   cerr << "blast command: " << blastcmd << endl;
   int rv=system(blastcmd.c_str());
   cerr << "blast done with return value: " << rv << endl;
   if (rv > 0) {
      cerr << "Failed to run blast\n";
      exit(1);
   }
   cerr << "tabular blast done. " << blastOutFile << "\n";
   return blastOutFile;
}

