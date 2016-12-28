#include "consensus.h"
#include "dnaqualstore.h"
#include <strformat.h>
#include <iterator>
#include <cstdio>
#include <cstdlib>

using namespace std;

string binarydir="/remote/RSU/sw-cache/metag/bin";

string runswarm(const string &idcntFile, int nt);
void writeConsensusFasta(const vector<pair<string,int> > &con, const string &fasFile);
void alignSingleton(const vector<pair<string,int> > &cons, const string &alnfile);
/**
 * @return the blast result file name.
 */
string runBlast(const string &fastaFile);
map<string,string> slurpSilvaTaxonomy();
void writeTaxMapping(const string &blastFile, const vector<pair<string,int> > &cons, 
      const map<string,string> &tax, const string &mappingFile);

int main(int argc, char* argv[]) {
   string fastqfile;
   int numthread=8;
   //string fastqfile="/home/kzhou/work/metag/A_ccs.shortid.fastq";
   int i=1;
   while (i<argc) {
      if (!strcmp(argv[i], "-i")) fastqfile=argv[++i];
      else if (!strcmp(argv[i], "-t")) numthread=atoi(argv[++i]);
      else {
         fastqfile=argv[i];
      }
      ++i;
   }
   DNAQualstore store(fastqfile);
   // produce otuput for swarm input
   string fstem = fastqfile.substr(0, fastqfile.rfind('.'));
   string fastaCntFile = fstem + "cnt.fasta";
   store.writeFasta(fastaCntFile);
   string swarmResultFile=runswarm(fastaCntFile, numthread);
   Consensus::setStore(&store);
   vector<vector<string> > cluster=readCluster(swarmResultFile);
   // The number of input sequence is the int
   vector<pair<string,int> > consensusResult;
   cout << "testing Consensus object\n";
   string conFile, baseFile, qualFile;
   conFile = fstem + ".consensus.txt";
   baseFile = fstem + ".base.txt";
   qualFile = fstem + ".qual.txt";
   cerr << "doing consensus building ...\n";
   int every=2;
   for (int i=0; i<cluster.size(); ++i) {
      Consensus cons(cluster[i]);
      cons();
      consensusResult.push_back(make_pair(cons.getConsensus(), cons.numseq()));
      cons.printResult(conFile, baseFile, qualFile);
      if (cons.numseq() > 1) {
         cerr << "working on cluster " << i << " ...\n";
      }
      else if ((i+1) % every == 0) {
         cerr << "working on cluster " << i << " ...\n";
         every = ceil(pow(every, 1.5));
      }
   }
   string conFasFile = fstem + "consensus.fas";
   writeConsensusFasta(consensusResult, conFasFile);
   string singlealnFile("singleton.aln");
   alignSingleton(consensusResult, singlealnFile);
   string blresult=runBlast(conFasFile);
   map<string, string> taxdic=slurpSilvaTaxonomy();
   string taxmapFile="swarm_out_taxmap.tab";
   writeTaxMapping(blresult, consensusResult, taxdic, taxmapFile);
   return 0;
}

/**
 * @param blastFile input of the tabular blast
 * @param cons consensus and its size. Only need size.
 * @param tax dictionary of blast subj to silva taxonomy
 * @param mappingFile is the output file.
 */
void writeTaxMapping(const string &blastFile, const vector<pair<string,int> > &cons, 
      const map<string,string> &tax, const string &mappingFile) 
{
   cerr << "writing otu taxonomy mapping to file ...\n";
   // the header has been updated in newer version 
   vector<string> header={"clusterid", "clsize", "hitid", "organism_name", "percent_ident", "consensesu_len", "hit_len", "align_len", "qstart", "qend", "sstart", "send", "bitscore"};
   ifstream inf(blastFile);
   ofstream ouf(mappingFile);
   copy(header.begin(), header.end(), ostream_iterator<string>(ouf, "\t"));
   ouf << endl;
   string line;
   getline(inf, line);
   vector<string>::iterator it;
   //map<string,string>::const_iterator mit;
   while (!inf.eof()) {
      vector<string> row = split(line, '\t');
      int clusterId=stoi(row[0].substr(2));
      int clsize=cons[clusterId].second;
      const string taxName=tax.find(row[1])->second;
      it=row.begin();
      it += 3;
      row.insert(it, taxName);
      it=row.begin();
      ++it;
      row.insert(it, to_string(clsize));
      copy(row.begin(), row.end(), ostream_iterator<string>(ouf, "\t"));
      ouf << endl;
      getline(inf, line);
   }
   cerr << "swarm tax otu mapping written to " << mappingFile << endl;
}


map<string,string> slurpSilvaTaxonomy() {
   cerr << "reading silva taxonomy into memory ...\n";
   string filePath="/refseq/silva123.header";
   map<string,string> id2tax;
   ifstream inf(filePath);
   string line;
   getline(inf, line);
   while (!inf.eof()) {
      string id=line.substr(0, line.find(' '));
      string tax=line.substr(line.find(' ')+1);
      id2tax[id]=tax;
      getline(inf, line);
   }
   return id2tax;
}

/**
For tabular output:
    qseqid sseqid pident qlen slen length qstart qend sstart send  bitscore
*/
string runBlast(const string &fastaFile) {
   cerr << "starting to run tabular blast ...\n";
   int percentIdentityCut=91;
   string blastOutFile = fastaFile.substr(0, fastaFile.rfind('.')) + "swarmbln.tab";
   vector<string> blastParameter={"-db", "silva123.fas", "-num_threads", "8", 
      "-out", blastOutFile, "-evalue", "0.001", "-word_size", "14", 
      "-perc_identity", to_string(percentIdentityCut), 
      "-query", fastaFile, "-max_target_seqs", "1", 
      "-outfmt", "'6 qseqid sseqid pident qlen slen length qstart qend sstart send  bitscore'"};
   string binarydir="/remote/RSU/sw-cache/metag/bin";
   string blastnpath=binarydir + "/blastn";
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
   cerr << "tabular blast done\n";
   return blastOutFile;
}

string runswarm(const string &idcntFile, int nt) {
   string outfile=idcntFile.substr(0, idcntFile.rfind('.')) + ".swarm.out";
   string swarmCmd=binarydir + "/swarm -d 10 -t " + to_string(nt) + " -o " + outfile + " " + idcntFile;
   cerr << "running swarm ...\n" << swarmCmd << endl;
   int rv=system(swarmCmd.c_str());
   if (rv != 0) {
      cerr << "failed to run command: " << swarmCmd << endl;
   }
   return outfile;
}
void alignSingleton(const vector<pair<string,int> > &cons, const string &alnfile) {
   cerr << "align singleton to large clusters ..\n";
   ofstream ouf(alnfile);
   int s=0;
   while (s < cons.size() && cons[s].second > 1) ++s;
   cerr << "singleton index at " << s << endl;
   // s points to the first singleton
   DNA seed("cluster0", cons[0].first, "clusterSize=" + to_string(cons[0].second));
   SimpleScoreMethod sm(11, -10, -37, -3);
   Dynaln<SimpleScoreMethod> aligner(sm), aligner2(sm);
   aligner.setSeq1(seed);
   aligner2.setSeq1(seed);
   int every=2;
   int good=0;
   for (int i=s; i<cons.size(); ++i) {
      if ((i-s+1) % every == 0) {
         cerr << "aligning " << i << " " << aligner.getIdentity() << endl;
         every = ceil(pow(every, 1.5));
      }
      DNA single("cl"+to_string(i), cons[i].first);
      aligner.setSeq2(single);
      aligner.runlocal();
      DNA rcdna;
      bool aln2better=false;
      double identity=aligner.getIdentity();
      if (identity < 0.80 && aligner.getCov1() < 0.8 && aligner.getCov2() < 0.8) {
         ouf << "forward align not good enough, trying reverse complement\n";
         rcdna=single.revcompCopy();
         aligner2.setSeq2(rcdna);
         aligner2.runlocal();
         if (aligner2.getScore()>aligner.getScore()) {
            aln2better=true;
            identity=aligner2.getIdentity();
         }
      }
      if (identity > 0.87) ++good;
      if (identity > 0.70) {
         if (aln2better) 
            aligner2.printAlign(ouf,90);
         else
            aligner.printAlign(ouf, 90);
      }
   }
   cerr << "alignments written to " << alnfile << endl;
   cerr << good << " matched well to large cluster\n";
   cerr << double(good)/(cons.size()-s) << " fraction good\n";
}

void writeConsensusFasta(const vector<pair<string,int> > &con, const string &fasFile) {
   ofstream ouf(fasFile);
   int i=0;
   while (i<con.size()) {
      ouf << DNA("cl" + to_string(i), con[i].first, 
            "cluster size=" + to_string(con[i].second));
      i++;
   }
}

