#ifndef CONSTAXMAP_H
#define CONSTAXMAP_H

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <iterator>
#include <strformat.h>
#include <list>
#include <tuple>

#include <bioseq.h>
#include "boxpile.h"
#include "piechart.h"
#include "textpile.h"

using namespace std;
using namespace orpara;

/**
 * Parameter used by runblast and slurptaxonomy
 * This should be used by main and passed as
 * parameter to blast and sluprtaxonomy.
 * This is a singleton class so there is no need
 * to pass around this parameter.
 * Just get the singleton object and use it.
 */
class ConsensusParameter {
   public: 
      /**
       * default constructor
       * /remote/DataAnalysis/metag is full.
       * Old reference data directory full, now switching to a new
       * data directory.
       * default binary is using Hacienda cluster.
       * /remote/RSU/sw-cache/metag/bin
       *
       *  old dataPath("/remote/Overflow/DataAnalysis/metag"),
       *  Use METAG_HOME and METAG_DATA to set up default paths
       *  for the application's reference data directory
       *  and binary path
       */
      ConsensusParameter() 
         : dataPath(),
           binaryPath(),
           refdb("silvaplus")
      { 
         //setDefaultBindir();
         probeFileSystem(); 
      }
      /**
       * Taxonomy file path is simplly the blast reference
       * database with .header suffix
       * Header files has three filed
       * id<TAB>species<TAB>strain "TAB delimited"
       */
      string getTaxFilePath() const { return dataPath + "/refseq/" + refdb + ".header"; }
      /** 
       * @return the reference database
       */
      const string& getReferenceDatabase() const { return refdb; }
      /**
       * Alias to getReferenceDatabase()
       * Any sarch program can be use other than blast.
       */
      const string& getBlastdb() const { return refdb; }
      const string& getRefdb() const { return refdb; }
      /**
       * blastn full path
       */
      string getBlastnPath() const { return binaryPath + "/blastn"; }
      /**
       * return the singleton instance pointer.
       */
      static ConsensusParameter& getInstance() { return onlyOne; }
      ConsensusParameter(const ConsensusParameter& cp)=default;
      ConsensusParameter& operator=(const ConsensusParameter& cp)=default;
      /**
       * Where is the blastn binary
       * @param bp is the directory name for blastn
       */
      void setBinaryPath(const string &bp) { binaryPath=bp; } 
      /**
       * Which directory contains the blast database
       *  Moved to due to out of space 
       *   /remote/Overflow/DataAnalysis/metag
       *   There are three databases to choose from 
       *   rdpplus, silvaplus, silvardp (no addon)
       */
      void setDataPath(const string &dp) { dataPath=dp; }
      /**
       * Right now we are using blast. You can use any
       * program for the reference serach.
       */
      void setRefdb(const string &dbname) { refdb=dbname; }
      /**
       * use the default METAG_HOME/bin
       * This should be the blastn binary location
       */
      void setDefaultBindir();

   private:
      string binaryPath;
      /**
       * This directory contains the silva123.header file
       * usuall $METAG_DATA/refseq/silva1223.header
       * It is a collection of the fasta header of the silva123.fas file.
       * For portability, METAG_DATA environment varable should be used.
       */
      string dataPath;
      /**
       * one of rdpplus, silvaplus, silvardp
       */
      string refdb;
      /**
       * automatically set up the binary path
       * and data path from environmental variable: METAG_DATA
       * and METAG_HOME
       */
      void probeFileSystem();
      static ConsensusParameter onlyOne;
};

/** 
 * There is no way to make the blast to run through SGE
 * -outfmt '6 bla bla foo bar'  SGE agressive remove severyal layers
 *  of quotes. This is for the default parameter of 6 only
 *  headers
 *  'qseqid sseqid pident length mismatch gapopen qstart qend sstart send
 *     evalue bitscore'
 *
 *  The evalue is ignored. Useless for very similar sequences.
 *  The Blastnrow id defined in blastnrow.h it differs from this
 *  one by a little bit.
 */
class Blastnrowd {
   public:
      Blastnrowd() { }
      /**
       * The evalue element is in the input.
       * It will be ignored (discarded).
       */
      Blastnrowd(const vector<string> &row);

      /** helper used by the constructor
       */
      void setTaxon();
      void setCount(int cnt) { queryCount=cnt; }
      /**
       * Set the consensus's depth
       */
      void setDepth(int dep) { consDepth=dep; }
      int getConsensusDepth() const { return consDepth; }
      friend ostream& operator<<(ostream& ous, const Blastnrowd &r);
      int getBitscore() const { return bitscore; }
      /**
       * has a nice scientific name
       */
      bool named() const { return isScientificSpeciesName(taxon.first); }
      /**
       * @return the species string.
       */
      string getSpecies() const { return taxon.first; }
      string::size_type getSpeciesLength() const { return getSpecies().length(); }
      /**
       * @return the species,strain pair. Strain could be empty.
       */
      pair<string,string> getTaxon() const { return taxon; }
      string getMappingName() const {
         if (named()) return getSpecies();
         else return taxon.first + " " + taxon.second;
      }
      /**
       * sort by cluster number from small to large
       */
      bool operator<(const Blastnrowd& r) const;
      /**
       * @return the query name in the format cl###
       */
      const string& getQuery() const { return qseqid; }
      /**
       * @return cluster id as unique identifier
       */
      int getQidNumber() const { return stoi(qseqid.substr(2)); }
      float getIdentity() const { return pident; }
      // the subject length is not available through the default
      // -outfmt 6
      //float getSubjcov() const { return float(send-sstart+1)/slen; }
      int getAlignLength() const { return length; }
      int getQueryCount() const { return queryCount; }
      void merge(const Blastnrowd &row) { queryCount += row.getQueryCount(); 
         consDepth += row.getConsensusDepth(); }
      /**
       * @return a subset of the fields for plotting purpose
       */
      tuple<string, int, string, float> getReduced() const { 
         return make_tuple(getQuery(), getConsensusDepth(), getSpecies(), getIdentity()); 
      }
      /**
       * @return a tuple of cluster_id, depth, taxon combined string
       * If taxon.second (strain) is empty then the taxon string
       * will be the species part taxon.first.
       */
      tuple<string, int, string> getSpeciesDepth() const; 
      vector<string> toStringVector();

      /** making a copy
       * if given as L value
       */
      static void setTaxDictionary(const map<string, pair<string,string> > &id2tax) { 
         taxdic=id2tax; }
      static void setTaxDictionary(map<string,pair<string,string> > &&id2tax) { taxdic=std::move(id2tax); }
      static int getDictionarySize() { return taxdic.size(); }
      /**
       * Will not otuput the END-of-LINE mark so that the caller
       * can append more field to the header. 
       * So you have to output endl to the ouf stream.
       */
      static void writeHeader(ostream &ouf);
      static const vector<string>& getHeader() { return header; }

   private:
      /**
       * Query sequence id, now formated cl##.
       * First two letter fixed, then followed by number from 1 to 
       * whatever.
       */
      string qseqid;
      /**
       * Subject (target) sequence id, or target sequence id
       */
      string sseqid;
      /**
       * Percent identity such as 95.85
       */
      float pident;
      // missing qlen and slen which is present in the 
      // sister class Blastnrow
      int length, mismatch, gapopen, qstart, qend, sstart, send,  bitscore;
      /**
       * Look up result from Silva, rdp, or anyother taxonomy
       * Made of species, strain pair. If the classification
       * is above species levle, then the correspoidng levle will
       * be used.
       */
      pair<string,string> taxon;
      /**
       * number of sequences represented by this one
       */
      int queryCount;
      /**
       * Query sequences consensus vertical depth.
       */
      int consDepth;
      /**
       * Taxonomy dictionary from Silva
       * SequenceId => (species, strain) Taxon string
       */
      //static map<string,string> taxdic;
      static map<string, pair<string, string> > taxdic;
      static vector<string> header;
};


/**
 * for the sorting blast row. Only need bitscore column
 * Othere variations are tolerated.
 */
template<class T>
class sortByBitscore {
   public:
      bool operator()(const T &r1, const T &r2) const {
         return r1.getBitscore() > r2.getBitscore();
      }
};

/**
 * For tabular output:
 *  qseqid sseqid pident qlen slen length qstart qend sstart send  bitscore
 *  There are too many this will take forever to finish on one node.
 *  in one it has typical 13K clusters for noisy large projects.
 *  For cleanner input thinks can be better. You need the algorithm
 *  to get rid of chimeras, none 16S sequences.
 *  The blast only 16S sequences.
 * @param lazy T/F if lazy true, then it will use existing
 *    blast result.
 * @return the blast result file name. use blastnrow.h 
 * for holding the rows of this output.
 * Running blast this way is too slow.
 * This has different headers from the default blast run through
 * the SGE engine where it is having problems with quotes.
 * Note: currently using all available CPU on the SGE slave.
 *   TODO: add a parameter to ConsensParameter to control
 *   the CPU usage.
 *
 * Right now blast is relying on the .ncbirc in the 
 * user's home directory to find the blast database.
 * ~/.ncbirc: file content
 * [BLAST]
 * BLASTDB=/path/to/refseq:/such/as/dataroot/metag/refseq
 * =======================
 */
string runBlast(const string &fastaFile, bool lazy=false);

/** The source is hardcoded. 
 * TODO: make it a parameter for the source file.
 * No longer used. Use loadTaxdic()
 */
map<string,string> slurpSilvaTaxonomy();

/**
 * New version. The job of parsing and processing
 * is done by external modules. This one simply load
 * from the external source. It could be a database
 * or a file. Right now it is a file.
 * it expects SEQID => { species, strain }
 * The strain can be empty.
 * This function will replace slurpSilvaTaxonomy.
 * The taxonomy can be from any source. No longer
 * Silva specific.
 */
map<string, pair<string,string> > loadTaxdic();

template<class T>
vector<T> pickShortestNamedRow(const vector<T> &hits) {
   vector<T> tmp;
   if (hits.size() > 1) { // pick the shortest 
      size_t minj = 0;
      string::size_type minlen=hits[0].getSpeciesLength();
      for (size_t j=1; j<hits.size(); ++j) {
         if (hits[j].getSpeciesLength() < minlen) {
            minj = j;
            minlen = hits[j].getSpeciesLength();
         }
      }
      tmp.push_back(hits[minj]);
   }
   else {
      tmp.push_back(hits[0]);
   }
   return tmp;
}

// template functions must be in the header
/**
 * Picke the best match or hit.
 * If the hit with the highest score is named, then stop 
 * searching for other hits with identical scores.
 * If the hit with the highest score is not named
 * and the first named hit is not more than 10 bitscore
 * less than the top hit, then the top hit will not
 * be included in the report.
 *
 * T is the blastnrow type
 */
template<class T>
vector<T> pickBestRow(vector<T> &hits) { 
   //cerr << hits.size() << " candidate hits for pickBestRow\n";
   vector<T> tmp;
   sort(hits.begin(), hits.end(), sortByBitscore<T>());
   string species; // tmp variable
   if (hits[0].named()) { 
      //cerr << "picking only one top named hit\n";
      //cerr << hits[0] << " is named\n";
      set<string> discovered;
      species=hits.front().getSpecies();
      tmp.push_back(hits.front());
      discovered.insert(species);
      //cerr << "top hit species: " << species << endl;
      tmp.push_back(hits.front());
      int bits=hits.front().getBitscore();
      size_t j=1;
      while (j<hits.size() && hits[j].getBitscore() == bits) {
         if (hits[j].named()) {
            species=hits[j].getSpecies();
            if (discovered.find(species) == discovered.end()) {
               //cerr << "identical score hit species: " << species << endl;
               tmp.push_back(hits[j]);
               discovered.insert(species);
            }
         }
         ++j;
      }
      if (tmp.size() > 1) { // pick the shortest 
         tmp = pickShortestNamedRow(tmp);
      }
   }
   else { // if there is at least one name then use it.
      bool hasNamed=false;
      bool saveUnnamed=false;
      //cerr << "first not named: " << hits[0].getSpecies() << endl;
      for (size_t i=1; i<hits.size(); ++i) { // search the rest
         if (hits[i].named()) {
            hasNamed=true;
            tmp.push_back(hits[i]);
            //cerr << "found one named " << hits[i].getSpecies() << endl;
            // only if the largest bitscore is > 10 points
            if ((hits[i].getBitscore() + 10) < hits.front().getBitscore()) {
               saveUnnamed=true;
               //cerr << "no need to save unnamed\n";
            }
            for (size_t j=i+1; j<hits.size() && hits[j].getBitscore() == hits[i].getBitscore(); ++j) {
               //cerr << hits[j].getSpecies() << " " << hits[j].getTaxon().second << endl;
               if (hits[j].named()) {
                  //cerr << " is named\n";
                  tmp.push_back(hits[j]);
               }
               //else {
               //   cerr << " is unnamed, not saving it\n";
               //}
            }
            //cerr << "after collecting named hits\n";
            //writeHits(tmp, cerr);
            if (tmp.size()>1) {
               //cerr << "more than one taxon picked\n";
               tmp = pickShortestNamedRow(tmp);
            }
            break;
         }
      }
      if (!hasNamed || saveUnnamed) {
         tmp.push_back(hits[0]);
      }
   }
   //cerr << "final rows picked\n";
   //writeHits(tmp, cerr);
   return tmp;
}

/**
 * @param hits is the input, after the run it will be the
 *       bad ones (discarded hits)
 * @param idencut identity cut, 96.9 is the default
 * @param covcut subject coverage cut 0.8 is the default
 *       the coverage cutoff is determined by the amplicon size 
 *       relative to the full length gene.
 * @return the good hits
 */
template<class T>
vector<T> filterHit(vector<T> &hits, float idencut=96.9, int lencut=1200) {
   vector<T> good, bad;
   for (size_t i=0; i<hits.size(); ++i) {
      if (hits[i].getIdentity() < idencut || hits[i].getAlignLength() < lencut) {
         bad.push_back(hits[i]);
      }
      else good.push_back(hits[i]);
   }
   hits=bad;
   return good;
}

template<class T>
list<T> mergeHit(const vector<T> &hits) {
   list<T> tmp(hits.begin(), hits.end());
   typename list<T>::iterator it, itt, del;
   it=tmp.begin();
   while (it != tmp.end()) {
      itt=it;
      ++itt;
      if (itt == tmp.end()) break;
      while (itt != tmp.end()) {
         //if (it->getTaxon().first == itt->getTaxon().first) {
         if (it->getMappingName() == itt->getMappingName()) {
            it->merge(*itt);
            del=itt;
            ++itt;
            tmp.erase(del);
         }
         else ++itt;
      }
      ++it;
   }
   return tmp;
}

template<class T>
void writeHits(const vector<T> &hit, ostream &ous) {
   for (size_t i=0; i<hit.size(); ++i) {
      ous << hit[i] << endl;
   }
}

template<class T>
vector<vector<string> > buildUniqueSortedHits(list<T> &merged) {
   // remove duplicated hits, higher by depth
   list<T> res;
   T last;
   typename list<T>::iterator it = merged.begin();
   res.push_back(*it);
   last = *it;
   ++it;
   while (it != merged.end()) {
      if (it->getQuery() == last.getQuery()) {
         if (it->getConsensusDepth() >= last.getConsensusDepth()) { // use this one
            res.back() = *it;
         }
      }
      else {
         res.push_back(*it);
      }
      last = *it;
      ++it;
   }
   res.sort([](T &a, T &b)->bool{ return a.getConsensusDepth() > b.getConsensusDepth(); });
   vector<vector<string> > vvs;
   for (it = res.begin(); it != res.end(); ++it) {
      vvs.push_back(it->toStringVector());
   }
   return vvs;
}

/**
 * write the result to a file
 * @param blastFile input of the tabular blast
 * @param id2cnt the number of sequences for each consensus.
 * @param mappingFile is the output file for mapped OTU and details of mapping.
 *        Two other files that are based on the mapped file are BAD and MERGED
 *        files: *.bad.tab and *.goodmerged.tab that differ only in the suffix.
 * @param idencut identity cut off
 * @param covcut subject horizontal coverage cutoff, default 0.8
 *    This is accuming full length amplicon. This parameter needs
 *    to be adjusted according to the amplicon size.
 * @return consens sequence id mapped to refdb
 */
template<class T>
set<string> writeTaxMapping(const string &blastFile, 
      const map<string, pair<int,int> > &id2cnt, 
      const string &mappingFile, float idencut=96.9, int lencut=1200)
{
   cerr << "writing otu taxonomy mapping to file ...\n";
   //T::setTaxDictionary(slurpSilvaTaxonomy());
   T::setTaxDictionary(loadTaxdic());
   ifstream inf(blastFile);
   ofstream ouf(mappingFile);
   T::writeHeader(ouf);
   ouf << endl;
   string line;
   getline(inf, line); // blast has no header
   vector<string>::iterator it;
   vector<T> all;
   vector<string> row = split(line, '\t');
   while (!inf.eof()) {
      string query=row[0];
      //cerr << "working on query " << query << endl;
      vector<T> hits; // all hits for one query
      while (!inf.eof() && query==row[0]) {
         hits.push_back(T(row));
         getline(inf, line);
         if (inf.eof()) break;
         row = split(line, '\t');
      }
      //cerr << "collecting the best hits\n";
      vector<T> bestHit=pickBestRow(hits);
      pair<int,int> clsize=id2cnt.find(query)->second;
      for (size_t i=0; i<bestHit.size(); ++i) {
         bestHit[i].setCount(clsize.first);
         bestHit[i].setDepth(clsize.second);
      }
      all.insert(all.end(), bestHit.begin(), bestHit.end());
   }
   //cerr << "all queries processed\n";
   sort(all.begin(), all.end()); 
   writeHits(all, ouf);
   cout << "taxonomy otu mapping written to " << mappingFile << endl;
   ouf.close();
   vector<T> good = filterHit(all, idencut, lencut);
   // all now contains bad
   set<string> mappedConsensus;
   for (unsigned i=0; i < good.size(); ++i) {
      mappedConsensus.insert(good[i].getQuery());
   }
   string filteredOut = mappingFile.substr(0, mappingFile.rfind('.')) + ".bad.tab";
   ouf.open(filteredOut);
   T::writeHeader(ouf);
   ouf << endl;
   writeHits(all, ouf);
   cout << "mapping filtered out written to " << filteredOut << endl;
   list<T> merged = mergeHit(good);
   ouf.close();
   string mergedFile = mappingFile.substr(0, mappingFile.rfind('.')) + ".goodmerged.tab";
   ouf.open(mergedFile);
   T::writeHeader(ouf);
   ouf << endl;
   // drawing graphs of various kinds
   vector<tuple<string,int,string,float> > boxdata;
   vector<tuple<string,int,string> > piedata;
   for (auto it = merged.begin(); it != merged.end(); ++it) {
      boxdata.push_back(it->getReduced());
      piedata.push_back(it->getSpeciesDepth());
      ouf << *it << endl;
   }
   BoxPile boxchart(boxdata);
   string boxfile="otubar.svg";
   boxchart.draw(boxfile);
   string boxfile2="otubar2.svg";
   boxchart.draw2(boxfile2);
   string boxfile3="otubar3.svg";
   boxchart.draw3(boxfile3);
   PieChart piechart(piedata);
   piechart.draw();
   cout << "Piechart done\n";
   cout << "merged good out mapping written to " << mergedFile 
      << " BoxChart saved in " << boxfile << " " << boxfile2
      << " " << boxfile3 << endl;
   TextPile::setHeader(T::getHeader());
   TextPile textpile(buildUniqueSortedHits(merged));
   textpile.draw("top10hits.png");
   return mappedConsensus;
}

/**
 * Version take vector<DNAQualCount>& in place of map<string,int>
 * @param consq consensus sequence in DNAQualCount format.
 */
template<class T>
void writeTaxMapping(const string &blastFile, 
      const vector<DNAQualCount> &consq, const string &mappingFile) {
   map<string,int> tmp;
   for (auto it=consq.begin(); it != consq.end(); ++it) {
      tmp.insert(make_pair(it->getName(), it->getCount()));
   }
   writeTaxMapping<T>(blastFile, tmp, mappingFile);
}

#endif
