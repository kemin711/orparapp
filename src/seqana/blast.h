#ifndef ALNSUMMARY_H
#define ALNSUMMARY_H
#include <string>
#include <cstring>
#include <iostream>

#ifndef LINE
#define LINE 90
#endif

//const int  COVERAGE 24  //sets the limit of coverage in one region
//defined in main program

//#define SCORE_THRESHOLD  40  //only scores larger than this is considered
//for extending an coverage
//#define E_THRESHOLD 0.3  

#define MATRIX BLOSUM62
#define Mexp -0.945

using namespace std;

extern const int PAM120[24][24]; 
extern const int BLOSUM62[24][24]; 
extern const int BLOSUM62_12[24][24];
extern const double aafq[20];

//for the header portion of blast results
class alnsummary
{
 public:
	int score() const {return s;}
	bool scoreHigherThan(double ss) const {return s>ss;}
	bool scoreLessThan(double ss) const {return s<ss;}
	void read(const char ln[], int ss);
	bool eLessThan(double ee) const {return e < ee;}
	bool eMoreThan(double ee) const {return e > ee;}
	string subjName() const ; //return subject Name

 private:
	string subj;  //target sequence name
	int s;  //score
	double e;  //E Value
	int n;  //number of matching segments
};

void split(char ln[], int &s, string &seq, int &e);
//for the alignment portion of blast results
//starting from the Score = line
class alignment
{
 public:
	int read(istream &ins, char ln[], string &bf);
	//will read upto BLASTX(N) or end
	//return 0 when reaching the end
	//old contents of bf will be flushed away
	int begin() const {return qs;}
	int end() const {return qe;}
	int sbjbegin() const {return ss;}
	int sbjend() const {return se;}
	float score() const {return sco;}
	double eValue() const {return exp;}
	int length() const {return len;}
	int identical_residue() const {return iden;}
	double calBias(double &efflen, int psw = 6) const;  
	//calculate the bias-ratio, efflen is the len minus XXX ZZZ BBB
	string querySeq() const {return qseq;}
	string matchSeq() const {return sbjseq;}
	float identity() const {return (float)iden/len;}
	bool isup() const {return qe>qs;}
	bool iscompl() const {return qs>qe;} //lower strand of DNA

friend class region;

 private:
 	//friend void split(char ln[], int &s, string &seq, int &e);
   //int expScore(); //calculate expected score
	float sco;  //score
	double exp;
	int iden; //number of identical residues
	int simil; //number of similar residues
	int len; //length of the alignment
	string qseq;
	string sbjseq;
	int qs;  //query start position
	int qe;  //query end position
	int ss;
	int se;
};

class region
{
 public:
	region(const alignment &al);  //copy constructor
 	void assign(int be, int en, float sc, double ex);
	int cover(region &r) const; //return 0 if no overlap with r
	//return 1 if overlap start, return 3 if overlap end
	//return 2 if completely covered; return 4 if not determined
	int add(const alignment &r, bool &accept, double scut, 
		double ecut, int cl); 
	//same as cover except action will be taken
	//to consolidate two regions of hit if they overlap, otherwise
	//no action will be taken
	bool isup()  const {return f>s;}
	bool iscompl() const {return s>f;} //lower strand of DNA
	bool sameDirection(alignment &r) const;
	int distance(const alignment &aa) const;
	//returns the distance between this region and alignment aa
	//returns -1 if overlaping

 private:
	int s, f; //start, finish
	float hs; //highest score; or the score for a new region
	double exp;
	int tc; //times covered
};

#endif
