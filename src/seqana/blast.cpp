//file blast.cpp
#include "blast.h"
#include <ctype.h>
#include <cstdlib>

using namespace std;

void split(char ln[], int &s, string &seq, int &e);

void alnsummary::read(const char ln[], int ss)
//the new blast program 2000 version have changed the output format
//no number of segments were seen with the summary line
//at lest in Tblastx 
{ //process one line of the summary
	const char *p = strchr(ln, ' ');
	subj.assign(ln, p-ln);
	char tmp[50];
	p = ln+ss;
	while (isspace(*p)) p++;
	strcpy(tmp, p);
	char *pp = strchr(tmp, ' ');
	*pp = '\0';
	s = atoi(tmp);
	//p = ++pp;
   char *xp = ++pp;
	while (isspace(*xp)) xp++;  //now p points to E-value
	pp = strchr(xp, ' ');
	if (pp) {
		*pp = '\0'; //for ungapped blast
		e = atof(xp);
		pp++;
		n = atoi(pp);
	}
	else e = atof(xp);  //for gaped or tblastx no n exist
}

string alnsummary::subjName() const
{
	string tmp;
	if (subj.substr(0, 3) == "gi|") {  
	//special naming for nr gi| protein
	//assuming everything start with gi|
		tmp = subj.substr(3);
		string::size_type i = tmp.find('|');
		if (i != string::npos) tmp.erase(i);
		tmp.insert(0, "nr:");
		return tmp;
	}
	return subj;
}

///////////////alignment////////////////////////////////////

int countIden(const char ln[])
{
//count identical and similar residues
//this function is mainly for the bugy older version of 
//blast.  This problem seem to have been corrected in the new version
//of blast 2000
	int cnt=0, i=0;
	while (ln[i] != '\0') if (isalpha(ln[i++])) cnt++;
	return cnt;
}

int alignment::read(istream &ins, char ln[], string &bf)
{ //text info from Score = ..... to the end of the alignment
//will be stored in bf
//the Score = line is initially stored in ln
//read up to BLASTN or BLASTX if reading the last alignment
//if reached the end of file, return 0;
//note: !!!!!!!!!!!!!!!!!
//the blast version older gave errorneous results 
//as to length and identity
	char *p, *pp;

	bf = ln;
	bf += "\n";
	//ln contains the  Score = 45.0 bits (92), Expect = 0.005

	p = ln + 9;
	while (*p == ' ') p++;
	pp = strchr(p, ' ');
	*pp = '\0';
	sco = atof(p);
	pp += 8;
	p = strchr(pp, '=') + 2;
	exp = atof(p);  //Expect
	ins.getline(ln, LINE); // Identities = ......
	bf += ln;
	bf += "\n";

//the blast version gave errorneous results 
//the 2000 version has corrected this problem
	p = ln + 14;
	pp = strchr(p, '/');
	*pp = '\0';
	iden = atoi(p);
	p = pp+1;
	pp = strchr(p, ' ');
	*pp = '\0';
	len = atoi(p);

/* the blastn does not have similarity info !! 
*/
	if (p = strchr(pp, '=')) {  //not parsing blastn
		p += 2;
		pp = strchr(p, '/');
		*pp = '\0';
		simil = atoi(p);
	}
	//ins.getline(ln, LINE);  // frame = -2 / -1  Blastx
	//or Strand = Plus / Plus in blastn
	//this line will be discared because it contains redundant 
	//information

	ins.getline(ln, LINE);
	while (strncmp(ln, "Query:", 6)) { //discard empty lines
		bf += "\n";
		ins.getline(ln, LINE); 
	} //read up to Query:
	bf += ln;
	bf += "\n";
	split(ln, qs, qseq, qe);
	ins.getline(ln, LINE); //line with ||||||||||||
	//iden += countIden(ln); //count identical and similar residues
	//only idencal residues are counted now

	bf += ln;
	bf += "\n";
	ins.getline(ln, LINE); //Sbjct: ......
	bf += ln;
	bf += "\n";
	split(ln, ss, sbjseq, se);
	ins.getline(ln, LINE); //empty line
	while (ln[0] != '>' && strncmp(ln, " Score =", 8) &&
   strncmp(ln, "BLAST", 5) && strncmp(ln, "TBLAST", 6) 
	&& strncmp(ln, "Query:", 6) && !ins.eof()) {
		bf += ln;
		bf += "\n";
		ins.getline(ln, LINE);
	}
	if (ins.eof()) return 0;
	while (!strncmp(ln, "Query:", 6)) {  //alignment more than one line
		bf += ln;
		bf += "\n";
		int tmps;
		string tmpseq;
		split(ln, tmps, tmpseq, qe);
		qseq += tmpseq;
		ins.getline(ln, LINE);  //line with ||||||
		iden += countIden(ln);
		bf += ln;
		bf += "\n";
		ins.getline(ln, LINE);
		bf += ln;
		bf += "\n";
		split(ln, tmps, tmpseq, se);
		sbjseq += tmpseq;
		ins.getline(ln, LINE);
		while (ln[0] != '>' && strncmp(ln, " Score =", 7) 
		&& strncmp(ln, "BLAST", 5) && strncmp(ln, "TBLAST", 6) 
		&& strncmp(ln, "Query:", 6) && !ins.eof()) {
			bf += ln;
			bf += "\n";
			ins.getline(ln, LINE);
		}
		if (ins.eof()) return 0;
	}
	//len = qseq.length();  //information obtained from 
	//header of each alignment; only works with the new
	//2000 version of blast
	return 1;
}

void split(char ln[], int &s, string &seq, int &e)
{
	char *p, *pp;
	p = strchr(ln, ' ') + 1;
	pp = strchr(p, ' ');
	*pp = '\0';
	s = atoi(p);
	p = pp + 1;
	while (isspace(*p)) p++;
	pp = strchr(p, ' ');
	*pp = '\0';
	seq = p;
	pp++;
	e = atoi(pp);
}

void seq2comp(const string &s, double v[])
{
	//char aa[] = "ARNDCQEGHILKMFPSTWYVBZX*";
	//same order as in matrix

	int i=0;
	for (i=0; i<24; i++) v[i]=0; //initialize vector to zero counts

	for (i=0; i<s.length(); i++) {
		switch(s[i]) {
		case 'A':
		case 'a': v[0]++; break;
		case 'R':
		case 'r': v[1]++; break;
		case 'N':
		case 'n': v[2]++; break; 
		case 'D':
		case 'd': v[3]++; break;
		case 'C':
		case 'c': v[4]++; break;
		case 'Q':
		case 'q': v[5]++; break;
		case 'E':
		case 'e': v[6]++; break;
		case 'G':
		case 'g': v[7]++; break;
		case 'H':
		case 'h': v[8]++; break;
		case 'I':
		case 'i': v[9]++; break;
		case 'L':
		case 'l': v[10]++; break;
		case 'K':
		case 'k': v[11]++; break;
		case 'M':
		case 'm': v[12]++; break;
		case 'F':
		case 'f': v[13]++; break;
		case 'P':
		case 'p': v[14]++; break;
		case 'S':
		case 's': v[15]++; break;
		case 'T':
		case 't': v[16]++; break;
		case 'W':
		case 'w': v[17]++; break;
		case 'Y':
		case 'y': v[18]++; break;
		case 'V':
		case 'v': v[19]++; break;
		case 'B':
		case 'b': v[20]++; break;
		case 'Z':
		case 'z': v[21]++; break;
		case 'X':
		case 'x': v[22]++; break;
		case '*': v[23]++; break;
		default: 
			cerr << v[i] << " is not an aa symbol, make sure the blast \n"
			<< " is run with -g F [no gap]\n";
			exit(1);
			break;
		}
	}
}

double alignment::calBias(double &efflen, int psw) const
{
//calculate bias index
//psw is the pseudo weight, 5 was found to be good
//lower weight tend to reject too many short true matches
//higher weight tend to accept too many bioased composition matches
//score <0.8 is found to be biased for most proteins

	double qv[24]={0}, dv[24]={0};

   //char aa[] = "ARNDCQEGHILKMFPSTWYVBZX*";
	//same order as in matrix
	seq2comp(qseq, qv);
	seq2comp(sbjseq, dv);
	efflen = len - qv[20] - qv[21] - qv[22];  //effective length
	//after discounting XXX BBBB ZZZ 
	int i;
	for (i=0; i<20; i++)  qv[i]=(qv[i]+psw*aafq[i])/(len+psw);
	for (i=0; i<20; i++)  dv[i]=(dv[i]+psw*aafq[i])/(len+psw);

	double eps=0;  //expected score
	for (i=0; i<20; i++) 
		for (int j=0; j<20; j++) eps += qv[i]*dv[j]*MATRIX[i][j];
	eps *= len;
	return (sco-eps)/(sco-len*Mexp);
}
	

//////////////////region//////////////////////
region::region(const alignment &al)
{
	s = al.qs;
	f = al.qe;
	hs = al.sco;
	tc = 1;
	exp = al.exp;
}

void region::assign(int be, int en, float sc, double ex)
{
	s = be;
	f = en;
	hs = sc;
	tc = 1;
	exp = ex;
}

bool region::sameDirection(alignment &r) const
{
	if (isup()) {
		if (r.isup()) return true;
		else return false;
	}
	else { //this is lower
		if (r.isup()) return false;
		else return true;
	}
}

int region::cover(region &r) const
{
   if (isup()) {
      if (r.isup()) {
			if (r.f < s || r.s > f) return 0;
			else {
				if (r.s > s && r.f > f) return 3;
				else if (r.s < s && r.f < f) return 1;
				else return 2;
			}
		}
		else return 0;
   }
   else { //this is lower
      if (r.isup()) return 0;
		else { //r is lower strand
			if (r.f > s || r.s < f) return 0;
			else {
				if (r.s > s && r.f > f) return 1;
				else if (r.s < s && r.f < f) return 3;
				else return 2;
			}
		}
   }
}

int region::distance(const alignment &aa) const
{
	region rr(aa);
	if (!cover(rr)) {
		int dist;
		if (isup()) {
			if (rr.isup()) {
				dist = f-rr.s;
				if (dist < 0) dist = rr.f - s;
				return dist;
			}
			else { //on different strands, and infinite long distance
				return 1000000; //very long distance in DNA sequence
			}
		}
		else {
			if (rr.isup()) return 1000000;
			else {
				dist = f-rr.s;
				if (dist < 0) return (rr.f-s);
			}
		}
	}
	else return -1;  //use -1 to represent coverage
}

int region::add(const alignment &r, bool &accept, double scut,
	double ecut, int cl)
{ 
	//score cut for segment should be very low
   if (isup()) {
      if (r.isup()) {
			if (r.qe < s || r.qs > f) {
				accept = false;
				return 0;
			}
			else {
				if (r.qs > s && r.qe > f) { //r overlap to the end
					if (tc < cl) {  //not exceeding covering limit
						if (r.sco > scut || r.exp < ecut) {
							if (hs - r.sco < 6 || r.sco > 50) {
								f = r.qe;
								tc++;
								accept =  true;
								return 3;
							}
						}
					}
					//else {  //exceeding covering limit
					//but in tblastx, segments following may have higher score!
					accept = false;
					return 3;
				}
				else if (r.qs < s && r.qe < f) {
					if (tc < cl) {
						if (r.sco > scut || r.exp < ecut) {
							if (hs- r.sco < 6 || r.sco > 50) {
								s = r.qs;
								tc++;
								accept = true;
								return 1;
							}
						}
					}
					accept = false;
					return 1;
				}
				else { //complete covers r
					if (tc < cl) {
						if (r.sco > scut || r.exp < ecut) {
							if (hs - r.sco < 6 || r.sco > 50) {
								tc++;
								accept = true;
								return 2;
							}
						}
					}
					accept = false;
					return 2;
				}
			}
		}
		else {
			accept = false;
			return 0;
		}
   }
   else { //this is lower.qstrand
      if (r.isup()) {
			accept = false;
			return 0;
		}
		else { //r is lower.qstrand
			if (r.qe > s || r.qs < f) {
				accept = false;
				return 0;
			}
			else {
				if (r.qs > s && r.qe > f)  {
					if (tc < cl) {
						if (r.sco > scut || r.exp < ecut) {
							if (hs - r.sco < 6 || r.sco > 50) {
								s = r.qs;
								tc++;
								accept =  true;
								return 1;
							}
						}
					}
					accept = false;
					return 1;
				}
				else if (r.qs < s && r.qe < f) {
					if (tc < cl) {
						if (r.sco > scut || r.exp < ecut) {
							if (hs - r.sco < 6 || r.sco > 50) {
								f = r.qe;
								tc++;
								accept =  true;
								return 3;
							}
						}
					}
					accept = false;
					return 3;
				}
				else { 
					if (tc < cl) {
						if (r.sco > scut || r.exp < ecut) {
							if (hs - r.sco < 6 || r.sco > 50) {
								tc++;
								accept = true;
								return 2;
							}
						}
					}
					accept = false;
					return 2;
				}
			}
		}
   }
}
