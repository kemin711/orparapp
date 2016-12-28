#include "match.h"
/////////// the match class /////////////////////////
ostream& operator<<(ostream &o, const match &m) {
	o << m.identity << '\t' << m.length << '\t' << m.mismatch << '\t' 
		<< m.gap << '\t' << m.qstart << '\t' << m.qend << '\t' 
		<< m.tstart << '\t' << m.tend << '\t' << m.E << '\t' << m.score;
	return o;
}

match::match(const match &m) {
		identity = m.identity;
		length = m.length;
		mismatch = m.mismatch;
		gap = m.gap;
		qstart = m.qstart;
		qend = m.qend;
		tstart = m.tstart;
		tend = m.tend;
		E = m.E;
		score = m.score;
}

match& match::operator=(const match &m) {
	if (this != &m) {
		identity = m.identity;
		length = m.length;
		mismatch = m.mismatch;
		gap = m.gap;
		qstart = m.qstart;
		qend = m.qend;
		tstart = m.tstart;
		tend = m.tend;
		E = m.E;
		score = m.score;
	}
	return *this;
}

match::match(istream &in) throw(inputException) {
	in >> identity >> length >> mismatch >> gap >> qstart >> qend
		>> tstart >> tend >> E >> score;
	if (in.fail()) {
		throw inputException();
	}
}

bool match::overlap(const match &m) const {
	static const int o=1; // overlap  margin --- --
	if ( (m.qstart >= qstart && m.qstart < qend-o) ||
			(m.qend > qstart+o && m.qend <= qend) ||
			 (m.tstart >= tstart && m.tstart < tend-o) || 
			 (m.tend > tstart+o && m.tend <= tend) || 
			 (m.qstart <= qstart && m.qend >= qend) ||
			 (m.tstart <= tstart && m.tend >= tend) ) return true; 
	return false;
}

/////////////////////////////////////////////////////
////////// hit class ////////////////////////////////

/* not used anymore 
hit::hit(istream &in) {
	in >> query >> target;
	matches.push_back(match(in));
	score =matches[0].getScore();
}
*/
hit::hit(string &q, string &t, istream &in) throw(end, inputException) {
	query = q;
	target = t;
	try {
		matches.push_back(match(in));
		matchCnt = 1;
		score  = matches[0].getScore();
		length = matches[0].getLength();
		identity = matches[0].getLenxid();
		minE = matches[0].getE();
		prodE = matches[0].getE();
		in >> q;
		if (in.eof()) throw end();
		in >> t;  // excception: reached the end of stream
		while (query == q && target == t) {
			matchCnt++;
			getNext(in);
			in >> q;
			if (in.eof()) throw end();
			in >> t;
		}
		if (in.fail()) throw inputException();
	}
	catch(inputException) {
		cerr << "input exception at: " << query << " " << target << endl;
		throw inputException();
	}
	identity = identity/length;
	if (minE < 6e-171)  minE = 0;
	if (prodE < 6e-171) prodE = 0;
}

bool hit::overlap(const match &m) const {
	for (int i=0; i<matches.size(); i++) 
		if (m.overlap(matches[i])) 
			return true;
	return false;
}

void hit::getNext(istream &in) {
	match tmp(in);
	if (!overlap(tmp)) {
		score  += tmp.getScore();
		length += tmp.getLength();
		identity += tmp.getLenxid();
		matches.push_back(tmp);
		if (tmp.getE() < minE) minE = tmp.getE();
		prodE *= tmp.getE();
	}
}

/* for output into relational database */
ostream& operator<<(ostream &o, const hit &h) {
	float selratio = (float)(h.getMatchCount())/h.getTotalMatch();
	o << h.query << '\t' << h.target << "\t" << h.getIdentity() << '\t'
		<< h.getLength() << '\t' << h.score << '\t'
		<< selratio << '\t' << h.minE << '\t' << h.prodE;

		//<< h.getMatchCount() << '\t' << h.getTotalMatch();
	return o;
}

/*
int hit::getLength() const {
	int len=0;
	for (int i=0; i<matches.size(); i++) {
		len += matches[i].getLength();
	}
	return len;
}
*/

void hit::dumpRaw(ostream &o) const {
	for (int i=0; i<matches.size(); i++) {
		o << query << '\t' << target << '\t' << matches[i] << endl;
	}
}
/*
float hit::getIdentity() const {
	float sum=0;
	for (int i=0; i<matches.size(); i++) {
		sum += matches[i].getLenxid();
	}
	return sum/getLength();
}
*/

char hit::fields[] = "query, target, average_identity, sum_match_length, sum_score, selratio, minE, prodE(products of all non-overlapping E)";
