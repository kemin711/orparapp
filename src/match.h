// match.h
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/* a hit is a summary of all the matches between each pair of 
 * query-target, only non-overlapping matches will be registered
 * for a hit.
 * According to NCBI blast -m 8 output has the following fields:
 * Query id, Subject id, % identity, alignment length, mismatches, gap
 * openings, query start, query end, sequence start, sequence end, e-value, bit
 * score
 *
 */

/* A match does not have an identifier: query-target.  It is all the 
 * numerical part of the blast -m 8 output
 */
class inputException {};

class match
{
	public:
		match(istream &in) throw(inputException);
		match(const match &m);
		match& operator=(const match &m);
		friend ostream& operator<<(ostream &o, const match &m);
		bool overlap(const match &m) const;
		float getScore() const { return score; }
		int getLength() const { return length; }
		float getLenxid() const { return length*identity; }
		double getE() const { return E; }

	private:
		float identity;
		int length, mismatch, gap, qstart, qend, tstart, tend;
		double E;   // exp from blast
		float score;
};

/* contains one or matches. A summary of all the matches between 
 * the query and the target in the database
 * */
class hit
{
	public:
		class end {};
		/* to read the first line */
		hit() : score(0), matches(), matchCnt(0) {}
		//hit(istream &in);
		// used q, and t to construct and will set the next q,t 
		hit(string &q, string &t, istream &in) throw(end, inputException);

		/* append the next match to this hit */
		void getNext(istream &in);
		bool sameQuery(const string &nq) const { return query==nq; }
		bool sameTarget(const string &nt) const { return target==nt; }
		bool sameHit(const string &q, const string &t) const;
		bool overlap(const match &m) const;
		/* r < 1.  this hit is 0.9x of h's sumScore */
		bool scoreAsBigAs(const hit &h, float r) const { return score > r*h.score; }
		/* output the hit result for loading into a table
		 * query,target,average_identity,sum_match_length,sum_score,selratio,
		 * minE, prodE(products of all non-overlapping E)
		 * */
		friend ostream& operator<<(ostream &o, const hit &h);
		int getLength() const { return length; }

		/* non-overlapping matches */
		int getMatchCount() const { return matches.size(); }
		int getTotalMatch() const { return matchCnt; }
		void dumpRaw(ostream &o) const;
		float getIdentity() const { return identity; }
		float getScore() const { return score; }
		// weighted average

		static char fields[];

	private:
		string query, target;
		vector<match> matches; //non-overlapping matches
		float score;   // sum of score from all non-overlapping matches
		int matchCnt;  // total match inputed
		int length;    // sum of length form all matches
		float identity; // identity averaged over all matches
		double minE;   // minimus of all the E values from all matches
		double prodE;  // product of all non-overlapping E's 
};

