#include "nucaln.h"
#include <iostream>

const char Nucaln::gapchar='-';

pair<string, string> Nucaln::getConsensus() {
   if (length() < len) {
      return make_pair(string(), string());
   }
   if (length() > len) {
      //cout << "length() " << length() << " len: " << len << endl;
      size_t x = length() - 24;
      if (x > 0 && bottom.substr(0,x) == string(x, 'A') && top.substr(0,x) == string(x, gapchar)) {
         top = top.substr(x);
         bottom = bottom.substr(x);
      }
      /*
      if (length() == 26 && bottom.substr(0,2) == "AA" 
            && top.substr(0,2) == string(2, gapchar)) {
         //cout <<"3' right flnak confused by the Poly T issue\n";
         top = top.substr(2);
         bottom = bottom.substr(2);
      }
      */
      else if (length() == 26 && top.substr(0,3) == "GGG"
            && bottom.substr(0,3) == "--G") {
         top = top.substr(2);
         bottom = bottom.substr(2);
      }
      else {
         findCandidate();
         if (hasCandidate()) {
            if (fixLength()) {
               //cout << "length problem can be fixed\n";
            }
            else {
               //cout << "Cannot fix length problem\n";
               return make_pair(string(), string());
            }
         }
         else { // no way to make sequence shorter!
            // return original
            return make_pair(top, bottom);
         }
      }
   }
   string tt(top), bb(bottom);
   int diff = 0;
   for (size_t i=0; i<tt.size(); ++i) {
      if (tt[i] == gapchar) {
         tt[i] = bb[i];
      }
      else if (bb[i] == gapchar) {
         bb[i] = tt[i];
      }
      else if (bb[i] != tt[i]) ++diff;
   }
   if (diff > 0) {
      return make_pair(tt, bb);
   }
   else {
      return make_pair(tt, "");
   }
}

string Nucaln::getUngappedTop() const {
   string tmp;
   for (size_t i=0; i<top.size(); ++i) {
      if (top[i] != gapchar) {
         tmp += top[i];
      }
   }
   return tmp;
}

string Nucaln::getUngappedBottom() const {
   string tmp;
   for (size_t i=0; i<bottom.size(); ++i) {
      if (bottom[i] != gapchar) {
         tmp += bottom[i];
      }
   }
   return tmp;
}

int Nucaln::getDiff() const {
   return length() - len;
}

bool Nucaln::fixLength() {
   int diff = getDiff();
   //cout << diff << " longer than expected\n apptemping to fix:\n";
   //cout << top << endl << bottom << endl;

   if (diff > 0 && getNumCandidate() >= diff) {
      for (int i = 0; i<diff; ++i) {
         CandidateRemoval x = candidate.top();
         if (x.getIndex() >= top.length()) {
            cerr << x.getIndex() << " ouside aln length: " << top.length() << endl;
            exit(1);
         }

         top.erase(x.getIndex(), 1);
         bottom.erase(x.getIndex(), 1);
         candidate.pop();
         updateCandidate(x.getIndex());
      }

      //cout << "modified aln:\n"
      //   << top << endl << bottom << endl;
      return true;
   }
   else return false;
}


void Nucaln::updateCandidate(const size_t fx) {
   // less efficient, we only have one or two element, which is fine
   vector<CandidateRemoval> tmp;
   while (!candidate.empty()) {
      CandidateRemoval c = candidate.top();
      if (c.getIndex() > fx) {
         c.lowerIndex();
      }
      tmp.push_back(c);
      candidate.pop();
   }
   for (size_t i=0; i<tmp.size(); ++i) {
      candidate.push(tmp[i]);
   }


   /* cannot walk through a priority queue!
   priority_queue<CandidateRemoval>::iterator it = candidate.begin();
   while (it != candidate.end()) {
      it->lowerIndex();
      ++it;
   }
   */
}

void Nucaln::findCandidate() {
   //priority_queue<CandidateRemoval> candidate;

   // scan the bottom strand for potential problem
   size_t b, t;
   for (t = 0; t < top.size(); ++t) {
      if (top[t] == gapchar) {
         b = t+1; // scan forward
         int count = 1;
         while (b < bottom.size() && bottom[t] == bottom[b]) {
            ++count;
            ++b;
         }
         if (count > 1) {
            //cout << count << " potential error on bottom strand\n";
            candidate.push(CandidateRemoval(false, t, count));
         }
         else { // scan backward
            b = t - 1;
            while (b > -1 && bottom[t] == bottom[b]) {
               ++count;
               --b;
            }
            if (count > 1)  {
               candidate.push(CandidateRemoval(false, t, count));
            }
         }
      }
   }

   // scan the top aln for potential problems
   for (b = 0; b < bottom.size(); ++b) {
      if (bottom[b] == gapchar) {
         t = b+1;
         int count = 1;
         while (t < top.size() && top[b] == top[t]) {
            ++count;
            ++t;
         }
         if (count > 1) {
            //cout << count << " potential error on top strand\n";
            candidate.push(CandidateRemoval(true, b, count));
         }
         else {
            t = b - 1;
            while (t > -1 && top[b] == top[t]) {
               ++count;
               --t;
            }
            if (count > 1) {
               candidate.push(CandidateRemoval(true, b, count));
            }
         }
      }
   }
}

ostream& operator<<(ostream &ous, const Nucaln &aln) {
   ous << aln.top << endl << aln.bottom << endl;
   return ous;
}

Nucaln::Nucaln(const Nucaln &nca)
   : top(nca.top), bottom(nca.bottom), candidate() {
}

Nucaln& Nucaln::operator=(const Nucaln &nca) {
   if (&nca != this) {
      top= nca.top;
      bottom= nca.bottom;
      candidate = nca.candidate;
   }
   return *this;
}

