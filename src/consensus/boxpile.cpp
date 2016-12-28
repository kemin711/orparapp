#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <plotter.h>
#include <iterator>
#include <algorithm>
#include <stdexcept>

#include <strformat.h>
#include "boxpile.h"

using namespace orpara;

//vector<tuple<string, int, string, float> > readOtumap(const string &file);
//void buildBar(const vector<tuple<string,int,string, float> > &om);
//int stackBox(Plotter &plt, const vector<tuple<string,int,string, float> > &om, size_t b,
//      int plotH, int leftX, const string &bColor, int totalRead);

// for simulating the original data type
// if same otu has two mappings 
//   if one has higher mapping the then use it
//   if both has identical mapping then use the latter (has scientific name)
void BoxPile::readOtumap(const string &file) {
   ifstream ifs(file);
   string line, clname; 
   getline(ifs, line); // header
   getline(ifs, line); // first data line
   // clid, depth, species, identity
   om.clear();
   //vector<tuple<string, int, string, float> > res;
   //vector<string> lastrow;
   while (!ifs.eof()) {
      //cout << line << endl;
      vector<string> row = split(line, "\t");
      /*
      if (!lastrow.empty() && row[0] == lastrow[0]) {
         if (stoi(row[2]) >= stoi(lastrow[2])) { // replace last one with this one
            om.back() = make_tuple(row[0], stoi(row[2]), row[4], stof(row[6]));
         }
         else { // saving the last one
            cout << "multiple mapping for cluster " << row[0] << endl
               << "ignoring the weaker mapping\n" << line << endl
               << "saving\n";
            copy(lastrow.begin(), lastrow.end(), ostream_iterator<string>(cout, " | "));
            cout << endl;
         }
      }
      else {
      */
      om.push_back(make_tuple(row[0], stoi(row[2]), row[4], stof(row[6])));
      getline(ifs, line);
      //lastrow = row;
   }
   cout << om.size() << " OTUs\n";
   removeDuplicatedData();
   sortData();
   //cout << "finished sorting\n";
   cleanom();
   //cout << "last element: " << get<2>(om.back()) << endl;
}

// file has no header
// Escherichia coli  0.9967304933067939   9699
// TAB-delimited row: species identity count
void BoxPile::readDirectaln(const string &file) {
   ifstream ifs(file);
   string line, clname; 
   getline(ifs, line); // first data line
   // species, identity, count
   om.clear();
   int i=1;
   while (!ifs.eof()) {
      //cout << line << endl;
      vector<string> row = split(line, "\t");
      //tuple content: clusterid, depth (or count), taxon, mapping identity
      om.push_back(make_tuple(to_string(i++), stoi(row[2]), row[0], stof(row[1])));
      getline(ifs, line);
   }
   cout << om.size() << " OTUs or Directaln targets\n";
   sortData(); // input is already sorted, just to be safe
}

int BoxPile::numberOfFields(const string &file) {
   ifstream inf(file);
   string line;
   getline(inf, line);
   int nf = 0;
   for (auto i=0; i<line.size(); ++i) 
      if (line[i] == '\t') ++nf;
   return nf+1;
}

void BoxPile::removeDuplicatedData() {
   vector<tuple<string,int,string,float> > res;
   res.push_back(om.front());
   size_t i = 1;
   tuple<string,int,string,float> last=om.front();
   int dep, lastdep;
   while (i < om.size()) {
      dep = get<1>(om[i]);
      lastdep = get<1>(last);
      if (get<0>(om[i]) == get<0>(last)) {
         if (get<1>(om[i]) >= get<1>(last)) { // replace last one with this one
            res.back() = om[i]; 
         }
         else { // saving the last one
            cout << "multiple mapping for cluster " << get<0>(last) << endl
               << "ignoring the weaker mapping\n" << get<1>(om[i]) << " " << get<2>(om[i]) << endl
               << "saving\n" << get<1>(last) << " " << get<2>(last) << endl;
         }
      }
      else {
         res.push_back(om[i]);
      }
      last = om[i];
      ++i;
   }
   if (res.size() < om.size()) om=res;
}

// 4.8. compiler does not handle lambda well!
// I have to use function object
void BoxPile::sortData() {
#if __GNUC__ > 4
   sort(om.begin(), om.end(), [](tuple<string,int,string,float> &a, tuple<string,int,string,float> &b)->bool { return get<1>(a) > get<1>(b); });
#else
   sort(om.begin(), om.end(), SortOtuByDepth());
#endif
}

//pair<int,string> getTotalDepth(const vector<tuple<string,int,string,float> > &otus, size_t b=0) {
pair<int,string> BoxPile::depthFrom(size_t b) {
   int total=0;
   string maxorg;
   for (auto i=b; i<om.size(); ++i) {
      total += get<1>(om[i]);
      if (get<2>(om[i]).length() > maxorg.length()) {
         maxorg=get<2>(om[i]);
      }
   }
   cout << "longest organism name: " << maxorg << endl;
   return pair<int, string>(total, maxorg);
}

void BoxPile::getRandomColor(int &R, int &G, int &B) {
   static int colormax=65536;
   R=rand()%colormax;
   G=rand()%colormax;
   B=rand()%colormax;
}

void BoxPile::setViewport() {
   float h = round(100*plotHeight/double(plotWidth))/10;
   string tmp = "a4,xsize=10in,ysize=" + to_string(h).substr(0,3) + "in";
   cout << "viewport setting: " << tmp << endl;
   if (viewport != 0) delete[] viewport;
   viewport = new char[tmp.size()+1];
   strcpy(viewport, tmp.c_str());
}

   
void BoxPile::draw(const string &outfile) {
   pair<int,string> sd_maxw = depthFrom();
   cout << "total depth: " << sd_maxw.first << " longest organism name: " 
      << sd_maxw.second << "\n";
   totalDepth = sd_maxw.first;
   // set up parameter
   PlotterParams params;
   //params.setplparam("PAGESIZE", (char*)"letter");
   params.setplparam("PAGESIZE", viewport);
   int drawing_w=int(double(plotWidth)*1.09);
   int drawing_h=int(double(plotHeight)*1.04);
   ofstream ouf(outfile);
   SVGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.erase();
   plt.fspace(-(0.1*plotWidth),-(0.1*plotHeight), drawing_w, drawing_h);
   plt.flinewidth(4);
   plt.pencolorname("black");
   // draw a border for the plot
   plt.fbox(-4, -4, plotWidth*1.04, plotHeight*1.02);
   cerr << "drawing width: " << drawing_w << " plot width: " << plotWidth << endl;
   plt.flinewidth(1);
   // use plotfont -T png --help-fonts to get a list
   //plt.fontname("HersheySans");
   plt.fontname("HersheySerif");
   int true_size = plt.fontsize(fontSize);
   string tmp = "999:" + sd_maxw.second + "|999|0.099|99.99";
   //cerr << "tmp " << tmp << endl;
   int maxwidth=plt.labelwidth(tmp.c_str());
   // draw boxes
   vector<string> borderColor={"red", "green", "blue"};
   size_t x=0;
   int toosmalli=0;
   while (x<3 && toosmalli > -1) {
      toosmalli=stackBox(plt, toosmalli, x*(1.02*maxwidth+30), borderColor[x]);
      ++x;
   }
   plt.closepl();
   cerr << "graphics written to " << outfile << endl;
}

void BoxPile::draw2(const string &outfile) {
   pair<int,string> sd_maxw = depthFrom();
   cout << "total depth: " << sd_maxw.first << " longest organism name: " 
      << sd_maxw.second << "\n";
   totalDepth = sd_maxw.first;
   // set up parameter
   PlotterParams params;
   //params.setplparam("PAGESIZE", (char*)"letter");
   params.setplparam("PAGESIZE", viewport);
   int drawing_w=int(double(plotWidth)*1.09);
   int drawing_h=int(double(plotHeight)*1.04);
   ofstream ouf(outfile);
   SVGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.erase();
   plt.fspace(-(0.1*plotWidth),-(0.1*plotHeight), drawing_w, drawing_h);
   plt.flinewidth(4);
   plt.pencolorname("black");
   // draw a border for the plot
   plt.fbox(-4, -4, plotWidth*1.04, plotHeight*1.02);
   cerr << "drawing width: " << drawing_w << " plot width: " << plotWidth << endl;
   plt.flinewidth(1);
   // use plotfont -T png --help-fonts to get a list
   //plt.fontname("HersheySans");
   plt.fontname("HersheySerif");
   int true_size = plt.fontsize(fontSize*1.6);
   string tmp = sd_maxw.second + " | 99.99";
   int maxwidth=plt.labelwidth(tmp.c_str());
   cerr << "maxwidth: " << maxwidth << " tmp " << tmp << endl;
   // draw box only and tex on the size
   // draw boxes
   plt.filltype(1); //1 for solid fill 0x8000 50% fill, higher the lighter
   size_t i;
   int currH=0; // height of the left bar
   int colR, colG, colB;
   int boxLeft=10;
   int legendLeft=250;
   // the y-axis label
   plt.fontsize(fontSize*2);
   plt.pencolorname("black");
   for (int p=0; p<=100; p+=10) {
      float y = (double)p*plotHeight/100;
      plt.fline(70, y, -3, y);
      plt.fmove(-9, y);
      plt.alabel('r', 'c', (to_string(p) + "%").c_str());
   }
   true_size = plt.fontsize(fontSize*1.6);
   plt.endpath();
   int rightShift=0;
   plt.flinewidth(0.25);
   double Hleg=0;
   for (i=0; i<om.size(); ++i) {
      getRandomColor(colR, colG, colB);
      plt.color(colR, colG, colB); // both pen and fill are the same
      double dH = get<1>(om[i])*plotHeight/(double)totalDepth;
      double H = (double)currH/totalDepth*plotHeight;
      plt.fbox(boxLeft, H, 60, H+dH);
      //cerr << "rightShift: " << rightShift << endl;
      // the legend box hight adjustment
      cerr << "Hleg: " << Hleg << " rightShift: " << rightShift << endl;
      if (Hleg+true_size > plotHeight) {
         Hleg = 0;
         ++rightShift;
         cerr << "reset Hleg: " << Hleg << "!!\n";
         cerr << "rightShift: " << rightShift << endl;
      }
      double shiftWidth = rightShift*maxwidth;
      double leftX=shiftWidth + legendLeft;
      if (rightShift > 0) leftX += 60;
      // draw legend box
      if (rightShift <= 1)
         plt.fbox(leftX, Hleg, leftX+25, Hleg+true_size);
      // draw legen text
      plt.fmove(leftX+40, Hleg);
      plt.pencolorname("black");
      string lab=get<2>(om[i]) + " | "
         + to_string(round(get<3>(om[i])*100)/100).substr(0,5);
      cerr << "label: " << lab << " i=" << i << " of total=" << om.size() << endl;
      if (rightShift > 1) {
         cerr << "rightShift=" << rightShift << " "
            << i << " labels reached limit, ingnoring the rest" << endl;
      }
      else {
         plt.alabel('l', 'b', lab.c_str());      
         //cerr << "label: " << lab << endl;
      }
      if (rightShift < 1) { // connecting line
         plt.fline(60, H+0.5*dH, legendLeft+1, (i+0.5)*true_size);
      }
      currH += get<1>(om[i]);
      Hleg += true_size;
   }
   plt.endpath();
   plt.closepl();
   cerr << "graphics written to " << outfile << endl;
}

// combine >= 1% together
void BoxPile::draw3(const string &outfile) {
   cout << "BoxPile::draw3() drawing the simplest version of the bar graph ...\n";
   pair<int,string> sd_maxw = depthFrom();
   cout << "total depth: " << sd_maxw.first << " longest organism name: " 
      << sd_maxw.second << "\n";
   totalDepth = sd_maxw.first;
   static double cutoff = 0.01; // TODO: should be an parameter 
   // combine low
   vector<tuple<string,int,string,float> > combined;
   size_t i=0;
   while (i < om.size() && (double)get<1>(om[i])/totalDepth >= cutoff) {
      //cout << get<2>(om[i]) << " " << (double)get<1>(om[i])/totalDepth << endl;
      combined.push_back(om[i]);
      ++i;
   }
   cout << combined.size() << " OTU has >= " << int(cutoff*100) << "% abundance\n";
   int lowCnt=0;
   cout << "low abundance taxons:\n";
   while (i<om.size()) {
      cout << get<2>(om[i]) << " " << get<1>(om[i]) << " " << (double)get<1>(om[i])/totalDepth << endl;
      lowCnt += get<1>(om[i]);
      ++i;
   }
   combined.push_back(make_tuple(string("cl<1% combined"), lowCnt, string("Others: Taxa < 1%"), (float)lowCnt/totalDepth)); 
   cout << lowCnt << " total <1% counts\n";
   // set up parameter
   PlotterParams params;
   //params.setplparam("PAGESIZE", (char*)"letter");
   params.setplparam("PAGESIZE", viewport);
   int drawing_w=int(double(plotWidth)*1.09);
   int drawing_h=int(double(plotHeight)*1.04);
   ofstream ouf(outfile);
   SVGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.erase();
   plt.fspace(-(0.1*plotWidth),-(0.1*plotHeight), drawing_w, drawing_h);
   plt.pencolorname("black");
   // draw a border for the plot
   plt.flinewidth(4);
   plt.fline(0, 0, plotWidth, 0);
   plt.fline(0, 0, 0, plotHeight*1.02);
   cerr << "drawing width: " << drawing_w << " plot width: " << plotWidth << endl;
   plt.flinewidth(1);
   // use plotfont -T png --help-fonts to get a list
   //plt.fontname("HersheySans");
   plt.fontname("HersheySerif");
   int true_size = plt.fontsize(fontSize*1.6);
   string tmp = sd_maxw.second + " | 99.99";
   int maxwidth=plt.labelwidth(tmp.c_str());
   cerr << "maxwidth: " << maxwidth << " tmp " << tmp << endl;
   // draw box only and tex on the size
   // draw boxes
   plt.filltype(1); //1 for solid fill 0x8000 50% fill, higher the lighter
   int currH=0;
   int colR, colG, colB;
   int boxLeft=10;
   int legendLeft=250;
   // the y-axis label
   plt.fontsize(fontSize*2);
   plt.pencolorname("black");
   for (int p=0; p<=100; p+=10) {
      float y = (double)p*plotHeight/100;
      plt.fline(70, y, -3, y);
      plt.fmove(-9, y);
      plt.alabel('r', 'c', (to_string(p) + "%").c_str());
   }
   true_size = plt.fontsize(fontSize*1.6);
   plt.endpath();
   double Hleg=0;
   for (i=0; i<combined.size(); ++i) {
      getRandomColor(colR, colG, colB);
      plt.color(colR, colG, colB); // both pen and fill are the same
      double dH = get<1>(combined[i])*plotHeight/(double)totalDepth;
      double H = (double)currH/totalDepth*plotHeight;
      plt.fbox(boxLeft, H, 60, H+dH);
      // the legend box
      plt.fbox(legendLeft, Hleg, legendLeft+25, Hleg+true_size);
      plt.fmove(legendLeft+40, Hleg);
      plt.pencolorname("black");
      string lab=get<2>(combined[i]);
      plt.alabel('l', 'b', lab.c_str());      
      currH += get<1>(combined[i]);
      Hleg += true_size;
   }
   plt.endpath();
   plt.closepl();
   cerr << "Simplified graphics written to " << outfile << endl;
}

string extractClusterId(const string &clname) {
   string::size_type i = 0;
   string num;
   while (i < clname.size()) {
      if (isdigit(clname[i])) {
         while (i<clname.size() && isdigit(clname[i])) {
            num += clname[i];
            ++i;
         }
         return num;
      }
      ++i;
   }
   throw runtime_error(clname + " does not contain digit");
}



/**
 * Default start from 0. You can start from a none zero index.
 * @return the otu index in om when it is too small to see.
 */
int BoxPile::stackBox(Plotter &plt, size_t b, int leftX, const string &bColor)
{
   pair<int,string> sd_maxw = depthFrom(b);
   cout << "starting from index " << b << " total depth: " << sd_maxw.first 
      << " longest organism name: " << sd_maxw.second << "\n";
   int denominator = 1000;
   int numdigits=5;
   if ((double)totalDepth/get<1>(om[b]) > 1000) {
      denominator=10000;
      numdigits=6;
   }
   int sumdepth=sd_maxw.first;
   // add more boxes
   plt.flinewidth(1);
   plt.fontname("HersheySerif");
   int true_size = plt.fontsize(fontSize);
   string tmp = "999:" + sd_maxw.second + "|999|0.099|99.99";
   int maxwidthCurrent=plt.labelwidth(tmp.c_str());
   // draw boxes
   plt.flinewidth(0.25);
   plt.filltype(0xDD00); // 0x8000 50% fill, higher the lighter
   size_t i;
   int lowi=-1;
   int currH=0;
   int colR, colG, colB;
   for (i=b; i<om.size(); ++i) {
      getRandomColor(colR, colG, colB);
      plt.fillcolor(colR, colG, colB);
      double dH = get<1>(om[i])*plotHeight/(double)sumdepth;
      double H = (double)currH/sumdepth*plotHeight;
      plt.pencolorname(bColor.c_str());
      plt.fbox(leftX, H, leftX+1.02*maxwidthCurrent, H+dH);
      plt.pencolorname("black");
      string lab=extractClusterId(get<0>(om[i])) + ":"
         + get<2>(om[i]) + "|" + to_string(get<1>(om[i])) + "|"
         + to_string(round(get<1>(om[i])/(double)totalDepth*denominator)/denominator).substr(0,numdigits) + "|"
         + to_string(round(get<3>(om[i])*100)/100).substr(0,5);
      plt.alabel('c', 'c', lab.c_str());
      if (dH < true_size*0.8 && lowi == -1) {
         cerr << get<0>(om[i]) << " abundance too low to be seen on graph\n";
         lowi=i;
      }
      currH += get<1>(om[i]);
   }
   cout << lowi << " index for low abundance index\n";
   plt.endpath();
   return lowi;
}

// remove square bracket
void BoxPile::cleanom() {
   removeSquare();
   removeDuplicatedData();
   sortData();
}

void BoxPile::removeSquare() {
   for (auto i=0; i<om.size(); ++i) {
      string org=get<2>(om[i]);
      string::size_type x = org.find(']');
      bool findsquare=false;
      if (x != string::npos) {
         org.erase(x, 1);
         findsquare=true;
      }
      x = org.find('[');
      if (x != string::npos) {
         org.erase(x, 1);
         findsquare=true;
      }
      if (findsquare) get<2>(om[i])=org;
   }
}

