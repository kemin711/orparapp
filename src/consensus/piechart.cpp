#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <plotter.h>
#include <iterator>
#include <algorithm>
// the following won't bring out PI_M
//#define _USE_MATH_DEFINES
#include <cmath> // does not seem to have PI_M

#include <strformat.h>
#include "piechart.h"

using namespace orpara;



PieChart::~PieChart() {
   if (viewport != 0) delete[] viewport;
}

void PieChart::initTree() {
   char *metagDataDir=getenv("METAG_DATA");
   string nodeFile = metagDataDir + string("/taxonomy/nodes.dmp");
   taxtree.buildTree(nodeFile);
   string nameFile = metagDataDir + string("/taxonomy/names.dmp");
   Taxon::loadTaxName(nameFile);
}

void PieChart::updateTaxonCount() {
   taxtree.countTaxon(om);
   genusCount = taxtree.getGenusCount();
   cout << "\ngenus hit_count\n";
   for (auto i=0; i<genusCount.size(); ++i) {
      cout << genusCount[i].first << " " << genusCount[i].second << endl;
   }
   familyCount = taxtree.getFamilyCount();
   cout << "\nfamily hit_count\n";
   for (auto i=0; i<familyCount.size(); ++i) {
      cout << familyCount[i].first << " " << familyCount[i].second << endl;
   }
#if __GNUC__ > 4
   sort(genusCount.begin(), genusCount.end(), 
         [](pair<string,int> &a, pair<string,int> &b)->bool { return a.second > b.second; });
   sort(familyCount.begin(), familyCount.end(), 
         [](pair<string,int> &a, pair<string,int> &b)->bool { return a.second > b.second; });
#else
   sort(genusCount.begin(), genusCount.end(), CompareSecondValue());
   sort(familyCount.begin(), familyCount.end(), CompareSecondValue());
#endif
}

// for simulating the original data type
// if same otu has two mappings 
//   if one has higher mapping the then use it
//   if both has identical mapping then use the latter (has scientific name)
void PieChart::readOtumap(const string &file) {
   om = taxtree.readMapping(file);
   cout << om.size() << " OTUs\n";
   // no need to sort OTU, need to sort Genus
   updateTaxonCount();
}

void PieChart::readDirectaln(const string &file) {
   ifstream inf(file);
   string line;
   getline(inf, line);
   int i = 1;
   while (!inf.eof()) {
      vector<string> row = split(line, "\t");
      om.push_back(make_tuple(to_string(i), stoi(row[2]), row[0]));
      getline(inf, line);
   }
   cout << om.size() << " OTUs\n";
   updateTaxonCount();
}

//pair<int,string> getTotalDepth(const vector<tuple<string,int,string,float> > &otus, size_t b=0) {
pair<int,string> PieChart::measureBound(const vector<pair<string,int> > &data) {
   int total=0;
   string maxstr;
   for (auto i=0; i<data.size(); ++i) {
      total += data[i].second;
      if (data[i].first.length() > maxstr.length()) {
         maxstr=data[i].first;
      }
   }
   cout << "longest taxon name: " << maxstr << endl;
   return pair<int, string>(total, maxstr);
}

void PieChart::getRandomColor(int &R, int &G, int &B) {
   static int colormax=65536;
   R=rand()%colormax;
   G=rand()%colormax;
   B=rand()%colormax;
}
   
void PieChart::draw() {
   setViewport();
   chart(genusCount, "Genus", "Genus.pie.svg");
   chart(familyCount, "Family", "Family.pie.svg");
   //chart(genusCount, "Genus", "Genus.pie.png");
   //chart(familyCount, "Family", "Family.pie.png");
}

void PieChart::setViewport() {
   float h = round(100*plotHeight/double(plotWidth))/10;
   string tmp;
   if (outfmt == SVG) {
      tmp = "a4,xsize=10in,ysize=" + to_string(h).substr(0,3) + "in";
   }
   else {
      tmp=to_string(plotWidth) + "x" + to_string(plotHeight);
   }
   cout << "viewport setting: " << tmp << endl;
   if (viewport != 0) delete[] viewport;
   viewport = new char[tmp.size()+1];
   strcpy(viewport, tmp.c_str());
}

pair<double, char> textAnglePolicy1(double angle, float &szscale) {
   if (angle>0 && angle <= M_PI/2) {
      szscale=3;
      return make_pair(0, 'l');
   }
   if (angle>M_PI/2 && angle <= 3*M_PI/2) {
      szscale=2;
      return make_pair(0, 'r');
   }
   if (angle < M_PI*11/6) {
      szscale = 1.5;
      return make_pair(330, 'l');
   }
   szscale = 1;
   return make_pair(angle*180/M_PI, 'l');
}

/**
 * Default start from 0. You can start from a none zero index.
 * @return the otu index in om when it is too small to see.
 */
void PieChart::chart(const vector<pair<string, int> > &taxcnt,
      const string &title, const string &outfile) 
{
   PlotterParams params;
   ofstream ouf(outfile);
   Plotter *plt;
   cout << "user space: " << plotWidth << " x " << plotHeight << endl;
   int drawing_w=int(double(plotWidth)*1.09);
   int drawing_h=int(double(plotHeight)*1.09);
   if (outfmt == SVG) {
      params.setplparam("PAGESIZE", viewport);
      plt = new SVGPlotter(cin, ouf, cerr, params);
   }
   else { // for PNG
      params.setplparam("BITMAPSIZE", viewport);
      plt = new PNGPlotter(cin, ouf, cerr, params);
   }
   if (plt->openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt->erase();
   plt->fspace(-(0.5*drawing_w),-(0.5*drawing_h), 0.5*drawing_w, 0.5*drawing_h);
   plt->flinewidth(2);
   plt->pencolorname("red");
   double xw=0.5*plotWidth;
   double yh=0.5*plotHeight;
   plt->fbox(-xw,-yh, xw, yh);
   pair<int,string> sd_maxw = measureBound(taxcnt);
   cout << "total depth: " << sd_maxw.first 
      << " longest taxon name: " << sd_maxw.second << "\n";
   int sumdepth=sd_maxw.first;
   plt->colorname("black");
   plt->flinewidth(1);
   plt->fontname("HersheySerif");
   int true_size = plt->fontsize(4*fontSize);
   plt->fmove(-xw + 2*true_size, yh - 2*true_size);
   plt->alabel('l', 't', title.c_str());
   true_size = plt->fontsize(fontSize);
   string tmp = sd_maxw.second + " " + to_string(sumdepth);
   int maxwidthLabel = plt->labelwidth(tmp.c_str());
   if (maxwidthLabel + radius > plotWidth || maxwidthLabel + radius > plotHeight) {
      cerr << "some labels will outside the box!\n";
   }
   // draw pies
   plt->flinewidth(0.50);
   plt->filltype(1); // use solid fill 0x8000 50% fill, higher the lighter
   double theta=0, dTheta=0; // in radian
   double x0, y0, x1, y1, x2, y2;
   int currD=0;
   int colR, colG, colB;
   plt->pencolorname("black");
   float txtscale=1;
   double blurY = -yh+true_size;
   for (auto i=0; i<taxcnt.size(); ++i) {
      getRandomColor(colR, colG, colB);
      //plt->color(colR, colG, colB);
      plt->fillcolor(colR, colG, colB);
      double dTheta = 2*M_PI*double(taxcnt[i].second)/sumdepth; // radian
      theta = 2*M_PI*(double)currD/sumdepth; // this get rid of rounding error
      //cout << "theta=" << theta << " delTheta=" << dTheta << endl;
      angle2xy(theta, x0, y0);
      angle2xy(theta+dTheta, x1, y1);
      plt->fline(x0, y0, 0, 0);
      // if dTheta > PI then you need to draw two pieces!
      if (dTheta > M_PI) {
         angle2xy(theta+M_PI, x2, y2);
         cerr << "draw two pices!\n";
         // draw the Pi part
         plt->farc(0, 0, x0, y0, x2, y2);
         plt->fmove(0,0);
         plt->farc(0, 0, x2, y2, x1, y1);
      }
      else { // one arc is enough
         plt->farc(0, 0, x0, y0, x1, y1);
      }
      plt->fline(x1, y1, 0, 0);
      plt->endpath(); // don't need to call this
      // draw label on the side
      string lab=taxcnt[i].first + " " + to_string(taxcnt[i].second);
      double halfAngle = theta+0.5*dTheta;
      plt->fmove((radius+2*true_size)*cos(halfAngle), (radius+2*true_size)*sin(halfAngle));
      pair<double, char> txtpos = textAnglePolicy1(halfAngle, txtscale);
      plt->ftextangle(txtpos.first); // in degree
      plt->fontsize(txtscale*fontSize);
      plt->alabel(txtpos.second, 'c', lab.c_str());
      if (2*radius*sin(0.5*dTheta) < true_size*0.8) {
         cerr << taxcnt[i].first << " abundance too low to be seen on graph\n";
         plt->fmove(xw-2*true_size, blurY);
         plt->ftextangle(0);
         plt->alabel('r', 'b', lab.c_str());
         blurY += true_size;
      }
      currD += taxcnt[i].second;
   }
   plt->flinewidth(2);
   plt->pencolorname("black");
   plt->filltype(0);
   plt->fcircle(0,0, radius);
   plt->closepl();
   cerr << "pie chart written to " << outfile << endl;
   delete plt;
}


