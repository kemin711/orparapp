#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <plotter.h>
#include <iterator>
#include <algorithm>

#include <strformat.h>
#include "textpile.h"

using namespace orpara;


// we don't know the header before hand
vector<string> TextPile::header={};

//vector<tuple<string, int, string, float> > readOtumap(const string &file);
//void buildBar(const vector<tuple<string,int,string, float> > &om);
//int stackBox(Plotter &plt, const vector<tuple<string,int,string, float> > &om, size_t b,
//      int plotH, int leftX, const string &bColor, int totalRead);

// for simulating the original data type
// if same otu has two mappings 
//   if one has higher mapping the then use it
//   if both has identical mapping then use the latter (has scientific name)
void TextPile::readData(const string &file) {
   ifstream ifs(file);
   string line, clname; 
   vector<string> row, lastrow;
   getline(ifs, line); // header
   row = split(line, "\t");
   if (row.size() == 3) {
      header = {"Species", "AVG_identity", "numreads"};
      cerr << "Use default Directaln header: ";
      copy(header.begin(), header.end(), ostream_iterator<string>(cerr, " "));
      cerr << endl;
   }
   else {
      header=row;
      if (header[header.size()-1].empty()) {
         header.resize(header.size()-1);
      }
      cerr << header.size() << " header fields\n";
      getline(ifs, line); // first data line
   }

   data.clear();
   while (!ifs.eof()) {
      row = split(line, "\t");
      if (row.size() != header.size()) {
         cerr << header.size() << " header:\n";
         copy(header.begin(), header.end(), ostream_iterator<string>(cerr, " | "));
         cerr << endl;
         cerr << row.size() << " data:\n";
         copy(row.begin(), row.end(), ostream_iterator<string>(cerr, " | "));
         throw runtime_error("data not the same size as header");
      }
      if (!lastrow.empty() && row[0] == lastrow[0]) {
         if (stoi(row[2]) >= stoi(lastrow[2])) { // replace last one with this one
            data.back() = row;
         }
         else { // saving the last one
            cout << "multiple mapping for cluster " << row[0] << endl
                 << "ignoring the weaker mapping\n" << line << endl
                 << "saving\n";
            copy(lastrow.begin(), lastrow.end(), ostream_iterator<string>(cout, " | "));
            cout << endl;
         }
      }
      else { data.push_back(row); }
      getline(ifs, line);
      lastrow = row;
   }
   sort(data.begin(), data.end(), [](const vector<string> &a, const vector<string> &b)->bool { return stoi(a[2]) > stoi(b[2]); });
   cout << data.size() << " OTUs collected from " << file << "\n";
   setViewport();
}

void TextPile::getRandomColor(int &R, int &G, int &B) {
   static int colormax=65536;
   R=rand()%colormax;
   G=rand()%colormax;
   B=rand()%colormax;
}

// should measure the column width of each data/header
// get max width
// then set up the maximum width of each field Plus
// spacer for each one of two spaces.
// use png format
void TextPile::setViewport() {
   cerr << numtoshow << " rows to show\n";
   if (maxstr.size() < header.size()) {
      maxstr = header;
   }
   //for (auto i=0; i<data.size(); ++i) {
   for (auto i=0; i<numtoshow && i<data.size(); ++i) {
      for (auto j=0; j<data[i].size(); ++j) {
         if (data[i][j].length() > maxstr[j].length())
            maxstr[j] = data[i][j];
      }
   }
   cout << "maximum length of each field\n";
   copy(maxstr.begin(), maxstr.end(), ostream_iterator<string>(cout, ", "));
   //plotHeight = (data.size()+1)*fontSize*1.2;
   plotHeight = (numtoshow+3)*fontSize*1.2;
   // add 3 spaces for padding
   //plotWidth = accumulate(maxstr.begin(), maxstr.end(), maxstr.size()*2,
   plotWidth = accumulate(maxstr.begin(), maxstr.end(), 0,
         [](size_t a, const string &b)->int { return a + b.length(); } )*fontSize*0.8;
   string tmp = to_string(plotWidth) + "x" + to_string(plotHeight);
   cout << "viewport setting: " << tmp << endl;
   if (viewport != 0) delete[] viewport;
   viewport = new char[tmp.size()+1];
   strcpy(viewport, tmp.c_str());
}
   
void TextPile::draw(const string &outfile) {
   PlotterParams params;
   // use png format
   params.setplparam("BITMAPSIZE", viewport);
   int drawing_w=int(double(plotWidth)*1.09);
   int drawing_h=int(double(plotHeight)*1.04);
   ofstream ouf(outfile);
   PNGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.erase();
   plt.fspace(-(0.1*plotWidth),-(0.1*plotHeight), drawing_w, drawing_h);
   plt.pencolorname("black");
   // draw a border for the plot
   //plt.fbox(-4, -4, plotWidth*1.04, plotHeight*1.02);
   cerr << "drawing width: " << drawing_w << " plot width: " << plotWidth << endl;
   plt.flinewidth(1);
   // use plotfont -T png --help-fonts to get a list
   plt.fontname("HersheySans");
   //plt.fontname("HersheySerif");
   //plt.fontname("NewCenturySchlbk-Roman");
   float true_size = plt.fontsize(fontSize);
   vector<float> maxlen(maxstr.size());
   for (auto i=0; i<maxstr.size(); ++i) {
      maxlen[i] = plt.flabelwidth(maxstr[i].c_str());
   }
   double widthOfText = accumulate(maxlen.begin(), maxlen.end(), 0);
   widthOfText += maxlen.size()*2*true_size;
   int colR, colG, colB;
   double H = 0.5*true_size;
   float left;
   //for (int i=data.size()-1; i>-1; --i) {
   int total = numtoshow;
   if (data.size() < numtoshow) 
      total = data.size();
   plt.filltype(0xDD00);
   plt.pencolorname("black");
   for (int i=total-1; i>-1; --i) {
      getRandomColor(colR, colG, colB);
      plt.fillcolor(colR, colG, colB);
      //plt.color(colR, colG, colB);
      plt.fbox(-true_size, H-0.5*true_size*1.2, widthOfText, H+0.5*true_size*1.2);
      left=0;
      for (auto f=0; f<header.size(); ++f) {
         plt.fmove(left, H);
         plt.alabel('l', 'c', data[i][f].c_str());
         left += (maxlen[f] + 2*true_size);
      }
      H += true_size*1.2; // moving up one line at a time
   }
   plt.pencolorname("black");
   plt.fline(-true_size, H, left, H);
   left=0;
   H += true_size*1.2;
   plt.fontname("HersheySerif-Bold");
   for (auto f=0; f<header.size(); ++f) {
      plt.fmove(left, H);
      plt.alabel('l', 'c', header[f].c_str());
      left += (maxlen[f] + 2*true_size);
   }
   plt.flinewidth(4);
   H += true_size*1.2;
   plt.filltype(0);
   plt.fbox(-true_size, -4, left, H);
   plt.closepl();
   cerr << "Tabular text written to " << outfile << endl;
}

