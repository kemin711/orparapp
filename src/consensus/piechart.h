#ifndef PIECHART_H
#define PIECHART_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <algorithm>
#include <plotter.h>
#include <iterator>

#include <strformat.h>
#include <ncbitaxonomy.h>

using namespace std;
using namespace orpara;

enum OutputFormat { SVG, PNG };

class CompareSecondValue {
   public:
      bool operator()(const pair<string,int> &a, const pair<string,int> &b) {
         return a.second > b.second; 
      }
};

/**
 * There is no need for repeated action. So a single object
 * will just do one thing. Right now SVG output format.
 * This class can be easily modified to output in multiple
 * formats.
 * TODO: PieChart and BoxPile should be dreived from
 * a common base class.  This is left for future release.
 * There are enough similarity between the two.
 */
class PieChart {
   public:
      PieChart() : om(), plotHeight(1200), plotWidth(1800), radius(300),
         fontSize(18), taxtree(), genusCount(), familyCount(),
         viewport(0), outfmt(SVG) { initTree(); }
      PieChart(const vector<tuple<string, int, string> > &input) : om(input), plotHeight(1200), plotWidth(1800),
         radius(300), fontSize(18), taxtree(), genusCount(), familyCount(),
         viewport(0), outfmt(SVG) { initTree(); updateTaxonCount(); }

      ~PieChart();

      /**
       * Method to draw a single bar graph
       * The output file name is fixed for now: otu_bargraph.svg
       * @param outfile output file name for the bar graph
       */
      //void draw(const string &outfile);
      void draw();
      void chart(const vector<pair<string, int> > &taxcnt, const string &title, const string &outfile);

      void setInput(vector<tuple<string, int, string> > &&input) { om = std::move(input); updateTaxonCount(); }

      /**
       * Add one OTU mapping result
       * cluster_id, depth, species + strain
       */
      void add(const tuple<string,int,string> &otu) { om.push_back(otu); }

      /**
       * Helper and testing function used at development stage.
       * Will fill the om member.
       * The input can be set either from the constructore or the 
       * setter function.
       * Command line tools will also need this function.
       */
      void readOtumap(const string &file);
      /**
       * Get input from directaln result
       */
      void readDirectaln(const string &file);
      void initTree();
      /**
       * Given the angle, set the x,y coordinate
       */
      void angle2xy(double radian, double &x, double &y) { x=radius*cos(radian); y=radius*sin(radian); }
      void setViewport();

   private:
      /**
       * Fill up genus and family hit count
       */
      void updateTaxonCount();
      /**
       * Add up counts from a particular index b.
       */
      //pair<int,string> depthFrom(size_t b=0);
      static pair<int,string> measureBound(const vector<pair<string,int> > &data);
      /**
       * OTU mapping to taxonomy abbreviated result
       * clusterid, depth, species, strain
       * Clusterid is not used in the operation more for 
       * standardization.
      */
      vector<tuple<string,int,string> > om;
      /**
       * default height 1200
       */
      int plotHeight;
      /**
       * default width 1800
       */
      int plotWidth;
      int radius;
      /**
       * Default 18, this can be changed
       */
      int fontSize;
      /**
       * Business logic for counting is stored in this class
       * Also used the input method of taxtree to read the
       * out mapping file.
       */
      NCBITaxonomy taxtree;
      vector<pair<string, int> > genusCount;
      vector<pair<string, int> > familyCount;
      char *viewport;
      OutputFormat outfmt;

      /**
       * Helper method
       */
      static void getRandomColor(int &R, int &G, int &B);
};

#endif
