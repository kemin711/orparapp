#ifndef BOXPILE_H
#define BOXPILE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <plotter.h>
#include <iterator>
#include <algorithm>

#include <strformat.h>

using namespace std;
using namespace orpara;

class SortOtuByDepth {
   public:
      bool operator()(const tuple<string, int, string, float> &a, const tuple<string,int,string,float> &b) {
         return get<1>(a) > get<1>(b); 
      }
};

/**
 * Given the OTU count and mapping identity, this class
 * is able to produce three bar charts:
 * 1. Detailed box pile
 * 2. Stacked bar of proportonal hight
 * 3. Stacked bar same as 2 but will combin all boxes < 1% into one
 *    box
 * There is no need for repeated action. So a single object
 * will just do one thing. Right now SVG output format.
 * This class can be easily modified to output in multiple
 * formats.
 */
class BoxPile {
   public:
      BoxPile() : om(), plotHeight(1300), plotWidth(2000), totalDepth(0),
         fontSize(16), viewport(0) { setViewport(); }
      /** 
       * @param input a vector of tuple of
       *   clusterid a string identifier for the cluster 
       *   depth the vertical depth of the consensus or count
       *   Species name
       *   mapping identity to reference sequence
       */
      BoxPile(const vector<tuple<string, int, string, float> > &input) : om(input), plotHeight(1200), plotWidth(2000),
         totalDepth(0), fontSize(16), viewport(0) { cleanom(); setViewport(); }
      ~BoxPile() { if (viewport != 0) delete[] viewport; }

      /**
       * Method to draw a single bar graph
       * The output file name is fixed for now: otu_bargraph.svg
       * @param outfile output file name for the bar graph
       */
      void draw(const string &outfile);
      /**
       * Put figure legends on the right-hand side, and lines connecting 
       * box to the first column of legeds. Use maximum two columns of legends.
       */
      void draw2(const string &outfile);
      /*
       * combine <1% into one category
       */
      void draw3(const string &outfile);
      void setInput(vector<tuple<string, int, string, float> > &&input) { om = std::move(input); }
      /**
       * Add one OTU mapping result
       * The input is from blastnrow getReduced() function.
       */
      void add(const tuple<string,int,string,float> otu) { om.push_back(otu); }
      void setViewport();

      /**
       * Helper and testing function used at development stage.
       * Will fill the om member.
       * This is also used for command line program input.
       * This function also removes duplicated entries.
       */
      void readOtumap(const string &file);
      /**
       * Read from the directaln output
       * that has refid, identity, count
       */
      void readDirectaln(const string &file);
      int numberOfFields(const string &file);
      bool isOtumap(const string &file) {
         return numberOfFields(file) == 15; }
      bool isDirectaln(const string &file) {
         return numberOfFields(file) == 3; }

   private:
      void cleanom();
      /**
       * Remove the square bracket for certain organism names
       */
      void removeSquare();
      void removeDuplicatedData();
      void sortData();
      int stackBox(Plotter &plt, size_t b, int leftX, const string &bColor);
      /**
       * Add up counts from a particular index b.
       */
      pair<int,string> depthFrom(size_t b=0);
      /* OTU mapping to taxonomy abbreviated result
       * clusterid, depth (or count), taxon, mapping identity
      */
      vector<tuple<string,int,string, float> > om;
      /**
       * default height 1200
       */
      int plotHeight;
      /**
       * default width 1800
       */
      int plotWidth;
      /**
       * Sum of all the otu's depth or count (directaln)
       * Or sum of om<1>
       */
      int totalDepth;
      /**
       * Default 18, this can be changed
       */
      int fontSize;
      char *viewport;

      /**
       * Helper method
       */
      static void getRandomColor(int &R, int &G, int &B);
};

#endif
