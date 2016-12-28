#ifndef TEXTPILE_H
#define TEXTPILE_H

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

/**
 * Simple text drawing routine to substitute a PDF
 * writer.
 * Header is shared with blastnrow.
 * Future versions should consider connect these two
 * classes.
 * Right now taking top 10 rows. This number can be
 * configured to any other number.
 */
class TextPile {
   public:
      TextPile() : data(), plotHeight(0), plotWidth(0),
         fontSize(16), viewport(0), maxstr(), 
         numtoshow(10)  { }
      TextPile(const vector<vector<string> > &input) 
         : data(input), plotHeight(0), plotWidth(0), fontSize(16), 
         viewport(0), maxstr(header), numtoshow(10) 
         { setViewport(); }
      ~TextPile() { if (viewport != 0) delete[] viewport; }

      /**
       * Method to draw a single table in png format
       * @param outfile output file name for the image file 
       */
      void draw(const string &outfile);
      void setInput(vector<vector<string> > &&input) { data = std::move(input); }
      void setViewport();

      /**
       * Helper and testing function used at development stage.
       * Will fill the om member.
       * Auto detect if number of fields is 3 then
       * assumes no header and it is directaln result.
       * Otherwise assumes otu mapping result 14, or 15 columns
       * and has header.
       */
      void readData(const string &file);
      int getNumtoshow() const { return numtoshow; }
      static void setHeader(const vector<string> &h) { header=h; }

   private:
      vector<vector<string> > data;
      /**
       * default height 1200
       */
      int plotHeight;
      /**
       * default width 1800
       */
      int plotWidth;
      /**
       * Default 18, this can be changed
       */
      int fontSize;
      char *viewport;
      vector<string> maxstr;
      int numtoshow;

      /**
       * Helper method
       */
      static void getRandomColor(int &R, int &G, int &B);
      static vector<string> header;
};

#endif
