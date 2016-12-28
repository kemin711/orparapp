#include "textpile.h"

int main(int argc, char* argv[]) {
   string inputFile = "testpiechart_data.tab";
   string outFile = "testtextpile_out.png";
   TextPile textpile;
   textpile.readData(inputFile);
   textpile.draw(outFile);

   return 0;
}

