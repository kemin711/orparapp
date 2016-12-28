#include <string>
#include <api/BamAlignment.h>
#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/SamHeader.h>
#include <cmath>
#include <functional>

using namespace std;

void countAccuracy(const string &bamfile);
void filterByAccuracy(const string &bamfile, function<bool(float)> &filter, const string &ofile);
void usage() {
   cerr << "bamfilter --count-accuracy bamfile\n"
      << "  The above command will count the accuracy in the file\n"
      << "bamfilter -a '<0.991' bamfile\n"
      << "  The above command will save reads with predicted accuracy\n"
      << "  less than 0.991\n"
      << "Options:\n"
      << "  -c or --count-accuracy FLAG count the accuracy of the\n"
      << "     input bam file\n"
      << "  -i STRING input bam file that can also be given as the\n"
      << "     positional parameter without any option\n"
      << "  -a '[<>=]FLOAT' indicate the cutoff direction\n"
      << "     and value.  For example, -a '<0.98' will only\n"
      << "     save read with accuracy less than 0.98\n"
      << "     default is >0.95\n"
      << "     Note: the option value must be quoted to prevent\n"
      << "           the shell from interpreting it differently\n"
      << "  -o STRING output file name. Bam file to be saved to\n"
      << "  -h or --help print this message\n";
   exit(1);
}
string makeOutFile(const string &infile, char cutd, float cutv);

int main(int argc, char* argv[]) {
   //string bamFile="A_ccs2_1.bam";
   string bamFile, outFile;
   function<bool(float)> filterfn = bind2nd(greater<float>(), 0.95);
   char cutd = '>'; // cutoff direction
   float cutv = 0.95; // cutoff value
   bool docount=false;
   int i=1;
   while (i < argc) {
      if (!strcmp(argv[i], "-i")) bamFile=argv[++i];
      else if (!strcmp(argv[i], "-c")) docount=true;
      else if (!strcmp(argv[i], "-a")) {
         string accuracyCut = argv[++i];
         cutd=accuracyCut[0];
         cutv=stof(accuracyCut.substr(1));
         if (cutd == '<') {
            filterfn = bind2nd(less<float>(), cutv);
         }
         else if (cutd == '>') {
            filterfn = bind2nd(greater<float>(), cutv);
         }
         else if (cutd == '=') {
            filterfn = bind2nd(equal_to<float>(), cutv);
         }
         else {
            cerr << "wrong format to specify accuracy cutoff\n";
         }
      }
      else if (!strcmp(argv[i], "--help") || !strcmp(argv[i], "-h")) usage();
      else if (!strcmp(argv[i], "-o")) outFile=argv[++i];
      else {
         bamFile = argv[i];
      }
      ++i;
   }
   if (bamFile.empty()) usage();

   if (docount) countAccuracy(bamFile);
   else {
      if (outFile.empty()) {
         outFile = makeOutFile(bamFile, cutd, cutv);
      }
      filterByAccuracy(bamFile, filterfn, outFile);
   }
   return 0;
}

string makeOutFile(const string &infile, char cutd, float cutv) {
   string outname=infile;
   string::size_type i = infile.rfind('.');
   if (i != string::npos) {
      outname = infile.substr(0, i);
   }
   string tag="greater";
   if (cutd == '<') tag = "less";
   if (cutd == '=') tag = "equalto";
   tag += to_string(int(round(cutv*1000)));
   outname += "_" + tag + ".bam";
   return outname;
}

void filterByAccuracy(const string &bamfile, function<bool(float)> &filter, const string &ofile) {
   cerr << "Filtering bam files ...\n";
   BamTools::BamReader breader;
   if (!breader.Open(bamfile)) {
      throw runtime_error("Failed to open bam file: " + bamfile);
   }
   /// add some header info. proto type
   BamTools::SamHeader sh = breader.GetHeader();
   sh.SortOrder="queryname";
   BamTools::SamProgram thisProgram("bamfilter-1.0");
   thisProgram.Name="bamfilter";
   thisProgram.Version="1.0";
   thisProgram.CommandLine="filter -a [<=>]cutoff bamfile need to do more coding here";
   sh.Programs.Add(thisProgram);
   //////////////////////////////////
   BamTools::BamAlignment baln;
   BamTools::RefVector noref;
   BamTools::BamWriter bwriter;
   bwriter.Open(ofile, sh, noref);
   while (breader.GetNextAlignment(baln)) {
      float rq, rqround;
      if (baln.GetTag("rq", rq)) {
         if (filter(rq)) {
            // save bam alignment into output file
            bwriter.SaveAlignment(baln);
         }
      }
      else {
         throw runtime_error("no rq tag in sequence");;
      }
   }
   breader.Close();
   bwriter.Close();
   cerr << "Filtered result written to " << ofile << endl;
   cout << ofile << endl; // to be used by shell programming
}

void countAccuracy(const string &bamfile) {
   BamTools::BamReader breader;
   if (!breader.Open(bamfile)) {
      throw runtime_error("Failed to open bam file: " + bamfile);
   }
   BamTools::BamAlignment baln;
   map<float, int> rqcount;
   while (breader.GetNextAlignment(baln)) {
      float rq, rqround;
      if (baln.GetTag("rq", rq)) {
         rqround=round(rq*100)/100;
         //cout << rqround << endl;
         ++rqcount[rqround];
      }
   }
   breader.Close();
   int sum=0;
   for (auto it=rqcount.begin(); it != rqcount.end(); ++it) {
      cout << it->first << " => " << it->second << endl;
      sum += it->second;
   }
   cout << sum << " total sequences\n";
}
