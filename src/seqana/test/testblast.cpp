#include <iostream>

int main(int argc, char* argv[]) {
   string blastCommand="/remote/RSU/sw-cache/metag/bin/blastn -db silva123.fas -num_threads 8 -out C_ccssidconsensusswarmbln.tab -evalue 0.001 -word_size 14 -perc_identity 91 -query C_ccssidconsensus.fas -max_target_seqs 1 -outfmt '6 qseqid sseqid pident qlen slen length qstart qend sstart send  bitscore'";

   return 0;
}
