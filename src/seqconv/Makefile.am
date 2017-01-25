bin_PROGRAMS = bam2cleanfastq checkdupfasq  checkdupfasta \
    cpfasqpart cpfastapart fastareduceid fastq2fasta fastq2tab fastqpick \
	 fastqreduceid fastqstat fishing linearize revcomp rna2dna tandemcount trimfastq

AM_CXXFLAGS=$(orpara_CXXFLAGS)
AM_LDFLAGS=$(orpara_LDFLAGS)
LDADD=$(orpara_LIBS)

bam2cleanfastq_SOURCES = bam2cleanfastq.cpp
bam2cleanfastq_CXXFLAGS = $(AM_CXXFLAGS) $(bam_CXXFLAGS)
bam2cleanfastq_LDFLAGS = $(AM_LDFLAGS) $(bam_LDFLAGS)
bam2cleanfastq_LDADD = $(LDADD) $(bam_LIBS)

checkdupfasq_SOURCES = checkdupfasq.cpp
checkdupfasta_SOURCES = checkdupfasta.cpp
cpfasqpart_SOURCES = cpfasqpart.cpp
cpfastapart_SOURCES = cpfastapart.cpp
fastareduceid_SOURCES = fastareduceid.cpp
fastq2fasta_SOURCES = fastq2fasta.cpp
fastq2tab_SOURCES = fastq2tab.cpp
fastqpick_SOURCES = fastqpick.cpp
fastqreduceid_SOURCES = fastqreduceid.cpp
fastqstat_SOURCES = fastqstat.cpp

fishing_SOURCES = fishing.cpp 
fishing_CXXFLAGS = $(AM_CXXFLAGS) $(boost_CXXFLAGS)
fishing_LDFLAGS = $(AM_LDFLAGS) $(boost_LDFLAGS)
fishing_LDADD = $(LDADD) $(boost_LIBS)

linearize_SOURCES = linearize.cpp 
linearize_CXXFLAGS = $(AM_CXXFLAGS) $(boost_CXXFLAGS)
linearize_LDFLAGS = $(AM_LDFLAGS) $(boost_LDFLAGS)
linearize_LDADD = $(LDADD) $(boost_LIBS)

revcomp_SOURCES = revcomp.cpp
rna2dna_SOURCES = rna2dna.cpp
tandemcount_SOURCES = tandemcount.cpp
trimfastq_SOURCES = trimfastq.cpp 
