AM_CXXFLAGS = -I.. $(orpara_CXXFLAGS) 
AM_LDFLAGS = -L.. $(orpara_LDFLAGS)
LDADD = ../liborparappcommon.la $(orpara_LIBS)
bin_PROGRAMS = selortho GCcontent stemloop comparetoref picksubseq bam2tab tabulatepair validateprtsite polymerstat comparefastq bamfilter

selortho_SOURCES = selortho.cpp 

GCcontent_SOURCES =  GCcontent.cpp

stemloop_SOURCES = stemloop.cpp 
stemloop_CXXFLAGS = $(AM_CXXFLAGS) $(bam_CXXFLAGS) $(boost_CXXFLAGS)
stemloop_LDFLAGS = $(AM_LDFLAGS) $(bam_LDFLAGS) $(boost_LDFLAGS) 
stemloop_LDADD = $(LDADD) $(bam_LIBS) $(boost_LIBS)

comparetoref_SOURCES = comparetoref.cpp genoanalyzer.h genoanalyzer.cpp

picksubseq_SOURCES = picksubseq.cpp

bam2tab_SOURCES = bam2tab.cpp
bam2tab_CXXFLAGS = $(bam_CXXFLAGS) $(AM_CXXFLAGS)
bam2tab_LDFLAGS = $(bam_LDFLAGS) $(AM_LDFLAGS)
bam2tab_LDADD = $(LDADD) $(bam_LIBS)

tabulatepair_SOURCES = tabulatepair.cpp

validateprtsite_SOURCES = validateprtsite.cpp

polymerstat_SOURCES = polymerstat.cpp

comparefastq_SOURCES = comparefastq.cpp

bamfilter_SOURCES = bamfilter.cpp
bamfilter_CXXFLAGS = $(bam_CXXFLAGS) $(AM_CXXFLAGS)
bamfilter_LDFLAGS = $(bam_LDFLAGS) $(AM_LDFLAGS)
bamfilter_LDADD = $(LDADD) $(bam_LIBS)

SUBDIRS = . cdna seqsep fasta
