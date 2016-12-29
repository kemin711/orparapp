bin_PROGRAMS = uniqueseq fastabund fastaleninfo combineseq

uniqueseq_SOURCES = uniqueseq.cpp

fastabund_SOURCES = fastabund.cpp 
fastabund_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/consensus $(plotter_CXXFLAGS)
fastabund_LDFLAGS = -L$(top_srcdir)/src $(plotter_LDFLALGS)
fastabund_LDADD = $(top_srcdir)/src/consensus/liborparaconsensus.la $(plotter_LIBS) -lpthread

fastaleninfo_SOURCES = fastaleninfo.cpp

combineseq_SOURCES = combineseq.cpp

