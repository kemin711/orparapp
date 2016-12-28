bin_PROGRAMS = uniqueseq fastabund fastaleninfo combineseq

uniqueseq_SOURCES = uniqueseq.cpp
#uniqueseq_LDADD = -lkemin

fastabund_SOURCES = fastabund.cpp 
fastabund_CXXFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/consensus $(plotter_CXXFLAGS)
fastabund_LDFLAGS = -L$(top_srcdir)/src $(plotter_LDFLALGS)
fastabund_LDADD = $(top_srcdir)/src/consensus/liborparapp.la $(plotter_LIBS) -lpthread

fastaleninfo_SOURCES = fastaleninfo.cpp
#fastaleninfo_LDADD = -lkemin 

combineseq_SOURCES = combineseq.cpp
#combineseq_LDADD = -lkemin 

