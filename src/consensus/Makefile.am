#AUTOMAKE_OPTIONS = subdir-objects

headers = dnaqualstore.h consensus.h boxpile.h piechart.h textpile.h stalactitegraph.h constaxmap.h blastnrow.h 
sources = dnaqualstore.cpp constaxmap.cpp blastnrow.cpp consensus.cpp boxpile.cpp stalactitegraph.cpp piechart.cpp textpile.cpp 

pkginclude_HEADERS = $(headers)
lib_LTLIBRARIES = liborparaconsensus.la
liborparaconsensus = $(headers) $(sources)
liborparaconsensus = -version-info 0:0:0
