#AUTOMAKE_OPTIONS = subdir-objects

headers = stagergap.h nucaln.h MismatchMapping.h match.h stagergap.h
sources = nucaln.cpp MismatchMapping.cpp match.cpp stagergap.cpp

#pkginclude_HEADERS = $(headers)
noinst_LTLIBRARIES = liborparappcommon.la
liborparappcommon_la_SOURCES = $(headers) $(sources)

SUBDIRS = . consensus seqana seqconv test
