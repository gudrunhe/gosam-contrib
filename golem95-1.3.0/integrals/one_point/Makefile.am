noinst_LTLIBRARIES=libgolem95_integrals_one_point.la

include ../../MakefileCommon.am
AM_FCFLAGS+=\
	-I$(builddir)/../../module \
	-I$(builddir)/../../kinematic \
	-I$(builddir)/../../numerical \
	-I$(builddir)/../../../avh_olo-2.2.1

libgolem95_integrals_one_point_la_SOURCES= generic_function_1p.f90
libgolem95_integrals_one_point_la_FCFLAGS= $(AM_FCFLAGS)

nodist_pkginclude_HEADERS= generic_function_1p.mod
CLEANFILES=*.mod

include Makefile.dep
