noinst_LTLIBRARIES=libgolem95_integrals_one_point.la

include ../../MakefileCommon.am
AM_FCFLAGS+=\
	-I$(builddir)/../../module \
	-I$(builddir)/../../kinematic \
	-I$(builddir)/../../numerical \
	-I$(builddir)/../../../avh_olo-3.6.1

libgolem95_integrals_one_point_la_SOURCES= generic_function_1p.f90

nodist_pkginclude_HEADERS= generic_function_1p.mod
CLEANFILES=*.mod

include Makefile.dep
