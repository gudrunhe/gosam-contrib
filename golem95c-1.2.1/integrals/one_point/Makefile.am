noinst_LTLIBRARIES=libgolem95_integrals_one_point.la

libgolem95_integrals_one_point_la_SOURCES= generic_function_1p.f90
libgolem95_integrals_one_point_la_FCFLAGS=\
	-I$(builddir)/../../module \
	-I$(builddir)/../../kinematic \
	-I$(builddir)/../../numerical \
	-I$(builddir)/../../../avh_olo-2.2.1


nodist_pkginclude_HEADERS= generic_function_1p.mod
CLEANFILES=*.mod

include Makefile.dep
