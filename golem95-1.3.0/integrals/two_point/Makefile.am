noinst_LTLIBRARIES=libgolem95_integrals_two_point.la

include ../../MakefileCommon.am
AM_FCFLAGS+=\
	-I$(builddir)/../../module \
	-I$(builddir)/../../kinematic \
	-I$(builddir)/../../numerical \
	-I$(builddir)/../one_point

libgolem95_integrals_two_point_la_SOURCES= \
	function_2p0m_1mi.f90 function_2p_m1m2.f90 generic_function_2p.f90
libgolem95_integrals_two_point_la_FCFLAGS= $(AM_FCFLAGS)

nodist_pkginclude_HEADERS= function_2p_m1m2.mod function_2p0m_1mi.mod \
	generic_function_2p.mod
CLEANFILES=*.mod

include Makefile.dep
