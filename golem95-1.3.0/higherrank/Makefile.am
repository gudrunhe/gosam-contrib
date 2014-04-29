noinst_LTLIBRARIES=libgolem95_higherrank.la

include ../MakefileCommon.am
AM_FCFLAGS+=\
	-I$(builddir)/../module \
	-I$(builddir)/../kinematic \
	-I$(builddir)/../numerical \
	-I$(builddir)/../form_factor \
	-I$(builddir)/../integrals/one_point \
	-I$(builddir)/../integrals/two_point \
	-I$(builddir)/../integrals/three_point \
	-I$(builddir)/../integrals/four_point

#	-I$(builddir)/../../avh_olo-2.2.1

libgolem95_higherrank_la_SOURCES=generic_function_np.f90 form_factor_higher_ranks.f90

nodist_pkginclude_HEADERS= generic_function_np.mod form_factor_higher_ranks.mod

CLEANFILES=*.mod


%.mod: %.o %.f90
	@true


include Makefile.dep
