noinst_LTLIBRARIES=libgolem95_formfactor.la

libgolem95_formfactor_la_SOURCES= form_factor_1p.f90 form_factor_2p.f90 \
	form_factor_3p.f90 form_factor_4p.f90 form_factor_5p.f90 \
	form_factor_6p.f90
libgolem95_formfactor_la_FCFLAGS=\
	-I$(builddir)/../module \
	-I$(builddir)/../kinematic \
	-I$(builddir)/../numerical \
	-I$(builddir)/../interface \
	-I$(builddir)/../integrals/one_point \
	-I$(builddir)/../integrals/two_point \
	-I$(builddir)/../integrals/three_point \
	-I$(builddir)/../integrals/four_point \
	-I$(builddir)/../../avh_olo-2.2.1

nodist_pkginclude_HEADERS= form_factor_5p.mod form_factor_6p.mod form_factor_4p.mod \
	form_factor_3p.mod form_factor_2p.mod form_factor_1p.mod
CLEANFILES=*.mod

include Makefile.dep
