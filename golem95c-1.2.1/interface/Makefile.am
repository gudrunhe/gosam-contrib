noinst_LTLIBRARIES=libgolem95_interface.la

AM_FCFLAGS= \
	-I$(builddir)/../module \
	-I$(builddir)/../kinematic \
	-I$(builddir)/../form_factor \
	-I$(builddir)/../integrals/one_point \
	-I$(builddir)/../integrals/two_point \
	-I$(builddir)/../integrals/three_point \
	-I$(builddir)/../integrals/four_point \
	-I$(builddir)/../numerical

libgolem95_interface_la_SOURCES= \
	tool_lt_to_golem.f90 tensor_integrals.f90 \
	gb0.f90 gc0.f90 gd0.f90 ge0.f90 gf0.f90

libgolem95_interface_la_FCFLAGS=$(AM_FCFLAGS)
nodist_pkginclude_HEADERS= tool_lt_to_golem.mod tensor_integrals.mod
if COMPILE_TENSREC
libgolem95_interface_la_SOURCES+=tens_rec.f90 tens_comb.f90
nodist_pkginclude_HEADERS+=tens_rec.mod tens_comb.mod
endif

CLEANFILES=*.mod

include Makefile.dep
