noinst_LTLIBRARIES=libgolem95_kinematics.la


libgolem95_kinematics_la_SOURCES= inverse_matrice.f90 matrice_s.f90

#libgolem95_kinematics_la_FCFLAGS=

include ../MakefileCommon.am

AM_FCFLAGS+= -I$(builddir)/../module \
	-I$(builddir)/../../avh_olo-2.2.1


nodist_pkginclude_HEADERS= inverse_matrice.mod matrice_s.mod
CLEANFILES=*.mod

include Makefile.dep

%.mod: %.o %.f90
	@true
