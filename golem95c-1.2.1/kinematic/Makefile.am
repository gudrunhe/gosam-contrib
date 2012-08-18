noinst_LTLIBRARIES=libgolem95_kinematics.la


libgolem95_kinematics_la_SOURCES= inverse_matrice.f90 matrice_s.f90
libgolem95_kinematics_la_FCFLAGS=\
	-I$(builddir)/../module \
	-I$(builddir)/../../avh_olo-2.2.1

nodist_pkginclude_HEADERS= matrice_s.mod inverse_matrice.mod
CLEANFILES=*.mod

include Makefile.dep
