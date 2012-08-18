lib_LTLIBRARIES=libsamurai.la
libsamurai_la_SOURCES=\
   constants.f90 kinematic.f90 ltest.f90 madds.f90 mcgs.f90 mfunctions.f90 \
   mgetbase.f90 mgetc1.f90 mgetc2.f90 mgetc3.f90 mgetc4.f90 mgetc5.f90 \
   mgetqs.f90 mglobal.f90 mrestore.f90 msamurai.f90 mtens.f90 mtests.f90 \
   ncuts.f90 notfirst.f90 options.f90 precision.f90 save.f90 maccu.f90
nodist_pkginclude_HEADERS=\
   constants.mod kinematic.mod ltest.mod madds.mod mcgs.mod mfunctions.mod \
   mgetbase.mod mgetc1.mod mgetc2.mod mgetc3.mod mgetc4.mod mgetc5.mod \
   mgetqs.mod mglobal.mod mrestore.mod msamurai.mod mtens.mod mtests.mod \
   ncuts.mod notfirst.mod options.mod precision.mod save.mod maccu.mod

AM_FCFLAGS= \
	-I. \
	-I$(top_builddir)/avh_olo-2.2.1 \
	-I$(top_builddir)/golem95c-1.2.1/module

libsamurai_la_LIBADD=\
   $(LIBLOOPTOOLS)

if COMPILE_QL
libsamurai_la_LIBADD+=-L$(top_builddir)/qcdloop-1.9 -lqcdloop \
		      -L$(top_builddir)/ff-2.0 -lff
else
# nop
endif

if COMPILE_GOLEM95C
libsamurai_la_LIBADD+=-L$(top_builddir)/golem95c-1.2.1 -lgolem
else
# nop
endif

libsamurai_la_LIBADD+=-L$(top_builddir)/avh_olo-2.2.1 -lavh_olo

CLEANFILES=\
   constants.mod kinematic.mod ltest.mod madds.mod mcgs.mod mfunctions.mod \
   mgetbase.mod mgetc1.mod mgetc2.mod mgetc3.mod mgetc4.mod mgetc5.mod \
   mgetqs.mod mglobal.mod mrestore.mod msamurai.mod mtens.mod mtests.mod \
   ncuts.mod notfirst.mod options.mod precision.mod save.mod maccu.mod

include Makefile.dep
