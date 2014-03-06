lib_LTLIBRARIES = libninja.la

libninja_la_SOURCES = basis.cc  cuts_utils.cc  ninja.cc  Xcuts.cc \
	cuts.cc Xninja.cc spinors.cc \
	coefficient_level_subs.cc	\
	integral_library.cc avholo_interface.cc \
	rambo.cc
if AVHOLO
libninja_la_SOURCES += avholo_wrapper.F90
if NINJAGOSAMCONTRIB
AM_FCFLAGS = -I$(top_builddir)/avh_olo-2.2.1
else
AM_FCFLAGS = $(FCINCLUDE)
endif
endif
if LOOPTOOLS
libninja_la_SOURCES +=  looptools_interface.cc
endif
if GOSAM
libninja_la_SOURCES +=  gosam_interface.cc  ninjago.F90
endif

noinst_HEADERS = ninja_implem.cxx  polys_implem.cxx basis.hh \
    cuts.hh  cuts_utils.hh  Xcuts.hh  cuts_utils-inl.hh \
	coefficient_level_subs.hh	\
	integermath.hh  tmp_utils.hh  ninja_scoped_array.hh  cuts_vector.hh \
	integral_library_wrapper.hh s_mat_wrapper.hh avholo_decls.hh \
	ninja_hash_table.hh integral_cache.hh \
	gosam_interface.hh

library_includedir=$(includedir)
nodistdir=$(includedir)

libninja_la_LIBADD=

if GOSAM
if NINJAGOSAMCONTRIB
nodist_pkginclude_HEADERS = ninjago_module.mod
libninja_la_LIBADD += -L$(top_builddir)/avh_olo-2.2.1 -lavh_olo
else
nodist_library_include_HEADERS = ninjago_module.mod
endif
endif

libninja_la_LIBADD += $(AVHFLAGS) $(QUADRUPLEFLAGS)

ACLOCAL_AMFLAGS = -I m4
SUBDIRS = ninja
