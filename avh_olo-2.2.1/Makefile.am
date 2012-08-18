lib_LTLIBRARIES=libavh_olo.la
libavh_olo_la_SOURCES= avh_olo_boxc.f90 avh_olo_box.f90 avh_olo_bub.f90 \
   avh_olo_func.f90 avh_olo_kinds.f90 avh_olo_main.f90 avh_olo_print.f90 \
   avh_olo_tri.f90 avh_olo_units.f90 avh_olo_wrp01.f90 avh_olo_xkind.f90

noinst_HEADERS=\
   avh_olo_a0_a.h90 avh_olo_a0_b.h90 avh_olo_b0_a.h90 avh_olo_b0_b.h90 \
   avh_olo_b11_a.h90 avh_olo_b11_b.h90 avh_olo_c0_a.h90 avh_olo_c0_b.h90 \
   avh_olo_d0_a.h90 avh_olo_d0_b.h90
nodist_pkginclude_HEADERS=\
   avh_olo_bern.mod avh_olo_boxc.mod avh_olo_box.mod avh_olo_bub.mod \
   avh_olo_func.mod avh_olo_kinds.mod avh_olo_li2a.mod avh_olo_li2c2.mod \
   avh_olo_li2c.mod avh_olo_loga2.mod avh_olo_loga.mod avh_olo_logc2.mod \
   avh_olo_logc.mod avh_olo.mod avh_olo_print.mod avh_olo_tri.mod \
   avh_olo_units.mod avh_olo_xkind.mod

CLEANFILES=$(nodist_pkginclude_HEADERS)

include Makefile.dep

avh_olo_func.mod avh_olo_loga.mod avh_olo_bern.mod avh_olo_li2a.mod \
avh_olo_loga2.mod avh_olo_logc.mod avh_olo_li2c.mod avh_olo_logc2.mod \
avh_olo_li2c2.mod: avh_olo_func.o

avh_olo.mod: $(nodist_pkginclude_HEADERS:avh_olo.mod=)

%.mod: %.o %.f90
	@true
