lib_LTLIBRARIES=libqcdloop.la

noinst_HEADERS=qlconstants.f qlonshellcutoff.f

libqcdloop_la_SOURCES= \
   qlinit.f qlI4.f qlI3.f qlI2.f qlI1.f qlI4array.f qlzero.f qlI4fin.f \
   qlI4DNS41.f qlI4sub0m.f qlI4sub1m.f qlI4sub2m.f qlI4sub2ma.f qlI4sub2mo.f \
   qlI4sub3m.f qlI3fin.f qlI3sub.f qlI2fin.f qlfndd.f qlLi2omprod.f \
   qlLi2omrat.f qllnomrat4.f qlLi2omx.f qlLi2omx2.f qltrisort.f qlsnglsort.f \
   qltri1.f qltri2.f qltri3.f qltri4.f qltri5.f qltri6.f qlbox1.f qlbox2.f \
   qlbox3.f qlbox4.f qlbox5.f qlbox6.f qlbox7.f qlbox8.f qlbox9.f qlbox10.f \
   qlbox11.f qlbox12.f qlbox13.f qlbox14.f qlbox15.f qlbox16.f qlcLi2omx2.f \
   qlcLi2omx3.f auxCD.f qlspencer.f qlratreal.f qlratgam.f ddilog.f qllnrat.f \
   qlfunctions.f qlYcalc.f qlkfn.f qlxpicheck.f

libqcdloop_la_LIBADD=-L$(builddir)/../ff-2.0 -lff
AM_FFLAGS=-I$(srcdir) -I$(srcdir)/../ff-2.0

include Makefile.dep
