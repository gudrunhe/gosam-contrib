# vim: ts=3:sw=3
SUBDIRS=

if COMPILE_FF
SUBDIRS+= ff-2.0
else
# nop
endif

if COMPILE_QL
SUBDIRS+= qcdloop-1.95
else
# nop
endif

if COMPILE_OLO
SUBDIRS+= avh_olo-3.6.1
else
# nop
endif

if COMPILE_GOLEM95C
SUBDIRS+= golem95-1.3.3
else
# nop
endif

if COMPILE_SAMURAI
SUBDIRS+= samurai-2.9.1
else
# nop
endif

if COMPILE_NINJA
SUBDIRS+= ninja-1.2.0
else
# nop
endif

pkgconfigdir= $(libdir)/pkgconfig
pkgconfig_DATA=samurai.pc

dist_noinst_SCRIPTS= autogen.sh

dist_pkgdata_DATA=\
	$(top_srcdir)/ff-2.0/fferr.dat \
	$(top_srcdir)/ff-2.0/ffperm5.dat \
	$(top_srcdir)/ff-2.0/ffwarn.dat

pkgdata_DATA=gosam.conf
ACLOCAL_AMFLAGS= -I m4
EXTRA_DIST=m4
