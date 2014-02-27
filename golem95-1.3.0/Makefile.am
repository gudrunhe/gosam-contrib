SUBDIRS = module kinematic numerical integrals form_factor higherrank interface

lib_LTLIBRARIES= libgolem.la
libgolem_la_SOURCES=
libgolem_la_LIBADD=\
   $(builddir)/interface/libgolem95_interface.la \
   $(builddir)/form_factor/libgolem95_formfactor.la \
   $(builddir)/integrals/four_point/libgolem95_integrals_four_point.la \
   $(builddir)/integrals/three_point/libgolem95_integrals_three_point.la \
   $(builddir)/integrals/two_point/libgolem95_integrals_two_point.la \
   $(builddir)/integrals/one_point/libgolem95_integrals_one_point.la \
   $(builddir)/higherrank/libgolem95_higherrank.la \
   $(builddir)/numerical/libgolem95_numerical.la \
   $(builddir)/kinematic/libgolem95_kinematics.la \
   $(builddir)/module/libgolem95_module.la \
   -L$(builddir)/../avh_olo-2.2.1 -lavh_olo \
   $(LIBLOOPTOOLS)
