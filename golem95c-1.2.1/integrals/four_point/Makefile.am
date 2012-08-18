noinst_LTLIBRARIES=libgolem95_integrals_four_point.la

libgolem95_integrals_four_point_la_SOURCES= \
	function_4p1m.f90 function_4p2m_3mi_onshell.f90 function_4p2m_adj.f90 \
	function_4p2m_opp.f90 function_4p3m.f90 function_4p4m.f90 \
	function_4p_ql10.f90  function_4p_ql11.f90 function_4p_ql12.f90 \
	function_4p_ql13.f90 function_4p_ql14.f90 function_4p_ql15.f90 \
	function_4p_ql16.f90 function_4p_ql6.f90 function_4p_ql7.f90 \
	function_4p_ql8.f90 function_4p_ql9.f90	generic_function_4p.f90
libgolem95_integrals_four_point_la_FCFLAGS=\
	-I$(builddir)/../../module \
	-I$(builddir)/../../kinematic \
	-I$(builddir)/../../numerical \
	-I$(builddir)/../one_point \
	-I$(builddir)/../two_point \
	-I$(builddir)/../three_point \
	-I$(builddir)/../../../avh_olo-2.2.1

nodist_pkginclude_HEADERS= function_4p4m.mod function_4p_ql7.mod \
	generic_function_4p.mod function_4p2m_adj.mod function_4p_ql14.mod \
	function_4p_ql13.mod function_4p_ql15.mod function_4p_ql10.mod \
	function_4p_ql9.mod function_4p_ql8.mod \
	function_4p2m_3mi_onshell.mod function_4p_ql11.mod \
	function_4p_ql12.mod function_4p2m_opp.mod function_4p3m.mod \
	function_4p1m.mod function_4p_ql16.mod function_4p_ql6.mod
CLEANFILES=*.mod
