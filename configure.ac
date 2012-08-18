AC_INIT([GoSam Convenience Package],[1.0],[reiterth@mpp.mpg.de],[gosam-contrib],[http://projects.hepforge.org/golem/gosam-contrib/])
AC_PREREQ(2.65)

AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config.aux])

dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------
VERSION=1.0
PACKAGE=gosam-contrib


dnl -----------------------------------------------
dnl ---#[ Checks for programs:
dnl -----------------------------------------------
AM_INIT_AUTOMAKE($PACKAGE, $VERSION, no-define, no-dependencies)

AC_LANG_PUSH([Fortran])

AC_PROG_F77([gfortran ifort g77 f77 xlf frt pgf77 fort77 fl32 af77])
AC_PROG_FC([gfortran ifort g77 f77 xlf frt pgf77 fort77 fl32 af77])

dnl Check the flags needed to link fc programs with ld (i.e. cc)
AC_FC_LIBRARY_LDFLAGS
dnl Check for underscoring of external names
AC_FC_WRAPPERS
dnl We need to use .f90 and not .f to enable Automake FC support
dnl Some Intel ifc/ifort do not understand .f95.  :-/
AC_FC_SRCEXT([f90])
dnl AC_FC_FIXEDFORM
dnl AC_FC_LINE_LENGTH([unlimited],[],[])
AC_FC_FREEFORM
dnl Enable long line support if available
AC_FC_LINE_LENGTH([unlimited],[],[AC_MSG_WARN([Fortran compiler does not accept long source lines], 77)])
AC_LANG_POP

AC_LANG_PUSH([Fortran 77])
MY_F77_LINE_LENGTH([132],[],[AC_MSG_WARN([Fortran 77 compiler does not accept long source lines], 77)])
AC_LANG_POP

LT_INIT
LT_INIT

AM_SANITY_CHECK
dnl -----------------------------------------------
dnl ---#] Checks for programs:
dnl -----------------------------------------------

dnl -----------------------------------------------
dnl ---#[ Setup datadir variable:
dnl -----------------------------------------------
dnl @synopsis AC_DEFINE_DIR(VARNAME, DIR [, DESCRIPTION])
dnl
dnl This macro sets VARNAME to the expansion of the DIR variable,
dnl taking care of fixing up ${prefix} and such.
dnl
dnl VARNAME is then offered as both an output variable and a C
dnl preprocessor symbol.
dnl
dnl Example:
dnl
dnl    AC_DEFINE_DIR([DATADIR], [datadir], [Where data are placed to.])
dnl
dnl @category Misc
dnl @author Stepan Kasal <kasal@ucw.cz>
dnl @author Andreas Schwab <schwab@suse.de>
dnl @author Guido U. Draheim <guidod@gmx.de>
dnl @author Alexandre Oliva
dnl @version 2006-10-13
dnl @license AllPermissive

AC_DEFUN([AC_DEFINE_DIR], [
   prefix_NONE=
   exec_prefix_NONE=
   test "x$prefix" = xNONE && prefix_NONE=yes && prefix=$ac_default_prefix
   test "x$exec_prefix" = xNONE && exec_prefix_NONE=yes && exec_prefix=$prefix
dnl In Autoconf 2.60, ${datadir} refers to ${datarootdir}, which in turn
dnl refers to ${prefix}.  Thus we have to use `eval' twice.
   eval ac_define_dir="\"[$]$2\""
   eval ac_define_dir="\"$ac_define_dir\""
   AC_SUBST($1, "$ac_define_dir")
   AC_DEFINE_UNQUOTED($1, "$ac_define_dir", [$3])
   test "$prefix_NONE" && prefix=NONE
   test "$exec_prefix_NONE" && exec_prefix=NONE
])

AC_DEFINE_DIR([DATADIR], [datadir])
dnl -----------------------------------------------
dnl ---#] Setup datadir variable:
dnl -----------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Set Precision:
dnl ---------------------------------------------------------------------
AC_ARG_WITH([precision], [AS_HELP_STRING([--with-precision],
	    [set the precision of the library to either
	     'double' or 'quadruple'. @<:@default=double@:>@])],
	     [], [with_precision=double])
AS_IF(
	[test "x$with_precision" == xquadruple],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_precision" == xquad],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_precision" == xdouble],
		[AC_SUBST([fortran_real_kind], ["kind(1.0d0)"])],
	[test "x$with_precision" == xintermediate],
		[AC_SUBST([fortran_real_kind], ["selected_real_kind(18,4931)"])],
	[AC_MSG_FAILURE([--with-precision was given with an unrecognized
	 parameter])])
dnl ---------------------------------------------------------------------
dnl ---#] Set Precision:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure FF:
dnl ---------------------------------------------------------------------
AC_ARG_ENABLE([ff],
   AS_HELP_STRING([--disable-ff], [do not compile and install FF]))

AM_CONDITIONAL(COMPILE_FF, [test "x$enable_ff" != xno])

AS_IF([test "x$enable_ff" != xno],
		[AC_SUBST([conf_with_ff], [""])],
		[AC_SUBST([conf_with_ff], ["# "])])
AS_IF([test "x$enable_ff" != xno],
		[AC_SUBST([conf_wout_ff], ["# "])],
		[AC_SUBST([conf_wout_ff], [""])])
AS_IF([test "x$enable_ff" != xno],
		[AC_SUBST([case_with_ff], ["    "])],
		[AC_SUBST([case_with_ff], ["!AC!"])])
AS_IF([test "x$enable_ff" != xno],
		[AC_SUBST([case_wout_ff], ["!AC!"])],
		[AC_SUBST([case_wout_ff], ["    "])])
dnl ---------------------------------------------------------------------
dnl ---#] Configure FF:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure QCDLoop:
dnl ---------------------------------------------------------------------
AC_ARG_ENABLE([ql],
   AS_HELP_STRING([--disable-ql], [do not compile and install QCDLoop]))

AM_CONDITIONAL(COMPILE_QL, [test "x$enable_ql" != xno])

AS_IF([test "x$enable_ql" != xno],
		[AC_SUBST([conf_with_ql], [""])],
		[AC_SUBST([conf_with_ql], ["# "])])
AS_IF([test "x$enable_ql" != xno],
		[AC_SUBST([conf_wout_ql], ["# "])],
		[AC_SUBST([conf_wout_ql], [""])])
AS_IF([test "x$enable_ql" != xno],
		[AC_SUBST([case_with_ql], ["    "])],
		[AC_SUBST([case_with_ql], ["!AC!"])])
AS_IF([test "x$enable_ql" != xno],
		[AC_SUBST([case_wout_ql], ["!AC!"])],
		[AC_SUBST([case_wout_ql], ["    "])])
dnl ---------------------------------------------------------------------
dnl ---#] Configure QCDLoop:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure AVH OneLOop:
dnl ---------------------------------------------------------------------
AC_ARG_ENABLE([olo],
   AS_HELP_STRING([--disable-olo], [do not compile and install AVH OneLOop]))

AM_CONDITIONAL(COMPILE_OLO, [test "x$enable_olo" != xno])
AC_SUBST([avh_olo_real_kind],["$fortran_real_kind"])

AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([conf_with_olo], [""])],
		[AC_SUBST([conf_with_olo], ["# "])])
AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([conf_wout_olo], ["# "])],
		[AC_SUBST([conf_wout_olo], [""])])
AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([case_with_olo], ["    "])],
		[AC_SUBST([case_with_olo], ["!AC!"])])
AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([case_wout_olo], ["!AC!"])],
		[AC_SUBST([case_wout_olo], ["    "])])
AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([case_with_avh], ["    "])],
		[AC_SUBST([case_with_avh], ["!AC!"])])
AS_IF([test "x$enable_olo" != xno],
		[AC_SUBST([case_wout_avh], ["!AC!"])],
		[AC_SUBST([case_wout_avh], ["    "])])
dnl ---------------------------------------------------------------------
dnl ---#] Configure AVH OneLOop:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure Golem95C:
dnl ---------------------------------------------------------------------
AC_ARG_ENABLE([golem95],
   AS_HELP_STRING([--disable-golem95], [do not compile and install Golem95C]))

AM_CONDITIONAL(COMPILE_GOLEM95C, [test "x$enable_golem95" != xno])

AS_IF([test "x$enable_golem95" != xno],
		[AC_SUBST([conf_with_golem95], [""])],
		[AC_SUBST([conf_with_golem95], ["# "])])
AS_IF([test "x$enable_golem95" != xno],
		[AC_SUBST([conf_wout_golem95], ["# "])],
		[AC_SUBST([conf_wout_golem95], [""])])
AS_IF([test "x$enable_golem95" != xno],
		[AC_SUBST([case_with_golem], ["    "])],
		[AC_SUBST([case_with_golem], ["!AC!"])])
AS_IF([test "x$enable_golem95" != xno],
		[AC_SUBST([case_wout_golem], ["!AC!"])],
		[AC_SUBST([case_wout_golem], ["    "])])

AM_CONDITIONAL(COMPILE_TENSREC, [test "x" == "x"])
dnl ---------------------------------------------------------------------
dnl ---#] Configure Golem95C:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure Samurai:
dnl ---------------------------------------------------------------------
SAMURAIVERSION=2.1.1

AC_ARG_ENABLE([samurai],
   AS_HELP_STRING([--disable-samurai], [do not compile and install Samurai]))

AM_CONDITIONAL(COMPILE_SAMURAI, [test "x$enable_golem95" != xno])

AS_IF([test "x$enable_samurai" != xno],
		[AC_SUBST([conf_with_samurai], [""])],
		[AC_SUBST([conf_with_samurai], ["# "])])
AS_IF([test "x$enable_samurai" != xno],
		[AC_SUBST([conf_wout_samurai], ["# "])],
		[AC_SUBST([conf_wout_samurai], [""])])
AS_IF([test "x$enable_samurai" != xno],
		[AC_SUBST([case_with_samurai], ["    "])],
		[AC_SUBST([case_with_samurai], ["!AC!"])])
AS_IF([test "x$enable_samurai" != xno],
		[AC_SUBST([case_wout_samurai], ["!AC!"])],
		[AC_SUBST([case_wout_samurai], ["    "])])
AC_SUBST(SAMURAIVERSION)
dnl ---------------------------------------------------------------------
dnl ---#] Configure Samurai:
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Configure LoopTools (external):
dnl ---------------------------------------------------------------------
AC_ARG_WITH([looptools],
  [AS_HELP_STRING([--with-looptools],
  [enable linking to LoopTools])],
  [],
  [with_looptools=no])

LIBLOOPTOOLS=
AS_IF(
  [test "x$with_looptools" == xyes],
  [AC_CHECK_LIB([ooptools -lgfortran], [ltexi_],
     [AC_SUBST([LIBLOOPTOOLS], ["-looptools"])
      AC_DEFINE([HAVE_LT], [1], [Define if you have Looptools])
     ],
     [AC_MSG_FAILURE(
        [--with-looptools was given, but test for -looptools failed. \
 Consider using --with-looptools=path/libooptools.a .])],
     [-looptools])],
   [test "x$with_looptools" != xno],
   [AC_CHECK_FILE([$with_looptools], 
     [AC_SUBST([LIBLOOPTOOLS], ["$with_looptools"])
      AC_DEFINE([HAVE_LT], [1], [Define if you have Looptools])
     ],
     [AC_MSG_FAILURE(
        [--with-looptools was given, but location '$with_looptools' \
 is wrong.])])]
  )

AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([conf_with_lt], [""])],
		[AC_SUBST([conf_with_lt], ["#"])])
AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([conf_wout_lt], ["#"])],
		[AC_SUBST([conf_wout_lt], [""])])

AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([case_with_lt], ["    "])],
		[AC_SUBST([case_with_lt], ["!AC!"])])
AS_IF([test "x$with_looptools" != xno],
		[AC_SUBST([case_wout_lt], ["!AC!"])],
		[AC_SUBST([case_wout_lt], ["    "])])

dnl -----------------------------------------------
dnl Set the precision used by LoopTools
dnl -----------------------------------------------
AC_ARG_WITH([lt-precision], [AS_HELP_STRING([--with-lt-precision],
	    [set the precision used by LoopTools to either
	     'double' or 'quadruple'. @<:@default=double@:>@])],
	     [], [with_lt_precision=double])

AS_IF(
	[test "x$with_lt_precision" == xquadruple],
		[AC_SUBST([lt_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_lt_precision" == xquad],
		[AC_SUBST([lt_real_kind], ["selected_real_kind(32,50)"])],
	[test "x$with_lt_precision" == xdouble],
		[AC_SUBST([lt_real_kind], ["kind(1.0d0)"])],
	[AC_MSG_FAILURE([--with-lt-precision was given with an
	 unrecognized parameter])])
dnl ---------------------------------------------------------------------
dnl ---#] Configure LoopTools (external):
dnl ---------------------------------------------------------------------

dnl ---------------------------------------------------------------------
dnl ---#[ Generates Makefile's, configuration files and scripts:
dnl ---------------------------------------------------------------------
AC_CONFIG_FILES([\
        gosam.conf \
	samurai.pc \
        Makefile \
        ff-2.0/Makefile \
        ff-2.0/ffinit.f \
        qcdloop-1.9/Makefile \
        avh_olo-2.2.1/Makefile \
        golem95c-1.2.1/Makefile \
        golem95c-1.2.1/module/Makefile \
        golem95c-1.2.1/module/precision_golem.f90 \
        golem95c-1.2.1/integrals/Makefile \
        golem95c-1.2.1/integrals/two_point/Makefile \
        golem95c-1.2.1/integrals/one_point/Makefile \
        golem95c-1.2.1/integrals/three_point/Makefile \
        golem95c-1.2.1/integrals/four_point/Makefile \
        golem95c-1.2.1/integrals/four_point/generic_function_4p.f90 \
        golem95c-1.2.1/interface/Makefile \
        golem95c-1.2.1/numerical/Makefile \
        golem95c-1.2.1/kinematic/Makefile \
        golem95c-1.2.1/form_factor/Makefile \
        samurai-2.1.1/Makefile \
        samurai-2.1.1/madds.f90 \
        samurai-2.1.1/msamurai.f90 \
        samurai-2.1.1/precision.f90 \
])
AC_OUTPUT
dnl ---------------------------------------------------------------------
dnl ---#] Generates Makefile's, configuration files and scripts:
dnl ---------------------------------------------------------------------
