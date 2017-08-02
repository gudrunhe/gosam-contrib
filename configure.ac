AC_INIT([GoSam Convenience Package],[2.0],[gosam@hepforge.org],[gosam-contrib],[http://projects.hepforge.org/golem/gosam-contrib/])
AC_PREREQ(2.65)

AC_CONFIG_HEADERS([config.h ninja-1.2.0/ninja/ninja_config.h])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([config.aux])
AC_CONFIG_SRCDIR([ninja-1.2.0/ninja_implem.cxx])
m4_include([m4/m4_ax_cxx_compile_stdcxx_11.m4])
m4_include([m4/m4_ax_cxx_tryandadd_flag.m4])
m4_include([m4/m4_ax_cxx_quadmath.m4])


dnl -----------------------------------------------
dnl Package name and version number (user defined)
dnl -----------------------------------------------
VERSION=2.0
PACKAGE=gosam-contrib


dnl -----------------------------------------------
dnl ---#[ Checks for programs:
dnl -----------------------------------------------
AM_INIT_AUTOMAKE([no-define no-dependencies subdir-objects ])

AC_PROG_SED

AC_LANG_PUSH([Fortran])

AC_PROG_CXX
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

AM_CONDITIONAL(COMPILE_SAMURAI, [test "x$enable_samurai" != xno])

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
dnl ---#[ Configure Ninja:
dnl ---------------------------------------------------------------------
NINJAVERSION=1.1.0

AC_DEFINE_DIR([DATADIR], [datadir], [Where data are placed to.])

AC_ARG_ENABLE([ninja],
   AS_HELP_STRING([--disable-ninja], [do not compile and install Ninja]))
AM_CONDITIONAL(COMPILE_NINJA, [test "x$enable_ninja" != xno])

AS_IF([test "x$enable_ninja" != xno],
		[AC_SUBST([conf_with_ninja], [""])],
		[AC_SUBST([conf_with_ninja], ["# "])])
AS_IF([test "x$enable_ninja" != xno],
		[AC_SUBST([conf_wout_ninja], ["# "])],
		[AC_SUBST([conf_wout_ninja], [""])])
AS_IF([test "x$enable_ninja" != xno],
		[AC_SUBST([case_with_ninja], ["    "])],
		[AC_SUBST([case_with_ninja], ["!AC!"])])
AS_IF([test "x$enable_ninja" != xno],
		[AC_SUBST([case_wout_ninja], ["!AC!"])],
		[AC_SUBST([case_wout_ninja], ["    "])])

# From now on it is C++ stuff
AC_LANG_PUSH([C++])

ACX_PTHREAD

# Check C++11 compiler support
AX_CXX_COMPILE_STDCXX_11
AX_CXX_COMPILE_STDCXX_11_INIT_LIST
AX_CXX_COMPILE_STDCXX_11_RANDOM


# Check support for optional compiler flags
AC_ARG_ENABLE([autoflags],
    AS_HELP_STRING([--disable-autoflags],
    [do not add automatic compiler flags]))
AS_IF([test "x$enable_autoflags" != xno],
            [AX_CXX_TRYANDADD_FLAG("-fcx-fortran-rules")
            AX_CXX_TRYANDADD_FLAG("-complex-limited-range")]
)

# Choose whether the interface to the OneLoop library should be
# compiled
#  AC_ARG_WITH([avholo],
#      AS_HELP_STRING([--with-avholo@<:@=FLAGS@:>@],
#      [include an interface with the OneLoop integral library,\
#      specifying the corresponding flags for dynamic linking \
#      (if the Fortran module 'avh_olo' is not in a standard path, \
#      you should add its directory to the FCINCLUDE variable)]),
#      [],
#      [with_avholo=no])
#  AS_IF([test "x$with_avholo" == xyes],
#         [AC_SUBST([AVHFLAGS], ["-lavh_olo"])
#          AC_DEFINE([NINJA_USE_ONELOOP], [1],
#          [Define to 1 if you have OneLoop.])
#          [avholo=true]
#          ],
#         [test "x$with_avholo" != xno],
#         [AC_SUBST([AVHFLAGS], ["$with_avholo"])
#          AC_DEFINE([NINJA_USE_ONELOOP], [1],
#          [Define to 1 if you have OneLoop.])
#          [avholo=true]
#          ],
#          [avholo=false])
#  AM_CONDITIONAL(AVHOLO, test "x$avholo" == xtrue)

AS_IF([test "x$enable_olo" != xno],
       [AC_SUBST([AVHFLAGS], ["-lavh_olo"])
        AC_DEFINE([NINJA_USE_ONELOOP], [1],
        [Define to 1 if you have OneLoop.])
        [avholo=true]
        ],
       [test "x$with_avholo" != xno],
       [AC_SUBST([AVHFLAGS], ["$with_avholo"])
        AC_DEFINE([NINJA_USE_ONELOOP], [1],
        [Define to 1 if you have OneLoop.])
        [avholo=true]
        ],
        [avholo=false])

AM_CONDITIONAL(AVHOLO, [test "x$enable_olo" != xno])



# Choose whether an internal integral cache for the OneLoop library
# should be used
AC_ARG_ENABLE([avholo_cache],
    AS_HELP_STRING([--disable-avholo_cache],
    [do not include a cache of Master Integrals for OneLoop]))
AS_IF([test "x$enable_avholo_cache" != xno],
      [AC_DEFINE([NINJA_USE_ONELOOP_WITH_CACHE], [1],
      [Define to 1 if a cache of Master Integrals is implemented for OneLoop.])
      ],
      [])

# Choose whether the interface to the LoopTools library should be
# compiled
AC_ARG_WITH([looptools],
   AS_HELP_STRING([--with-looptools@<:@=FLAGS@:>@],
                  [include an interface with the LoopTools library, \
                  specifying the corresponding flags for static linking. \
                  (if the header file cloopttols.h is not in a standard \
                  path, you should add its directory to the \
                  CPPFLAGS variable)]),
                  [],
                  [with_looptools=no])
AS_IF([test "x$with_looptools" == xyes],
       [AC_SUBST([LTFLAGS], ["-looptools"])
        AC_DEFINE([NINJA_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true]
        ],
       [test "x$with_looptools" != xno],
       [AC_SUBST([LTFLAGS], ["$with_looptools"])
        AC_DEFINE([NINJA_USE_LOOPTOOLS], [1], [Define if you have LoopTools.])
        [looptools=true]
        ],
        [looptools=false])
AM_CONDITIONAL(LOOPTOOLS, test "x$looptools" == xtrue)
AS_IF([test "x$looptools" == xtrue],
    [[case $host_os in
      darwin*)
         exninjastatic=true
    esac]])
AM_CONDITIONAL(EXNINJASTATIC, test "x$exninjastatic" == xtrue)

# Always enable massless specialization
AC_DEFINE([NINJA_MASSLESS], [1],
          [Define to 1 massless specialization is enabled.])

# Choose whether to include higher-rank support
AC_ARG_ENABLE([higher_rank],
    AS_HELP_STRING([--disable-higher_rank],
                   [disable support for higher-rank numerators]))
AS_IF([test "x$enable_higher_rank" != xno],
      [AC_DEFINE([NINJA_X1RANK], [1],
      [Define to 1 if higher-rank support is enabled.])])


# Test quadruple precision.  Intermediate precision is not compatible
# with Ninja.
AS_IF(
        [test "x$with_precision" == xquadruple],
        AC_DEFINE([NINJA_QUADRUPLE], [1],
    [Define to 1 if Ninja is compiled in quadruple precision.])
    AC_SUBST([QUADRUPLEFLAGS], ["-lquadmath"])
    [quadruple=true],
        [test "x$with_precision" == xquad],
        AC_DEFINE([NINJA_QUADRUPLE], [1],
    [Define to 1 if Ninja is compiled in quadruple precision.])
    AC_SUBST([QUADRUPLEFLAGS], ["-lquadmath"])
    [quadruple=true],
    [test "x$with_precision" == xintermediate],
    AS_IF([test "x$enable_ninja" != xno],
    [AC_MSG_FAILURE([intermediate precision for Ninja is not available])]),
        [quadruple=false])
AM_CONDITIONAL(QUADRUPLE, test "x$quadruple" == xtrue)

# Choose whether to compile in quadruple precision
# AC_ARG_WITH([quadruple],
#     AS_HELP_STRING([--with-quadruple@<:@=FLAGS@:>@],
#     [compile the library in quadruple precision (requires GCC libquadmath), \
#     specifying the corresponding flags for the linker]),
#     [],
#     [with_quadruple=no])

# Choose whether to compile quadninja
AC_ARG_ENABLE([quadninja],
    AS_HELP_STRING([--enable-quadninja@<:@=FLAGS@:>@],
    [build quadninja, a copy of ninja in quadruple precision \
    (requires GCC libquadmath) alongside ninja, \
    specifying the corresponding flags for the linker]),
    [],
    [enable_quadninja=no])

# Check support for quadruple precision if needed
# AS_IF([test "x$with_quadruple" != xno] || [test "x$enable_quadninja" != xno],
#     [AX_CXX_COMPILE_LIBQUADMATH]
# )
# AS_IF([test "x$with_quadruple" != xno] || [test "x$enable_quadninja" != xno],
#     [AX_CXX_COMPILE_STDCXX_11_FLOAT128_RANDOM]
# )
AS_IF([test "x$enable_quadninja" != xno],[AX_CXX_COMPILE_LIBQUADMATH])
AS_IF([test "x$enable_quadninja" != xno],[AX_CXX_COMPILE_STDCXX_11_FLOAT128_RANDOM])

## Quadruple precision
#AS_IF([test "x$with_quadruple" == xyes],
#       AC_SUBST([QUADRUPLEFLAGS], ["-lquadmath"])
#       AC_DEFINE([NINJA_QUADRUPLE], [1],
#       [Define to 1 if Ninja is compiled in quadruple precision.])
#       [quadruple=true],
#       [test "x$with_quadruple" != xno],
#       AC_SUBST([QUADRUPLEFLAGS], ["$with_quadruple"])
#       AC_DEFINE([NINJA_QUADRUPLE], [1],
#       [Define to 1 if Ninja is compiled in quadruple precision.])
#       [quadruple=true],
#       [quadruple=false])
#AM_CONDITIONAL(QUADRUPLE, test "x$quadruple" == xtrue)

# Quadninja
AS_IF([test "x$enable_quadninja" == xyes],
       AC_SUBST([QUADRUPLEFLAGS], ["-lquadmath"])
       AC_DEFINE([QUADNINJA], [1],
       [Define to 1 if QuadNinja is built.])
       [quadninja=true],
       [test "x$enable_quadninja" != xno],
       AC_SUBST([QUADRUPLEFLAGS], ["$enable_quadninja"])
       AC_DEFINE([QUADNINJA], [1],
       [Define to 1 if QuadNinja is built.])
       [quadninja=true],
       [quadninja=false])
AM_CONDITIONAL(QUADNINJA, test "x$quadninja" == xtrue)

# Set up macros for QuadNinja
AS_IF([test "x$enable_quadninja" != xno],
    AS_IF([test "x$with_precision" == xquadruple],
    AC_MSG_ERROR(Option --with-precision=quadruple is not compatible with option --enable-quadninja))
    AC_DEFINE([QUADNINJA_CONFIG_H_INTERNAL], [1],
              [Macro for QuadNinja])
    AS_IF([test "x$with_avholo" != xno],
    AC_DEFINE([QUADNINJA_USE_ONELOOP], [1],
              [Macro for QuadNinja]))
    AS_IF([test "x$enable_avholo_cache" != xno],
    AC_DEFINE([QUADNINJA_USE_ONELOOP_WITH_CACHE], [1],
              [Macro for QuadNinja]))
    #AS_IF([test "x$with_looptools" != xno],
    #AC_DEFINE([QUADNINJA_USE_LOOPTOOLS], [NINJA_USE_LOOPTOOLS],
    #          [Macro for QuadNinja]))
    AC_DEFINE([QUADNINJA_MASSLESS], [1],
              [Macro for QuadNinja])
    AS_IF([test "x$enable_higher_rank" != xno],
    AC_DEFINE([QUADNINJA_X1RANK], [1],
              [Macro for QuadNinja]))
)

# Choose whether compiling interfaces: f90 module
AC_ARG_ENABLE([f90module],
    AS_HELP_STRING([--disable-f90module], [do not include Fortran-90 interface of Ninja]))
AS_IF([test "x$enable_f90module" != xno], [f90module=true],
 [f90module=false]
  )
AM_CONDITIONAL(F90MODULE, test "x$f90module" == xtrue)

AC_ARG_VAR([FCINCLUDE],[additional Fortran include flags, e.g. \
                        -I<include dir> if you have fortran modules \
                        installed in a nonstandard directory <include dir>])


# These are always defined for the contrib
AC_DEFINE([NINJA_CONFIG_H_INTERNAL], [1],
          [Define to 1 if this an internal is a source file.])
# AC_DEFINE([NINJA_X1RANK], [1],
#           [Define to 1 if higher-rank support is enabled.])
AM_CONDITIONAL(GOSAM, true)
AM_CONDITIONAL(NINJAGOSAMCONTRIB, true)

AC_DEFINE([NINJA_GOSAM],[1],[Define to 1 if GoSam is used.])

AC_DEFINE_UNQUOTED([VERSION],["$NINJAVERSION"],[Version of Ninja])

AC_SUBST(NINJAVERSION)
dnl ---------------------------------------------------------------------
dnl ---#] Configure Ninja:
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
        qcdloop-1.95/Makefile \
        avh_olo-3.6.1/Makefile \
        avh_olo-3.6.1/avh_olo_xkind.f90 \
        golem95-1.3.3/Makefile \
        golem95-1.3.3/module/Makefile \
        golem95-1.3.3/module/precision_golem.f90 \
        golem95-1.3.3/integrals/Makefile \
        golem95-1.3.3/integrals/two_point/Makefile \
        golem95-1.3.3/integrals/one_point/Makefile \
        golem95-1.3.3/integrals/three_point/Makefile \
        golem95-1.3.3/integrals/four_point/Makefile \
        golem95-1.3.3/integrals/four_point/generic_function_4p.f90 \
        golem95-1.3.3/higherrank/Makefile \
        golem95-1.3.3/interface/Makefile \
        golem95-1.3.3/numerical/Makefile \
        golem95-1.3.3/kinematic/Makefile \
        golem95-1.3.3/form_factor/Makefile \
        samurai-2.9.1/Makefile \
        samurai-2.9.1/mmisavholo.f90 \
        samurai-2.9.1/mmisqcdloop.f90 \
        samurai-2.9.1/mmisgolem.f90 \
        samurai-2.9.1/mmislooptools.f90 \
        samurai-2.9.1/mmasters.f90 \
        samurai-2.9.1/precision.f90 \
        ninja-1.2.0/Makefile \
        ninja-1.2.0/ninja/Makefile \
	ninja-1.2.0/quadninja/Makefile
])
AC_OUTPUT
dnl ---------------------------------------------------------------------
dnl ---#] Generates Makefile's, configuration files and scripts:
dnl ---------------------------------------------------------------------
