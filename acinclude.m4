
# MY_F77_LINE_LENGTH([LENGTH], [ACTION-IF-SUCCESS],
#		    [ACTION-IF-FAILURE = FAILURE])
# ------------------------------------------------
# Look for a compiler flag to make the Fortran (FC) compiler accept long lines
# in the current (free- or fixed-format) source code, and adds it to FFLAGS.
# The optional LENGTH may be 80, 132 (default), or `unlimited' for longer
# lines.  Note that line lengths above 254 columns are not portable, and some
# compilers (hello ifort) do not accept more than 132 columns at least for
# fixed format.  Call ACTION-IF-SUCCESS (defaults to nothing) if successful
# (i.e. can compile code using new extension) and ACTION-IF-FAILURE (defaults
# to failing with an error message) if not.  (Defined via DEFUN_ONCE to
# prevent flag from being added to FFLAGS multiple times.)
# You should call AC_FC_FREEFORM or AC_FC_FIXEDFORM to set the desired format
# prior to using this macro.
#
# The known flags are:
# -ffixed-line-length-N with N 72, 80, 132, or 0 or none for none.
# -ffixed-line-length-none: GNU gfortran
#       -qfixed=132 80 72: IBM compiler (xlf)
#                -Mextend: Cray
#            -132 -80 -72: Intel compiler (ifort)
#                          Needs to come before -extend_source because ifort
#                          accepts that as well with an optional parameter and
#                          doesn't fail but only warns about unknown arguments.
#          -extend_source: SGI compiler
#     -W NN (132, 80, 72): Absoft Fortran
#          +extend_source: HP Fortran (254 in either form, default is 72 fixed,
#			   132 free)
#                   -wide: Lahey/Fujitsu Fortran (255 cols in fixed form)
#                      -e: Sun Fortran compiler (132 characters)
AC_DEFUN_ONCE([MY_F77_LINE_LENGTH],
[AC_LANG_PUSH([Fortran 77])dnl
m4_case(m4_default([$1], [132]),
  [unlimited], [my_f77_line_len_string=unlimited
	               my_f77_line_len=0
                       my_f77_line_length_test='
      subroutine longer_than_132(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,'\
'arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,arg17,arg18,arg19)'],
  [132],            [my_f77_line_len=132
		       my_f77_line_length_test='
      subroutine longer_than_80(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,'\
'arg10)'],
  [80],             [my_f77_line_len=80
		       my_f77_line_length_test='
      subroutine longer_than_72(arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9)'],
  [m4_warning([Invalid length argument `$1'])])
: ${my_f77_line_len_string=$my_f77_line_len}
AC_CACHE_CHECK(
[for Fortran 77 flag needed to accept $my_f77_line_len_string column source lines],
	       [my_cv_f77_line_length],
[my_cv_f77_line_length=unknown
my_f77_line_length_FFLAGS_save=$FFLAGS
for ac_flag in none \
	        -ffixed-line-length-none \
	       -ffixed-line-length-$my_f77_line_len \
	       -qfixed=$my_f77_line_len -Mextend \
	       -$my_f77_line_len -extend_source \
	       "-W $my_f77_line_len" +extend_source -wide -e
do
  test "x$ac_flag" != xnone && FFLAGS="$my_f77_line_length_FFLAGS_save $ac_flag"
  AC_COMPILE_IFELSE([[$my_f77_line_length_test
      end subroutine]],
		    [my_cv_f77_line_length=$ac_flag; break])
done
rm -f conftest.err conftest.$ac_objext conftest.$ac_ext
FFLAGS=$my_f77_line_length_FFLAGS_save
])
if test "x$my_cv_f77_line_length" = xunknown; then
  m4_default([$3],
	     [AC_MSG_ERROR([Fortran does not accept long source lines], 77)])
else
  if test "x$my_cv_f77_line_length" != xnone; then
    FFLAGS="$FFLAGS $my_cv_f77_line_length"
  fi
  $2
fi
AC_LANG_POP([Fortran 77])dnl
])# MY_F77_LINE_LENGTH
