dnl Origin: https://stackoverflow.com/questions/1354996/need-an-autoconf-macro-that-detects-if-m64-is-a-valid-compiler-option
dnl Modified by Johann Felix v. Soden-Fraunhofen

dnl @synopsis CXX_FLAGS_CHECK [compiler flags] [result_variable]
dnl @summary check whether compiler supports given C++ flags or not
AC_DEFUN([CXX_FLAG_CHECK],
[dnl
   AC_MSG_CHECKING([if $CXX supports/ignores $1])
   AC_LANG_PUSH([C++])
  ac_saved_cxxflags="$CXXFLAGS"
  CXXFLAGS="-Werror $1"
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([])],
    [AC_MSG_RESULT([yes])
     $2=yes],
    [AC_MSG_RESULT([no])
     $2=no]
  )
  CXXFLAGS="$ac_saved_cxxflags"
  AC_LANG_POP([C++])
])
