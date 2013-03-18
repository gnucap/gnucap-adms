#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([gnucap-adms], [0.0.4], [felix@salfelder.org])
AC_CONFIG_SRCDIR([README])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wno-portability])
LT_INIT()

dnl ------------------- debugging -------------------
AC_MSG_CHECKING([if debug code should be compiled in])
AC_ARG_ENABLE([debug],
[  --enable-debug            Enable building of debug code. [[default: disabled]]],
[
if test "X$enable_debug" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_debug=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_debug=no
])

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_LANG(C++)

# Checks for libraries.

# Checks for header files.
GC_CPPFLAGS+=$(gnucap-conf --cppflags)

# Checks for typedefs, structures, and compiler characteristics.

dnl ----------------- admsXml ------------------
AC_CHECK_PROG(HAVE_GCCONF, gnucap-conf, true, false)

if test "x$HAVE_GCCONF" = "xfalse"; then
	AC_MSG_ERROR([cannot find gnucap-conf, not implemented])
fi

dnl ----------------- admsXml ------------------
AC_CHECK_PROG(HAVE_ADMSXML, admsXml, true, false)
if test "x$HAVE_GCCONF" = "xfalse"; then
	AC_MSG_ERROR([cannot find admsXml])
fi

AC_MSG_CHECKING([extern support])
ADMSXML_FLAGS=
ADMSXML_NOEXTERN=
# BUG in admsXml
# if admsXml -DUSE_EXTERN src/gnucap-adms.h 2>/dev/null >/dev/null; then
# doesnt work. using ugly check
if echo "extern module foo();endmodule" | admsXml /dev/stdin 2>/dev/null >/dev/null; then
	AC_MSG_RESULT([yes])
	AC_DEFINE([VERILOG_EXTERN], [1], [admsXml understands 'extern'])
	ADMSXML_FLAGS=-DUSE_EXTERN
	ADMSXML_NOEXTERN=[\[extern=\'no\'\]]
else
	AC_MSG_RESULT([no])
fi
AC_SUBST([ADMSXML_FLAGS])
AC_SUBST([ADMSXML_NOEXTERN])
AM_CONDITIONAL([VERILOG_EXTERN], [test x$ADMSXML_NOEXTERN != x] )

dnl ----------------- plugindir ------------------
AC_MSG_CHECKING([local])
AC_ARG_ENABLE([local],
[  --enable-local            install plugins locally. [[default: disabled]]],
[
if test "X$enable_local" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_local=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_local=no
])

if test "X$enable_local" = "Xno"; then
	gnucap_libdir=$(gnucap-conf --libdir)
else
	gnucap_libdir=$HOME/.gnucap/lib
fi

test "x$prefix" = xNONE &&  prefix=$ac_default_prefix
libdir=${prefix}/lib/gnucap-adms

AC_SUBST([gnucap_libdir])
AC_SUBST([GC_CPPFLAGS])
CXXFLAGS+=" -Wall -pedantic"

# Checks for library functions.

if test "$enable_debug" = "yes" ; then
	CPPFLAGS="$CPPFLAGS -DTRACE_UNTESTED"
else
	CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi

AC_CONFIG_FILES([Makefile src/Makefile])
AC_CONFIG_FILES([tests/Makefile examples/Makefile])
AC_CONFIG_FILES([src/gnucap-adms], [chmod +x src/gnucap-adms])
AC_CONFIG_FILES([src/gnucap.xml])
AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   CPPFLAGS:   $CPPFLAGS
   plugindir:  $gnucap_libdir
])
