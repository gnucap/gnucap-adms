#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([gnucap-adms], [0.0.1], [felix@salfelder.org])
AC_CONFIG_SRCDIR([README])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE
LT_INIT()

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
if admsXml src/gnucap-adms.h 2>/dev/null >/dev/null; then
	AC_MSG_RESULT([yes])
	AC_DEFINE([VERILOG_EXTERN], [1], [admsXml understands 'extern'])
else
	AC_MSG_RESULT([no])
fi

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

AC_CONFIG_FILES([Makefile src/Makefile examples/Makefile])
AC_CONFIG_FILES([src/gnucap-adms], [chmod +x src/gnucap-adms])
AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   CPPFLAGS:   $CPPFLAGS
   plugindir:  $gnucap_libdir
])
