#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([gnucap-adms], [0.0.8-rc2], [felix@salfelder.org])
AC_CONFIG_SRCDIR([README])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wno-portability silent-rules])
LT_INIT([disable-static])

AS_IF([ test x$program_suffix = xNONE ],,[_program_suffix=$program_suffix])
AM_CONDITIONAL([VPATH_BUILD], [test ! x$srcdir = "x."])

dnl ------------------- gnucap-uf? -------------------

AC_MSG_CHECKING([gnucap-uf availibility])
echo | gnucap$_program_suffix -a "lazy models-uf.so" -b /dev/stdin > /dev/null 2>/dev/null
MODELS_UF=${PIPESTATUS[[1]]}
AS_IF([ test $MODELS_UF -eq 0 ],[
	GNUCAP_UF=yes
		 ],[
	GNUCAP_UF=no
			 ])
AC_MSG_RESULT([$GNUCAP_UF])

dnl this is required for testing.
AM_CONDITIONAL([GNUCAP_UF], [test x"$GNUCAP_UF" != xno])

dnl ------------------- debug? -------------------

AC_MSG_CHECKING([gnucap debug plugin availibility])
echo | gnucap$_program_suffix -a "lazy c_debug.so" -b /dev/stdin > /dev/null 2>/dev/null
C_DEBUG=${PIPESTATUS[[1]]}
AS_IF([ test $C_DEBUG -eq 0 ],[
	C_DEBUG=yes
		 ],[
	C_DEBUG=no
			 ])
AC_MSG_RESULT([$C_DEBUG])

dnl this is required for testing.
AM_CONDITIONAL([C_DEBUG], [test x"$C_DEBUG" != xno])

dnl ------------------- debugging -------------------
AC_MSG_CHECKING([if debug code should be compiled in])
AC_ARG_ENABLE([debug],
[  --enable-debug          Enable building of debug code. [[default: disabled]]],
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
GC_CPPFLAGS+=$(gnucap-conf$_program_suffix --cppflags)

# Checks for typedefs, structures, and compiler characteristics.

dnl ----------------- admsXml ------------------
AC_CHECK_PROG(HAVE_GCCONF, gnucap-conf$_program_suffix, true, false)

if test "x$HAVE_GCCONF" = "xfalse"; then
	AC_MSG_ERROR([cannot find gnucap-conf$_program_suffix, not implemented])
fi

dnl ----------------- admsXml ------------------
AC_CHECK_PROG(HAVE_ADMSXML, admsXml, true, false)
if test "x$HAVE_ADMSXML" = "xfalse"; then
	AC_MSG_ERROR([cannot find admsXml])
fi

AC_MSG_CHECKING([extern support])
# FIXME: check gnucap capabilities...
AC_SUBST([ADMSXML_FLAGS])

dnl ----------------- qucs ------------------
AC_MSG_CHECKING([qucs])
AC_ARG_ENABLE([qucs],
[  --enable-qucs           try to build qucs models (experimental). [[default: disabled]]],
[
if test "X$enable_qucs" = "Xno" ; then
	AC_MSG_RESULT([no])
else
	AC_MSG_RESULT([yes])
	enable_qucs=yes
fi
],
[
	AC_MSG_RESULT([no])
	enable_qucs=no
])
AM_CONDITIONAL([QUCS], [test x"$enable_qucs" != xno])
dnl ----------------- plugindir ------------------
AC_MSG_CHECKING([local])
AC_ARG_ENABLE([local],
[  --enable-local          install plugins locally. [[default: disabled]]],
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
	gnucap_libdir=$(gnucap-conf$_program_suffix --libdir)
	gnucap_pkglibdir=$(gnucap-conf$_program_suffix --pkglibdir)
else
	# has this ever worked?
	gnucap_libdir=$HOME/.gnucap/lib
	gnucap_pkglibdir=$HOME/.gnucap/lib
fi

test "x$prefix" = xNONE &&  prefix=$ac_default_prefix
pkglibdir=${prefix}/lib/gnucap-adms
AC_SUBST([pkglibdir])

libdir=${prefix}/lib
gnucap_exec_prefix="$(gnucap-conf$_program_suffix --exec-prefix)"
GNUCAP_MAKE="$gnucap_exec_prefix/bin/gnucap.mk$_program_suffix"
AC_SUBST([GNUCAP_MAKE])

AC_SUBST([gnucap_libdir])
AC_SUBST([gnucap_pkglibdir])
AC_SUBST([GC_CPPFLAGS])
CXXFLAGS+=" -Wall -pedantic"
AC_DEFINE_UNQUOTED(ADMS_DATADIR, ["$prefix/share/gnucap-adms"], [path to xml files])
AC_DEFINE_UNQUOTED(ADMS_INCLUDEDIR, ["$prefix/include/gnucap-adms"], [path to include])
AC_DEFINE_UNQUOTED(ADMS_LIBDIR, ["$prefix/lib/gnucap-adms"], [libdir])
AC_DEFINE_UNQUOTED(ADMS_CPPFLAGS, ["-I$prefix/include/gnucap-adms"], [pass to cpp])
AC_DEFINE_UNQUOTED(ADMS_LDFLAGS, ["-L$prefix/lib/gnucap-adms -Wl,-rpath,$prefix/lib/gnucap-adms"], [pass to ld])
AC_DEFINE_UNQUOTED(ADMS_LIBS, ["-lgnucap_adms"], [link against])

WHICH_MAKE=`which make`
AC_SUBST([WHICH_MAKE])

AC_SUBST([nothing])
# Checks for library functions.

if test "$enable_debug" = "yes" ; then
	CPPFLAGS="$CPPFLAGS -DTRACE_UNTESTED"
else
	CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([qucs/Makefile])
AC_CONFIG_FILES([tests/Makefile])
AC_CONFIG_FILES([src/gnucap-adms], [chmod +x src/gnucap-adms])
AC_CONFIG_FILES([gnucap-adms.mk], [chmod +x gnucap-adms.mk])
AC_CONFIG_FILES([src/gnucap.xml])
AC_OUTPUT

AC_MSG_RESULT([
** Configuration summary for $PACKAGE $VERSION:

   prefix:     $prefix
   CPPFLAGS:   $CPPFLAGS
   plugindir:  $gnucap_libdir
   pkglibdir:  $pkglibdir
   uf:         $GNUCAP_UF
   qucs:       $enable_qucs
])
