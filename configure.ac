#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.67])
AC_INIT([gnucap-adms], [0.0.7-rc2], [felix@salfelder.org])
AC_CONFIG_SRCDIR([README])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE([-Wno-portability silent-rules])
LT_INIT([disable-static])

AS_IF([ test x$program_suffix = xNONE ],,[_program_suffix=$program_suffix])

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
GC_CPPFLAGS+=$(gnucap-conf$_program_suffix --cppflags)

# Checks for typedefs, structures, and compiler characteristics.

dnl ----------------- admsXml ------------------
AC_CHECK_PROG(HAVE_GCCONF, gnucap-conf$_program_suffix, true, false)

if test "x$HAVE_GCCONF" = "xfalse"; then
	AC_MSG_ERROR([cannot find gnucap-conf$_program_suffix, not implemented])
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
	gnucap_libdir=$(gnucap-conf$_program_suffix --libdir)
else
	gnucap_libdir=$HOME/.gnucap/lib
fi

test "x$prefix" = xNONE &&  prefix=$ac_default_prefix
pkglibdir=${prefix}/lib/gnucap-adms
libdir=${prefix}/lib

AC_SUBST([gnucap_libdir])
AC_SUBST([pkglibdir])
AC_SUBST([GC_CPPFLAGS])
CXXFLAGS+=" -Wall -pedantic"
AC_DEFINE_UNQUOTED(ADMS_DATADIR, ["$prefix/share/gnucap-adms"], [path to xml files])
AC_DEFINE_UNQUOTED(ADMS_INCLUDEDIR, ["$prefix/include/gnucap-adms"], [path to include])
AC_DEFINE_UNQUOTED(ADMS_LIBDIR, ["$prefix/lib/gnucap-adms"], [libdir])
AC_DEFINE_UNQUOTED(ADMS_CPPFLAGS, ["-I$prefix/include/gnucap-adms"], [pass to cpp])
AC_DEFINE_UNQUOTED(ADMS_LDFLAGS, ["-L$prefix/lib/gnucap-adms -Wl,-rpath,$prefix/lib/gnucap-adms"], [pass to ld])
AC_DEFINE_UNQUOTED(ADMS_LIBS, ["-lgnucap_adms"], [link against])

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
   pkglibdir:  $pkglibdir
])
