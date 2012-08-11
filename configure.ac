#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AC_CONFIG_SRCDIR([discipline.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE(gnucap-adms, 2012-08-11)
LT_INIT()

# Checks for programs.
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_CXX

# Checks for libraries.

# Checks for header files.
CPPFLAGS+=$(gnucap-conf --cppflags)

# Checks for typedefs, structures, and compiler characteristics.

gnucap_libdir=$(gnucap-conf --libdir)
AC_SUBST([gnucap_libdir])

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([gnucap-adms], [chmod +x gnucap-adms])
AC_OUTPUT
