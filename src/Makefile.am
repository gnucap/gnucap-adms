pkgdata_DATA = gnucap_0.xml gnucap_1.xml gnucap_2.xml implicit.xml

dist_bin_SCRIPTS = gnucap-adms
DEV=cap res opamp
EXTRA_DIST = $(DEV:%=%.va) gnucap.xml split.py \
             discipline.h implicit.xml md5.h

include ../adms.mk

gnucap_adms_include_HEADERS = \
    constants.h \
    discipline.h \
    gnucap-adms.h \
    gcuf_compat.h \
    e_adms.h \
    e_admsrc.h \
    m_part.h

gnucap_adms_includedir = $(pkgincludedir)

# libdir = @gnucap_libdir@
admsdir = $(pkglibdir)
adms_LTLIBRARIES = libgnucap_adms.la
libgnucap_adms_la_SOURCES = e_adms.cc

GNUCAP_ADMS_MAKE=$(shell echo gnucap-adms.mk | sed -e '$(transform)')
AM_CPPFLAGS = $(GC_CPPFLAGS)

plugindir = @gnucap_pkglibdir@
plugin_LTLIBRARIES = lang_adms.la
lang_adms_la_LDFLAGS = -module -shared
lang_adms_la_SOURCES = \
	lang_adms.cc md5.c c_attach.cc \
	verilog_include.cc signon.cc
lang_adms_la_CPPFLAGS = $(AM_CPPFLAGS) \
    -DGNUCAP_ADMS_IMPLICIT=\"$(GNUCAP_ADMS_IMPLICIT)\" \
    -DGNUCAP_ADMS_MAKE=\"$(GNUCAP_ADMS_MAKE)\"

%_0.xml %_1.xml %_2.xml: %.xml
	python $(srcdir)/split.py

BUILT_SOURCES = gnucap_0.xml

# hack. somethings wrong with deps.
res.cc cap.cc: gnucap_0.xml

clean-manually:
	rm -f *.cc *.o *~ *.c gnucap.xml *.hxx *.so gnucap_*.xml

# fix the serious ones first.
PLUGIN_CXXFLAGS = -Wno-unused-variable -Wno-unused-but-set-variable

plugin_LTLIBRARIES += d_adms_vbr.la
d_adms_vbr_la_LDFLAGS = $(ADMSLDFLAGS)
d_adms_vbr_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
d_adms_vbr_la_SOURCES = d_adms_vbr.cc
d_adms_vbr_la_CPPFLAGS = $(AM_CPPFLAGS) -DADD_VERSION

plugin_LTLIBRARIES += d_adms_i.la
d_adms_i_la_LDFLAGS = $(ADMSLDFLAGS)
d_adms_i_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
d_adms_i_la_SOURCES = d_adms_i.cc
d_adms_i_la_CPPFLAGS = $(AM_CPPFLAGS) -DADD_VERSION

plugin_LTLIBRARIES += adms_res.la
adms_res_la_LDFLAGS = $(ADMSLDFLAGS)
adms_res_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_res_la_SOURCES = res.va
adms_res_la_LIBADD = libgnucap_adms.la

plugin_LTLIBRARIES += adms_cap.la
adms_cap_la_LDFLAGS = $(ADMSLDFLAGS)
adms_cap_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_cap_la_SOURCES = cap.va
adms_cap_la_LIBADD = libgnucap_adms.la

plugin_LTLIBRARIES += adms_opamp.la
adms_opamp_la_LDFLAGS = $(ADMSLDFLAGS)
adms_opamp_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_opamp_la_SOURCES = opamp.va
adms_opamp_la_LIBADD = libgnucap_adms.la

if C_DEBUG
debug_plugin =
else
debug_plugin = c_debug.la
endif
c_debug_la_SOURCES = c_debug.cc
c_debug_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
c_debug_la_LDFLAGS = -shared -module

plugin_LTLIBRARIES += ${debug_plugin}

# automatically by "make depend"?
#$(SOURCES:%.cc=%.o): %.o: %.cc %.h %_tr.hxx %_ac.hxx %_top.hxx

CLEANFILES = .adms* cap.cc res.cc opamp.cc .*.adms *.hxx gnucap_*.xml .*.xml
CLEANFILES+= cap.h res.h opamp.h constants.h

# stated in adms.mk. does not work?
%.cc %.h: %.va $(top_srcdir)/src/implicit.xml gnucap_0.xml
