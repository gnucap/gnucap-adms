plugindir = @gnucap_pkglibdir@
plugin_LTLIBRARIES =

include ../adms.mk
AM_CPPFLAGS = $(GC_CPPFLAGS) -I $(top_srcdir)/src

# fix the serious ones first.
PLUGIN_CXXFLAGS = -Wno-unused-variable -Wno-unused-but-set-variable
lgadms=../src/libgnucap_adms.la

plugin_LTLIBRARIES += adms_res.la
adms_res_la_LDFLAGS = $(ADMSLDFLAGS)
adms_res_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_res_la_SOURCES = res.va
adms_res_la_LIBADD = $(lgadms)

plugin_LTLIBRARIES += adms_cap.la
adms_cap_la_LDFLAGS = $(ADMSLDFLAGS)
adms_cap_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_cap_la_SOURCES = cap.va
adms_cap_la_LIBADD = $(lgadms)

plugin_LTLIBRARIES += adms_opamp.la
adms_opamp_la_LDFLAGS = $(ADMSLDFLAGS)
adms_opamp_la_CXXFLAGS = $(PLUGIN_CXXFLAGS)
nodist_adms_opamp_la_SOURCES = opamp.va
adms_opamp_la_LIBADD = $(lgadms)
