pkgdata_DATA = gnucap_0.xml gnucap_1.xml gnucap_2.xml

dist_bin_SCRIPTS=gnucap-adms
DEV=cap res opamp
EXTRA_DIST = $(DEV:%=%.va) gnucap.xml split.py \
             discipline.h

gnucap_adms_include_HEADERS = discipline.h gnucap-adms.h e_adms.h
gnucap_adms_includedir = $(pkgincludedir)

plugindir = @gnucap_libdir@
# libdir = @libdir@/gnucap-adms
lib_LTLIBRARIES = libe_adms.la
libe_adms_la_LDFLAGS = -module -shared
libe_adms_la_SOURCES = e_adms.cc

AM_CPPFLAGS=$(GC_CPPFLAGS)

plugin_LTLIBRARIES =

%_0.xml %_1.xml %_2.xml: %.xml
	python @srcdir@/split.py

i=0 1 2
ADMS_GC_XML=$(i:%=gnucap_%.xml)
ADMS_ES=$(i:%=-e gnucap_%.xml)
%.cc %.h %_tr.hxx %_top.hxx %_ac.hxx: %.va $(ADMS_GC_XML)
	admsXml -I$(top_srcdir)/src $< $(ADMS_ES) -o $*

.PRECIOUS: %.cc %.h %.hxx

clean-manually:
	rm -f *.cc *.o *~ *.c gnucap.xml *.hxx *.so gnucap_*.xml

.va.cc:
PLUGIN_LDFLAGS=-module -shared

plugin_LTLIBRARIES += adms_res.la
adms_res_la_LDFLAGS = $(PLUGIN_LDFLAGS)
nodist_adms_res_la_SOURCES = res.va
adms_res_la_LIBADD = -le_adms

plugin_LTLIBRARIES += adms_cap.la
adms_cap_la_LDFLAGS = $(PLUGIN_LDFLAGS)
nodist_adms_cap_la_SOURCES = cap.va
adms_cap_la_LIBADD = -le_adms

plugin_LTLIBRARIES += adms_opamp.la
adms_opamp_la_LDFLAGS = $(PLUGIN_LDFLAGS)
nodist_adms_opamp_la_SOURCES = opamp.va
adms_opamp_la_LIBADD = -le_adms

# automatically by "make depend"?
#$(SOURCES:%.cc=%.o): %.o: %.cc %.h %_tr.hxx %_ac.hxx %_top.hxx

CLEANFILES = .adms* cap.cc res.cc opamp.cc .*.adms *.hxx gnucap_*.xml .*.xml
CLEANFILES+= cap.h res.h opamp.h constants.h
