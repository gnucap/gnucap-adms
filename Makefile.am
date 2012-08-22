pkgdata_DATA = gnucap_0.xml gnucap_1.xml gnucap_2.xml

dist_bin_SCRIPTS=gnucap-adms

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
	python split.py

i=0 1 2
ADMS_GC_XML=$(i:%=gnucap_%.xml)
ADMS_ES=$(i:%=-e $(srcdir)/gnucap_%.xml)
%.cc %.h %_tr.hxx %_top.hxx %_ac.hxx: %.va $(ADMS_GC_XML)
	admsXml $< $(ADMS_ES) -o $*

.PPRECIOUS: %.cc %.h %.hxx

clean-manually:
	rm -f *.cc *.o *~ *.c gnucap.xml *.hxx *.so gnucap_*.xml

PLUGIN_LDFLAGS=-module -shared

plugin_LTLIBRARIES += adms_res.la
adms_res_la_LDFLAGS = $(PLUGIN_LDFLAGS)
adms_res_la_SOURCES = res.cc
adms_res_la_LIBADD = -le_adms

plugin_LTLIBRARIES += adms_cap.la
adms_cap_la_LDFLAGS = $(PLUGIN_LDFLAGS)
adms_cap_la_SOURCES = cap.cc
adms_cap_la_LIBADD = -le_adms

# automatically by "make depend"?
#$(SOURCES:%.cc=%.o): %.o: %.cc %.h %_tr.hxx %_ac.hxx %_top.hxx
