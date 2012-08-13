pkgdata_DATA = gnucap_0.xml gnucap_1.xml gnucap_2.xml

dist_bin_SCRIPTS=gnucap-adms

gnucap_adms_include_HEADERS = discipline.h gnucap-adms.h e_adms.h
gnucap_adms_includedir = $(pkgincludedir)

plugindir = @gnucap_libdir@
plugin_LTLIBRARIES = e_adms.la
e_adms_la_LDFLAGS = -module -shared
e_adms_la_SOURCES = e_adms.cc

%_0.xml %_1.xml %_2.xml: %.xml
	python split.py

i=0 1 2
ADMS_GC_XML=$(i:%=gnucap_%.xml)
ADMS_ES=$(i:%=-e $(srcdir)/gnucap_%.xml)
%.cc %.h %_tr.hxx %_top.hxx %_ac.hxx: %.va $(ADMS_GC_XML)
	admsXml $< $(ADMS_ES) -o $*

.PPRECIOUS: %.cc %.h %.hxx

clean:
	rm -f *.cc *.h *.o *~ *.c gnucap.xml *.hxx *.so gnucap_*.xml

plugin_LTLIBRARIES += adms_res.la
adms_res_la_LDFLAGS = -module -shared
adms_res_la_SOURCES = res.cc

plugin_LTLIBRARIES += adms_cap.la
adms_cap_la_LDFLAGS = -module -shared
adms_cap_la_SOURCES = cap.cc

# automatically by "make depend"?
#$(SOURCES:%.cc=%.o): %.o: %.cc %.h %_tr.hxx %_ac.hxx %_top.hxx
