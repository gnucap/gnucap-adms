# adms rules
# makefile fragment to be included from Makefile.am

ADMSXML=admsXml
GNUCAP_ADMS_IMPLICIT=$(pkgdatadir)/implicit.xml
ADMS_IMPLICIT=$(top_srcdir)/src/implicit.xml

gnucap.xml: gnucap.xml.in

i=0 1 2
ADMS_GC_XML = $(i:%=$(top_builddir)/src/gnucap_%.xml)
ADMS_ES=$(ADMS_GC_XML:%=-e %)
%.cc %.h %_tr.hxx %_top.hxx %_ac.hxx: %.va $(ADMS_GC_XML) $(ADMS_IMPLICIT)
	adms_implicit_transforms=$(ADMS_IMPLICIT) $(ADMSXML) -I$(top_srcdir)/src $< $(ADMS_ES) -o $*

.PRECIOUS: %.cc %.h %.hxx

.va.cc:
ADMSLDFLAGS=-module -shared
