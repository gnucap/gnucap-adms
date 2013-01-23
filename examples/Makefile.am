
GNUCAP_ADMS=../src/gnucap-adms
GNUCAP=gnucap

EXFILES=rc.va
CHECK=opamp

EXTRA_DIST=opamp.gc opamp.ref lamp.va lamp.gc

pkgexamplesdir=${pkgdatadir}/examples
dist_pkgexamples_DATA=${EXFILES}

i=0 1 2
ADMS_ES=$(i:%=-e ../adms/gnucap_%.xml)
ADMS_GC_XML=../adms/gnucap_0.xml

# before installation
# %.cc %_tr.hxx %_top.hxx %_ac.hxx: %.va $(ADMS_GC_XML)
# 	admsXml -I../adms $< $(ADMS_ES) -o $*

%.so: MCXXFLAGS=--shared -fPIC $(CXXFLAGS)
%.so: MCPPFLAGS=-I ../src $(CPPFLAGS)

#%.so: %.cc
#$(CXX) $(MCXXFLAGS) $(MCPPFLAGS) $(MLDFLAGS) $< -o $@

%.so: %.va
	$(GNUCAP_ADMS) -v -c $<

%.so: %.cc
	$(GNUCAP_ADMS) -v $<

%.cc %_tr.hxx %_top.hxx %_ac.hxx: %.va
	$(GNUCAP_ADMS) -v $<

%.ii: %.va
	$(GNUCAP_ADMS) -v -i $<

%.so: %_tr.hxx %_top.hxx %_ac.hxx %.cc

%.reg: %.test %.ref
	diff $+

%.test: %.gc
	$(GNUCAP) $< | grep '^[a-zA-Z0-9]*=' > $@

check: $(CHECK:%=%.reg)

.PHONY=%.reg
