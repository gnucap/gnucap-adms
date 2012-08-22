##
# this makefile is intended to work _after_ installation

EXFILES=rc.va

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
	gnucap-adms -v -c $<

%.so: %.cc
	gnucap-adms -v $<

%.cc %_tr.hxx %_top.hxx %_ac.hxx: %.va
	gnucap-adms -v $<

%.ii: %.va
	gnucap-adms -v -i $<

%.so: %_tr.hxx %_top.hxx %_ac.hxx %.cc
