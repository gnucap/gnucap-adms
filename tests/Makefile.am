
GNUCAP_ADMS=../src/gnucap-adms
GNUCAP=LD_LIBRARY_PATH=../src/.libs gnucap

EXFILES=if.va
if VERILOG_EXTERN
EXT_CHECK=if
if.test: if.so
endif
CHECK=$(EXT_CHECK)

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

%.cc %.h: %.va
	$(GNUCAP_ADMS) -v $<

%.ii: %.va
	$(GNUCAP_ADMS) -v -i $<

%.so: %_tr.hxx %_top.hxx %_ac.hxx %.cc

%.reg: %.test %.ref
	diff $+

%.test: %.gc
	$(GNUCAP) $< | grep '^[a-zA-Z0-9]*=' > $@

check: $(CHECK:%=%.reg)

clean: clean-examples

clean-examples:
	rm -f *.out *.so *.cc *.hxx *.h

.PHONY: %.reg clean-examples
