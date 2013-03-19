
GNUCAP_ADMS=../src/gnucap-adms
GNUCAP=LD_LIBRARY_PATH=../src/.libs gnucap

EXFILES=rc.va
if VERILOG_EXTERN
EXT_CHECK=rc ager
rc.test: rc.so
endif
CHECK=opamp sensor_chip callfunctions $(EXT_CHECK)

SENSOR_SO=op.so ldo.so vbg.so sensor.so

EXTRA_DIST=opamp.gc opamp.ref \
           lamp.va lamp.gc \
           ager.gc ager.va ager.ref \
           rc.va rc.gc \
           $(SENSOR_SO:%.so=%.va) \
           $(SENSOR_SO:%.so=%vm.core.gc) \
           sensor_chip.gc sensor_chip.ref sensor_chip.ckt \
           callfunctions.va callfunctions.gc

AGER_PS=ager_tr.2.4.ps \
        ager_tr2.2.4.ps \
        ager_tr_sin.2.4.ps \
        ager_tw_tr.1.2.4.ps

plot: $(AGER_PS) lamp_tr.1.3.ps

pkgexamplesdir=${pkgdatadir}/examples
dist_pkgexamples_DATA=${EXFILES}

i=0 1
ADMS_ES=$(i:%=-e ../src/gnucap_%.xml)
ADMS_GC_XML=../src/gnucap_0.xml

# ugh, also creates some hxx files!
%.cc %.h: %.va $(ADMS_GC_XML)
	admsXml @ADMSXML_FLAGS@ -I../src $< $(ADMS_ES) -o $*

%.so: MCXXFLAGS=--shared -fPIC $(CXXFLAGS)
%.so: MCPPFLAGS=-I ../src $(CPPFLAGS)

#%.so: %.cc
#$(CXX) $(MCXXFLAGS) $(MCPPFLAGS) $(MLDFLAGS) $< -o $@

#%.so: %.va
#	$(GNUCAP_ADMS) -v -c $<

%.so: %.cc %.h
	CPPFLAGS="$(CPPFLAGS) -I../src" LDADD="-L../src/.libs" $(GNUCAP_ADMS) -v $<

# ugh, also creates some hxx files!
#%.cc %.h: %.va
#	$(GNUCAP_ADMS) -v $<

%.ii: %.va
	$(GNUCAP_ADMS) -v -i $<

%.so: %.h %.cc

callfunctions.out: callfunctions.so callfunctions.gc
	$(GNUCAP) callfunctions.gc 2> $@ | tee 2$@ | grep too\ negative ;\
	[ 0 -ne $$? ]

callfunctions.reg: callfunctions.out
	[ 4 -eq $$( grep -c positive 2$< ) ]
	[ 4 -eq $$( grep -c negative $< ) ]
	[ 1 -eq $$( grep -c too\ negative $< ) ]

sensor_so: $(SENSOR_SO)

sensor_chip.test: $(SENSOR_SO)
	
%.reg: %.test %.ref
	diff $+

%.test: %.gc
	$(GNUCAP) $< | grep '^[a-zA-Z0-9]*=' > $@

plot: $(SENSOR_PS)

check: $(CHECK:%=%.reg)

clean: clean-ex

clean-ex:
	rm -f *.so *.hxx *.h *.cc

.PHONY=%.reg clean-ex sensor_so plot

include ../gnuplot.mk
