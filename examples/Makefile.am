# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
GNUCAP = LD_LIBRARY_PATH=../src/.libs $(shell echo gnucap | sed -e '$(transform)')

EXFILES=rc.va
if VERILOG_EXTERN
EXT_TESTS = rc.sh ager.sh
EXT_LA = rc.la ager.la
rc.test: rc.so
endif

SENSOR_SO = op.so ldo.so vbg.so sensor.so

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = \
    op.la vbg.la ldo.la sensor.la callfunctions.la \
    lamp.la motor.la $(EXT_LA)
AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libe_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

EXTRA_DIST= \
    lamp.va ager.va rc.va motor.va \
    $(SENSOR_SO:%.so=%.va) \
    $(SENSOR_SO:%.so=%vm.core.gc) \
    sensor_chip.ckt \
    callfunctions.va callfunctions.gc callfunctions.sh

TEST_EXTENSIONS = .sh
AM_TESTS_FD_REDIRECT = 9>&2
GC_TESTS = opamp.sh sensor_chip.sh lamp.sh motor.sh $(EXT_TESTS)
CUSTOM_TESTS = callfunctions.sh
TESTS = $(GC_TESTS) $(CUSTOM_TESTS)

EXTRA_DIST+= $(GC_TESTS:%.sh=%.gc) $(GC_TESTS:%.sh=%.ref)

$(GC_TESTS): %.sh: %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "exec 2>&9" >> $@
	@echo "diff -rup $< $*.out || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

$(TESTS): %.sh: %.out

# %.out %_tr.out: %.gc
#	$(GNUCAP) $< > $*.out

PREDIFF= \
    '1,/^\#·See·the·file/d' \
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../plugins/.libs/%) $(GNUCAP_ARGS) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

AGER_PS=ager_tr.2.4.ps \
        ager_tr2.2.4.ps \
        ager_tr_sin.2.4.ps \
        ager_tw_tr.1.2.4.ps

plot: $(AGER_PS) lamp_tr.1.3.ps

pkgexamplesdir=${pkgdatadir}/examples
dist_pkgexamples_DATA=${EXFILES}

i=0 1 2
ADMS_ES=$(i:%=-e ../src/gnucap_%.xml)
ADMS_GC_XML=../src/gnucap_0.xml

# ugh, also creates some hxx files!
%.cc %.h: %.va $(ADMS_GC_XML)
	admsXml @ADMSXML_FLAGS@ -I$(top_srcdir)/src $< $(ADMS_ES) -o $*

%.so: MCXXFLAGS=--shared -fPIC $(CXXFLAGS)
%.so: MCPPFLAGS=-I../src $(CPPFLAGS)

#%.so: %.cc
#$(CXX) $(MCXXFLAGS) $(MCPPFLAGS) $(MLDFLAGS) $< -o $@

# %.so: %.cc %.h
#	CPPFLAGS="$(CPPFLAGS) -I../src" LDADD="-L../src/.libs" CXXFLAGS="$(CXXFLAGS)" $(GNUCAP_ADMS) -v $<

# ugh, also creates some hxx files!
#%.cc %.h: %.va
#	$(GNUCAP_ADMS) -v $<

# %.ii: %.va
# 	$(GNUCAP_ADMS) -v -i $<

%.so: %.h %.cc

%.out %.ret 2%.out: %.so %.gc
	$(GNUCAP) $(lastword $+) 2> 2$*.out | tee $*.out | grep too\ negative ;\
	[ 0 -ne $$? ]

lamp.out: lamp.so
motor.out: motor.so
ager.out: ager.so
rc.out: rc.so

callfunctions.log: callfunctions.out 2callfunctions.out

sensor_so: $(SENSOR_SO)

sensor_chip.out: $(SENSOR_SO) $(PWD)/sensor_chip.ckt $(PWD)/opvm.core.gc $(PWD)/vbgvm.core.gc $(PWD)/ldovm.core.gc \
                 vbg.so ldo.so sensor.so $(PWD)/sensorvm.core.gc

$(PWD)/%.ckt: %.ckt
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.gc: %.gc
	[ -f $@ ] || $(LN_S) $< .
	
%.reg: %.ref %.test
	diff $< $*.test

plot: $(SENSOR_PS)

clean-local: clean-ex

clean-ex:
	rm -f *.so *.hxx *.h *.cc

CLEANFILES = .*.adms *.out $(GC_TESTS) *.hxx .*.xml

.PHONY: $(TESTS:%.sh=%.out) clean-ex sensor_so plot
.PRECIOUS: %.cc %.h %.so %.la

@am__include@ @top_srcdir@/gnuplot.mk
