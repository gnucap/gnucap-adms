# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
GNUCAP = $(shell echo gnucap | sed -e '$(transform)')

AM_TESTS_ENVIRONMENT = \
    export PATH='../src:$(PATH)'\
           LD_LIBRARY_PATH='../src/.libs'\
           GNUCAP='$(GNUCAP)'\
           GNUCAP_ADMS_DATADIR='$(PWD)/../src'\
           GNUCAP_ADMS_INCLUDEDIR='$(abs_top_srcdir)/src'\
           GNUCAP_ADMS_IMPLICIT='$(abs_top_srcdir)/src/implicit.xml'\
           GNUCAP_ADMS_CPPFLAGS='-I$(abs_top_srcdir)/src'\
           GNUCAP_ADMS_LDFLAGS='-L$(abs_top_builddir)/src/.libs'\
           GNUCAP_ADMS_LIBS='-lgnucap_adms'\
           GNUCAP_ADMS_TMPDIR='$(PWD)/.admstmp'\
           REDIRECT='exec 2>&9'\
           srcdir='$(srcdir)';

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
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

EXTRA_DIST= \
    lamp.va ager.va rc.va motor.va \
    $(SENSOR_SO:%.so=%.va) \
    $(SENSOR_SO:%.so=%vm.core.gc) \
    sensor_chip.ckt \
    sensor_chip.gc sensor_chip.ref \
    callfunctions.va callfunctions.gc \
    opamp.gc opamp.ref \
    inline.gc inline.ref

sensor_chip.log: $(PWD)/sensor_chip.ckt

TEST_EXTENSIONS = .sh
AM_TESTS_FD_REDIRECT = 9>&2
GC_TESTS = \
    op.sh \
    sensor.sh \
    lamp.sh motor.sh $(EXT_TESTS) \
    vbg.sh
CUSTOM_TESTS = \
    callfunctions.sh \
    sensor_chip.sh \
    inline.sh \
    opamp.sh
TESTS = $(GC_TESTS) $(CUSTOM_TESTS)

inline.log inline.out: cleantmp

cleantmp:
	rm -rf .admstmp

$(GC_TESTS:%.sh=%.log): %.log: %.so

callfunctions.out: callfunctions.so

EXTRA_DIST+= $(GC_TESTS:%.sh=%.gc) $(GC_TESTS:%.sh=%.ref)
FILTER=cat

callfunctions.sh: %.sh: %.out
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "exec 2>&9" >> $@
	@echo "diff -rup $< $*.out || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

sensor_chip.sh inline.sh opamp.sh $(GC_TESTS): %.sh: %.gc Makefile %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-$(GNUCAP)}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP \$$GNUCAP_ARGS $(LT) $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -rup $(lastword $+) - || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

# $(TESTS): %.sh: %.out

# %.out %_tr.out: %.gc
#	$(GNUCAP) $< > $*.out

PREDIFF= \
    '1,/^\#·See·the·file/d' \
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

inline.sh: POSTDIFF='1,/^\#Time/d'

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_TESTS_ENVIRONMENT)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../plugins/.libs/%) $(GNUCAP_ARGS) $< | \
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

callfunctions.log: callfunctions.out 2callfunctions.out
sensor_chip.log: vbg.so op.so ldo.so

sensor_so: $(SENSOR_SO)

vbg.out sensor_chip.log sensor_chip.out: \
    $(SENSOR_SO) $(PWD)/sensor_chip.ckt $(PWD)/opvm.core.gc $(PWD)/vbgvm.core.gc $(PWD)/ldovm.core.gc \
    vbg.so ldo.so sensor.so $(PWD)/sensorvm.core.gc

EXTRA_DIST+= opvm.core.gc vbgvm.core.gc ldovm.core.gc sensorvm.core.gc

$(PWD)/%.ckt: %.ckt
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.gc: %.gc
	[ -f $@ ] || $(LN_S) $< .
	
%.reg: %.ref %.test
	diff $< $*.test

plot: $(SENSOR_PS)

clean-local: clean-ex
	rm -rf .admstmp

clean-ex:
	rm -f *.so *.hxx *.h *.cc

CLEANFILES = .*.adms *.out $(GC_TESTS) *.hxx .*.xml \
    callfunctions.sh \
    inline.sh opamp.sh sensor_chip.sh


$(TESTS:%.sh=%.log): .P

.P:
	@:

.PHONY: .P clean-ex sensor_so plot
.PRECIOUS: %.cc %.h %.so %.la

@am__include@ @top_srcdir@/gnuplot.mk
