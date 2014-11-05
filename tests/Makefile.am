# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
# GNUCAP = LD_LIBRARY_PATH=../src/.libs
GNUCAP = $(shell echo gnucap | sed -e '$(transform)')
plugindir = @gnucap_libdir@
FILTER = $(plugindir)/gnucap-testfilter

sensor_chip.log: $(PWD)/sensor_chip.ckt

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = \
    $(EXT_LIBS) \
    if.la if2.la \
    vs.la cs.la cap.la \
    coil.la ccvs.la vccs.la \
    analogfunction.la \
    afsquare.la \
    square.la \
    int.la pid.la vcvs.la \
    diff.la \
    mul.la


AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .
$(PWD)/%.ckt: %.ckt
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.gc: %.gc
	[ -f $@ ] || $(LN_S) $< .
	

DIST_TESTS = $(SO_TESTS) $(EXT_TESTS)

if VERILOG_EXTERN
EXT_LA = rc.la ager.la motor.la
ager.log rc.log : %.log: %.so
EXT_TESTS = \
    subcoil.sh \
    rc.sh \
    ager.sh \
    ttint.sh
EXT_LIBS = \
    subcoil.la
endif


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

AM_TESTS_FD_REDIRECT = 9>&2
TEST_EXTENSIONS = .sh
TESTS = $(EXT_TESTS) $(GC_TESTS) $(GCUF_TESTS) $(INLINE_TESTS)
TESTS+= $(CUSTOM_TESTS)
INLINE_TESTS = abs.sh instanceparam.sh
GC_TESTS = \
    if.sh \
    if2.sh \
    coil.sh \
    vs.sh \
    cs.sh cap.sh \
    analogfunction.sh \
    afsquare.sh \
    square.sh \
    int.sh \
    pid.sh \
    diff.sh \
    mul.sh \
    ccvs.sh \
    vccs.sh \
    vcvs.sh
# from examples
GC_TESTS+ = $(SO_TESTS)

# will only work with -uf
GCUF_TESTS = rc_uic.sh ddc_rc.sh ddc.2.sh

POSTDIFF=
# remove cruft from model compiler
$(empty)abs.sh: POSTDIFF='1,/^testdata/d'
$(empty)instanceparam.sh: POSTDIFF='1,/^testdata/d'
$(empty)ttint.sh: POSTDIFF='0,/^testdata/d'
$(empty)inline.sh: POSTDIFF='0,/^outputstart/d'


EXTRA_DIST = \
    bug1.va \
    $(GC_TESTS:%.sh=%.gc) \
    $(GC_TESTS:%.sh=%.ref) \
    $(GCUF_TESTS:%.sh=%.gc) \
    $(GCUF_TESTS:%.sh=%.ref) \
    $(INLINE_TESTS:%.sh=%.gc) \
    $(INLINE_TESTS:%.sh=%.ref) \
    rc_uic.gc ddc_rc.gc \
    subcoil.gc subcoil.ref \
    lamp.va ager.va rc.va motor.va \
    $(SENSOR_SO:%.so=%.va) \
    $(SENSOR_SO:%.so=%vm.core.gc) \
    sensor_chip.ckt \
    sensor_chip.gc sensor_chip.ref \
    callfunctions.va callfunctions.gc \
    opamp.gc opamp.ref \
    inline.gc inline.ref
EXTRA_DIST+= $(DIST_TESTS:%.sh=%.gc) $(DIST_TESTS:%.sh=%.ref)
EXTRA_DIST+= \
    opvm.core.gc \
    vbgvm.core.gc \
    ldovm.core.gc \
    sensorvm.core.gc \
    ldo.gc \
    ldo.ref \
    motor.gc \
    motor.ref \
    net.va

SO_TESTS = \
    op.sh \
    sensor.sh \
    lamp.sh \
    vbg.sh
EXT_SO_TESTS = \
    motor.sh
CUSTOM_TESTS = \
    callfunctions.sh \
    sensor_chip.sh \
    inline.sh \
    opamp.sh

vbg.out sensor_chip.log sensor_chip.out: \
    $(SENSOR_SO) $(PWD)/sensor_chip.ckt $(PWD)/opvm.core.gc $(PWD)/vbgvm.core.gc $(PWD)/ldovm.core.gc \
    vbg.so ldo.so sensor.so $(PWD)/sensorvm.core.gc

callfunctions.out: callfunctions.so

callfunctions.sh: %.sh: %.out
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "exec 2>&9" >> $@
	@echo "diff -rup $< $*.out || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

SENSOR_SO = op.so ldo.so vbg.so sensor.so

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES += \
    op.la vbg.la ldo.la sensor.la callfunctions.la \
    lamp.la $(EXT_LA)
AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

LT = <

SHS = sensor_chip.sh inline.sh opamp.sh $(INLINE_TESTS) $(GCUF_TESTS) $(GC_TESTS) $(EXT_TESTS)
$(SHS): %.sh: %.gc Makefile %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-$(GNUCAP)}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP \$$GNUCAP_ARGS $(LT) $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -rup $(lastword $+) - || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

subcoil.log vccs.log ccvs.log cap.log: %.log: %.so
coil.log: coil.so
vs.log: vs.so
cs.log: cs.so
if.log: if.so
if2.log: if2.so

$(GC_TESTS:%.sh=%.log): %.log: %.so
$(GC_TESTS:%.sh=%.out): %.out: %.so

PREDIFF= \
    '1,/^\#·See·the·file/d' \
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

$(TESTS:%.sh=%.log): Makefile

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(AM_TESTS_ENVIRONMENT)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../plugins/.libs/%) $(GNUCAP_ARGS) $(LT) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

i=0 1 2
ADMS_ES=$(i:%=-e ../src/gnucap_%.xml)
ADMS_GC_XML=../src/gnucap_0.xml

# ugh, also creates some hxx files!
%.cc %.h: %.va $(ADMS_GC_XML)
	adms_implicit_transforms=$(top_srcdir)/src/implicit.xml admsXml @ADMSXML_FLAGS@ -I$(top_srcdir)/src $< $(ADMS_ES) -o $(notdir $*)

%.cc %.h: %.va
	$(GNUCAP_ADMS) -v $<

%.ii: %.va
	$(GNUCAP_ADMS) -v -i $<

%.so: %_tr.hxx %_top.hxx %_ac.hxx %.cc

# %.reg: %.test %.ref
# 	diff $+

# %.test: %.gc
# 	$(GNUCAP) $< | grep '^[a-zA-Z0-9]*=' > $@

check: $(CHECK:%=%.reg)

clean-local: clean-examples cleantmp

clean-examples:
	rm -f *.out *.so *.cc *.hxx *.h

abs.log abs.out: cleantmp
instanceparam.log instanceparam.out: cleantmp

cleantmp:
	rm -rf .admstmp

CLEANFILES = .*.adms *.out \
    $(GC_TESTS) \
    $(EXT_TESTS) \
    $(INLINE_TESTS) \
    *.hxx .*.xml \
    core \
    rc_uic.sh \
    ddc_rc.sh \
    ddc.2.sh

CLEANFILES+= $(check_LTLIBRARIES:%.la=%.h) \
             $(check_LTLIBRARIES:%.la=%.cc)
CLEANFILES+= .*.adms *.out *.hxx .*.xml \
    callfunctions.sh \
    inline.sh opamp.sh sensor_chip.sh


$(TESTS:%.sh=%.log): .P
$(TESTS:%.sh=%.out): .P

.P:
	@:

.PHONY: %.reg clean-examples .P

@am__include@ @top_srcdir@/gnuplot.mk
