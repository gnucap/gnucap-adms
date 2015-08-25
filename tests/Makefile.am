# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
# GNUCAP = LD_LIBRARY_PATH=../src/.libs
GNUCAP = $(shell echo gnucap | sed -e '$(transform)')
plugindir = @gnucap_libdir@
FILTER = $(plugindir)/gnucap-testfilter
GNUCAP_PLUGINS = d_adms_i.so d_adms_vbr.so

sensor_chip.log: $(PWD)/sensor_chip.ckt
aging_mos.log: $(PWD)/aging_mos.va
diff.2.log: $(PWD)/diff.va

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = \
    subcoil.la \
    subcoil.1.la \
    if.la if2.la \
    chainvs.la \
    vs.la \
    vs1.la \
    vs2.la \
    cs.la cap.la \
    coil.la \
    ccvs.la \
    cccs.la \
    res.la \
    short.la \
    collapse.la \
    vccs.la \
    analogfunction.la \
    afsquare.la \
    square.la \
    int.la pid.la vcvs.la \
    diff.la \
    mul.la \
    rc.la \
    rc_emb.la \
    ager.la \
    motor.la \
    vbg.la ldo.la sensor.la callfunctions.la \
    op.la \
    lamp.la

## if AGING
check_LTLIBRARIES += \
    aging_mos.la
## endif


AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
AM_CXXFLAGS = -Wno-unused-variable -Wno-unused-but-set-variable

LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || [ -L $@ ] || $(LN_S) .libs/$@ .
$(PWD)/%.ckt: %.ckt
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.gc: %.gc
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.va: %.va
	[ -f $@ ] || $(LN_S) $< .

DIST_TESTS = $(SO_TESTS)

ager.log rc.log : %.log: %.so

AM_TESTS_ENVIRONMENT = \
    export PATH='../src:$(PATH)'\
           LD_LIBRARY_PATH='../src/.libs'\
           GNUCAP='$(GNUCAP)'\
           GNUCAP_ADMS_CPPFLAGS='-I$(abs_top_srcdir)/src'\
           GNUCAP_ADMS_DATADIR='$(PWD)/../src'\
           GNUCAP_ADMS_INCLUDEDIR='$(abs_top_srcdir)/src'\
           GNUCAP_ADMS_IMPLICIT='$(abs_top_srcdir)/src/implicit.xml'\
           GNUCAP_ADMS_LDFLAGS='-L$(abs_top_builddir)/src/.libs'\
           GNUCAP_ADMS_LIBDIR='$(abs_top_builddir)/src/.libs'\
           GNUCAP_ADMS_LIBS='-lgnucap_adms'\
           GNUCAP_ADMS_MAKE='$(abs_top_builddir)/gnucap-adms.mk'\
           GNUCAP_ADMS_TMPDIR='$(PWD)/.admstmp'\
           GNUCAP_ADMS_XMLDIR='$(abs_top_builddir)/src'\
           REDIRECT='exec 2>&9'\
           srcdir='$(srcdir)';

AM_TESTS_FD_REDIRECT = 9>&2
TEST_EXTENSIONS = .sh
TESTS = $(GC_TESTS) $(GCUF_TESTS) $(INLINE_TESTS)
TESTS+= $(CUSTOM_TESTS)
INLINE_TESTS = abs.sh instanceparam.sh
# tests that require an .so
SO_TESTS = \
    op.sh \
    sensor.sh \
    lamp.sh \
    vbg.sh \
    subcoil.sh \
    subcoil.1.sh \
    rc.sh \
    short.sh \
    collapse.sh \
    rc_emb.sh \
    motor.sh
CUSTOM_TESTS = \
    callfunctions.sh \
    coil_auto.sh \
    coil_br.sh \
    diff.2.sh \
    sensor_chip.sh \
    vsstamp.sh \
    csstamp.sh \
    vcvsstamp.sh \
    ccvsstamp.sh \
    cccsstamp.sh \
    inline.sh \
    rc1.sh \
    tr.sh \
    int_cont.sh \
    motor1.sh \
    adms_res.sh \
    opamp.sh
GC_TESTS = \
    $(SO_TESTS) \
    if.sh \
    if2.sh \
    coil.sh \
    chainvs.sh \
    vs.sh \
    vs1.sh \
    vs2.sh \
    cs.sh cap.sh \
    analogfunction.sh \
    afsquare.sh \
    square.sh \
    int.sh \
    ldo.sh \
    pid.sh \
    diff.sh \
    mul.sh \
    cccs.sh \
    ccvs.sh \
    res.sh \
    vccs.sh \
    vcvs.sh
# will only work with -uf
GCUF_TESTS = \
    bvs.sh \
    ager.sh \
	 aging_mos.sh \
    rc_uic.sh \
    ddc_rc.sh \
    ddc.2.sh \
    ttint.sh

XFAIL_TESTS =

# not implemented yet
XFAIL_TESTS += collapse.sh

# not implemented yet
XFAIL_TESTS += diff.2.sh

# does not work.
XFAIL_TESTS += coil_br.sh

POSTDIFF=
# remove cruft from model compiler
$(empty)abs.sh: POSTDIFF='1,/^testdata/d'
$(empty)instanceparam.sh: POSTDIFF='1,/^testdata/d'
$(empty)ttint.sh: POSTDIFF='0,/^testdata/d'
$(empty)inline.sh: POSTDIFF='0,/^outputstart/d'
$(empty)aging_mos.sh: POSTDIFF='0,/^done·load_va/d'

# FIXME: more automatically...?
EXTRA_DIST = \
    aging_mos.va \
    bug1.va \
    $(GC_TESTS:%.sh=%.gc) \
    $(GC_TESTS:%.sh=%.ref) \
    $(GCUF_TESTS:%.sh=%.gc) \
    $(GCUF_TESTS:%.sh=%.ref) \
    $(INLINE_TESTS:%.sh=%.gc) \
    $(INLINE_TESTS:%.sh=%.ref) \
    rc_uic.gc ddc_rc.gc \
    lamp.va ager.va rc.va motor.va \
    $(SENSOR_SO:%.so=%.va) \
    $(SENSOR_SO:%.so=%vm.core.gc) \
    sensor_chip.ckt \
    callfunctions.va \
    $(DIST_TESTS:%.sh=%.gc) $(DIST_TESTS:%.sh=%.ref) \
    $(CUSTOM_TESTS:%.sh=%.gc) $(CUSTOM_TESTS:%.sh=%.ref) \
    opvm.core.gc \
    vbgvm.core.gc \
    ldovm.core.gc \
    sensorvm.core.gc \
    ldo.gc \
    ldo.ref \
    motor.gc \
    motor.ref \
    net.va

vbg.out sensor_chip.log sensor_chip.out: \
    $(SENSOR_SO) $(PWD)/sensor_chip.ckt $(PWD)/opvm.core.gc $(PWD)/vbgvm.core.gc $(PWD)/ldovm.core.gc \
    vbg.so ldo.so sensor.so $(PWD)/sensorvm.core.gc

%.out %.ret 2%.out: %.gc
	$(GNUCAP) < $(lastword $+) 2> 2$*.out | tee $*.out | grep too\ negative ;\
	[ 0 -ne $$? ]

callfunctions.out: callfunctions.so

SENSOR_SO = op.so ldo.so vbg.so sensor.so


LT = <

callfunctions.sh: LT=2> >(tee /dev/fd/2 | grep negative | while read line; do echo "[ERROR] \$$line"; done ) <

SHS = \
	$(CUSTOM_TESTS) \
   $(INLINE_TESTS) $(GCUF_TESTS) $(GC_TESTS)

$(SHS): %.sh: %.gc Makefile %.ref
	@echo "#!$(SHELL)" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-$(GNUCAP)}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP -a ../src/.libs/d_adms_i.so -a ../src/.libs/d_adms_vbr.so \$$GNUCAP_ARGS $(LT) $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -rup $(lastword $+) - || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

subcoil.log subcoil.1.log vccs.log ccvs.log cap.log: %.log: %.so
coil.log: coil.so
coil_auto.log: coil.so
coil_br.log: coil.so
vs.log: vs.so
cs.log: cs.so
if.log: if.so
if2.log: if2.so
motor1.log: motor.so
callfunctions.log: callfunctions.so

$(GC_TESTS:%.sh=%.log): %.log: %.so
$(GC_TESTS:%.sh=%.out): %.out: %.so

PREDIFF= \
    '1,/^core-lib·version/d' \
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

$(TESTS:%.sh=%.log): Makefile

check_custom: $(CUSTOM_TESTS:%.sh=%.log)

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(AM_TESTS_ENVIRONMENT)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../src/.libs/%) $(GNUCAP_ARGS) $(LT) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

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
inline.log inline.out: cleantmp

cleantmp:
	rm -rf .admstmp

CLEANFILES = .*.adms *.out \
    $(GC_TESTS) \
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

include ../adms.mk

@am__include@ @top_srcdir@/gnuplot.mk
