# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
# GNUCAP = LD_LIBRARY_PATH=../src/.libs
GNUCAP = $(shell echo gnucap | sed -e '$(transform)')
plugindir = @gnucap_pkglibdir@
if GNUCAP_UF
ALT_REF_DIR=uf
FILTER = $(plugindir)/gnucap-testfilter
else
FILTER = cat
endif
GNUCAP_PLUGINS = d_adms_i.so d_adms_vbr.so

sensor_chip.log: $(PWD)/sensor_chip.ckt
aging_mos.log: $(PWD)/aging_mos.va
diff.2.log: $(PWD)/diff.va

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = \
    $(UFCHECKLIB) \
    forloop.la \
    BSIM6.1.1.la \
    case.la \
    ddx.la \
    subcoil.la \
    subcoil.1.la \
    if.la \
    if2.la \
    chainvs.la \
    vs.la \
    vs1.la \
    vs2.la \
    cap.la \
    cap_rescale.la \
    chemical.la \
    cs.la \
    coil.la \
    ccvs.la \
    cccs.la \
    implicit.la \
    mid.la \
    res.la \
    res2.la \
    short.la \
    collapse.la \
    vccs.la \
    af.0.la \
    af.1.la \
    af.2.la \
    afbuf.la \
    afsquare.la \
    square.la \
    int.la pid.la \
	 vcvs.la \
    diff.la \
    hyst1.la \
    meas.la \
    limit.la \
    log.la \
    mul.la \
    rc.la \
    rc_emb.la \
    sin.la \
    motor.la \
    vbg.la ldo.la sensor.la callfunctions.la \
    op.la \
    lamp.la

if GNUCAP_UF
UFCHECKLIB = \
    ager.la \
    ager.sub.la \
    aging_mos.la \
	 d_ageQ_18.la \
	 mos1_a7.la \
    mos.1.la
else
	UFCHECKLIB =
endif

BSIM6_1_1_lo_CXXFLAGS =
@nothing@BSIM6_1_1_la-BSIM6.1.1.lo: CXXFLAGS = -O0 # -Wall -pedantic
@nothing@BSIM6.1.1.lo: CXXFLAGS = -O0 # -Wall -pedantic

# BUG. assertions fail in BSIM6 (uninitialized vars)
@nothing@BSIM6_1_1_la-BSIM6.1.1.lo: CPPFLAGS += -DNDEBUG
@nothing@BSIM6.1.1.lo: CPPFLAGS += -DNDEBUG

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
$(PWD)/%.v: %.v
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/%.va: %.va
	[ -f $@ ] || $(LN_S) $< .
$(PWD)/modelcard.%: modelcard.%
	[ -f $@ ] || $(LN_S) $< .

DIST_TESTS = $(SO_TESTS)

rc.log : %.log: %.so

EXLDL=$(abs_top_builddir)/examples/.libs
SRCLDL=$(abs_top_builddir)/src/.libs

AM_TESTS_ENVIRONMENT = \
    export PATH='../src:$(PATH)'\
           LD_LIBRARY_PATH='$(abs_builddir)/.libs:$(SRCLDL):$(EXLDL)${LD_LIBRARY_PATH:%=:%}'\
           ALT_REF_DIR='$(ALT_REF_DIR)'\
           GNUCAP='$(GNUCAP)'\
           GNUCAP_TESTFILTER='$(FILTER)'\
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
TEST_EXTENSIONS = .sh .sp .gc
SP_LOG_COMPILER = $(top_srcdir)/tests/sp_log_compiler
GC_LOG_COMPILER = $(top_srcdir)/tests/gc_log_compiler


BSIM6_TESTS = \
    include.2.sp \
    bsim6_inv_dc.gc \
    bsim6_ringosc_17.gc \
    bsim6_inv_tr.gc \
    bsim6_inverter_transient.gc

TESTS = \
	$(CUSTOM_TESTS) \
	$(SH_TESTS) \
	$(GC_TESTS) \
	$(BSIM6_TESTS) \
	$(INLINE_TESTS) \
	$(GCUF_TESTS)

INLINE_TESTS = abs.gc instanceparam.gc

# tests that require an .so
SO_TESTS = \
    res.gc \
    op.gc \
    sensor.gc \
    hyst1.gc \
    lamp.gc \
    mid.gc \
    vbg.gc \
    subcoil.gc \
    subcoil.1.gc \
    rc.gc \
    short.gc \
    collapse.gc \
    rc_emb.gc \
    motor.gc
CUSTOM_TESTS = \
    if.gc \
    include.1.gc \
    include.3.gc \
    callfunctions.gc \
    coil_auto.gc \
    coil_br.gc \
    diff.2.gc \
    paramset.gc \
    sensor_chip.gc \
    vsstamp.gc \
    csstamp.gc \
    vcvsstamp.gc \
    vccsstamp.gc \
    ccvsstamp.gc \
    cccsstamp.gc \
    inline.gc \
    rc1.gc \
    tr.gc \
    tr1.gc \
    int_cont.gc \
    motor1.gc \
    adms_res.gc \
    opamp.gc
GC_TESTS = \
    $(SO_TESTS) \
    if2.gc \
    implicit.gc \
    coil.gc \
    chainvs.gc \
    ddx.gc \
    vs.gc \
    forloop.gc \
    case.gc \
    vs1.gc \
    vs2.gc \
    cap.gc \
    cap_rescale.gc \
    cs.gc \
    af.0.gc \
    af.1.gc \
    af.2.gc \
    afbuf.gc \
    afsquare.gc \
    square.gc \
    int.gc \
    ldo.gc \
    pid.gc \
    diff.gc \
    limit.gc \
    log.gc \
    meas.gc \
    mul.gc \
    cccs.gc \
    ccvs.gc \
    res2.gc \
    sin.gc \
    vccs.gc \
    vcvs.gc
# will only work with -uf
if GNUCAP_UF
GCUF_TESTS = \
    bvs.gc \
    ager.gc \
    ager.sub.gc \
    aging_mos.gc \
    rc_uic.gc \
    ddc_rc.gc \
    ddc.2.gc \
    mos.1.gc \
    mos_tt.gc \
    ttint.gc
endif

XFAIL_TESTS =

# some hack in -uf makes this work...
if !GNUCAP_UF
XFAIL_DIFF = diff.gc
endif
XFAIL_TESTS += $(XFAIL_DIFF)

# .verilog does not work
XFAIL_TESTS += include.2.sp

# not implemented yet
XFAIL_TESTS += collapse.gc

# not implemented yet
XFAIL_TESTS += diff.2.gc

# does not work.
XFAIL_TESTS += coil_br.gc

# time restore mismatch
XFAIL_TESTS += tr1.gc

POSTDIFF=
# remove cruft from model compiler
@nothing@$(INLINE_TESTS:%=%.out): POSTDIFF='1,/;·testdata/d'
@nothing@ttint.gc.out: POSTDIFF='0,/^testdata/d'
@nothing@if.gc.out: POSTDIFF='s/··[a-z][a-z][a-z][a-z]·$$/·myif/' \
                            's/=·0\.3333.*$$/=·0\.3333/'
@nothing@if2.gc.out: POSTDIFF='s/=·0\.3333.*$$/=·0\.3333/'
@nothing@inline.gc.out: POSTDIFF='0,/outputstart/d'
@nothing@aging_mos.gc.out: POSTDIFF='0,/^done·load_va/d'
@nothing@ddx.gc.out: POSTDIFF='s/-9\....E-15/·0.·······/g' \
	                       's/-27\...E-30/·0.·······/' \
	                       's/9\....E-15/0.·······/g' \
                          's/385.1E-24/0.·······/'

# ignore strobes about RSSMOD
@nothing@$(BSIM6_TESTS:%=%.out): \
	POSTDIFF='/\[BSIM6\]·Although/d' '/^gnucap-adms·/d' \
	         's/·pmosHACK·$$/pmosHACK/' \
	         's/·nmosHACK·$$/nmosHACK/' \
	         '/[a-z]·$$/s/·$$//' \
	         '/·0·$$/s/·$$//' \
	         '/BD·$$/s/·$$//' \
	         '/^>>>/s/·$$//' \
	         '/^\./s/·$$//'

@nothing@include.1.gc.out: POSTDIFF='/^gnucap-adms·/d'
@nothing@include.2.sp.out: POSTDIFF='/^gnucap-adms·/d'
@nothing@include.3.gc.out: POSTDIFF='/^gnucap-adms·/d'
@nothing@mos.1.gc.out: POSTDIFF='/^gnucap-adms·/d'

@nothing@bsim6_ringosc_17.gc.out \
         bsim6_inv_tr.gc.out \
			bsim6_inverter_transient.gc.out: BSIM6.1.1.so

@nothing@mos.1.gc.out: $(PWD)/mos.1.va

# FIXME: more automatically...?
EXTRA_DIST = \
    aging_mos.va \
    mos.1.va \
    bug1.va \
    modelcard.nmos \
    modelcard.pmos \
    $(BSIM6_TESTS) \
    $(BSIM6_TESTS:%.gc=%.ref) \
    $(GC_TESTS) \
    $(GC_TESTS:%.gc=%.ref) \
    $(GCUF_TESTS) \
    $(GCUF_TESTS:%.gc=%.ref) \
    $(INLINE_TESTS) \
    $(INLINE_TESTS:%.gc=%.ref) \
    rc_uic.gc ddc_rc.gc \
    lamp.va ager.va rc.va motor.va \
    $(SENSOR_SO:%.so=%.va) \
    $(SENSOR_SO:%.so=%vm.core.gc) \
    sensor_chip.ckt \
    callfunctions.va \
    $(DIST_TESTS) $(DIST_TESTS:%.gc=%.ref) \
    $(CUSTOM_TESTS) $(CUSTOM_TESTS:%.gc=%.ref) \
    opvm.core.gc \
    vbgvm.core.gc \
    ldovm.core.gc \
    sensorvm.core.gc \
    ldo.gc \
    ldo.ref \
    motor.gc \
    motor.ref \
    net.va \
    sp_log_compiler

vbg.out sensor_chip.log sensor_chip.out: \
    $(SENSOR_SO) $(PWD)/sensor_chip.ckt $(PWD)/opvm.core.gc $(PWD)/vbgvm.core.gc $(PWD)/ldovm.core.gc \
    vbg.so ldo.so sensor.so $(PWD)/sensorvm.core.gc

%.out %.ret 2%.out: %.gc
	$(GNUCAP) < $(lastword $+) 2> 2$*.out | tee $*.out | grep too\ negative ;\
	[ 0 -ne $$? ]

callfunctions.out: callfunctions.so

SENSOR_SO = op.so ldo.so vbg.so sensor.so

res2.log: res.so
rc1.log: rc.so
hyst1.log: mul.so

LT = <

callfunctions.gc.out: LT=2> >(tee /dev/fd/2 | grep negative | while read line; do echo "[ERROR] $$line"; done ) <

# SHS = \
#    $(GCUF_TESTS)

# obsolete, hopefully
$(SHS): %.sh: %.gc Makefile %.ref
	@echo "#!$(SHELL)" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-$(GNUCAP)}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP -a c_debug.so -a ../src/.libs/d_adms_i.so -a ../src/.libs/d_adms_vbr.so \$$GNUCAP_ARGS $(LT) $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -w -rup $(lastword $+) - || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

if GNUCAP_UF
ager.sub.log ager.log aging_mos.log mos.1.log: %.log: %.so
ager.sub.gc.out ager.gc.out aging_mos.gc.out mos.1.gc.out: %.gc.out: %.so
endif

subcoil.log subcoil.1.log vccs.log ccvs.log cap.log: %.log: %.so

include.3.gc.out: $(PWD)/includefile.v

# SO_TESTS?
# (don't remember why this is/was needed)
$(SO_TESTS:%.gc=%.log): %.log: %.so
$(SO_TESTS:%.gc=%.gc.out): %.gc.out: %.so
coil.log: coil.so
coil_auto.log: coil.so
coil_br.log: coil.so
vs.log: vs.so
cs.log: cs.so
if.log: if.so
if2.log: if2.so
motor1.log: motor.so
callfunctions.log: callfunctions.so

$(GC_TESTS:%.gc=%.log): %.log: %.so
$(GC_TESTS:%.gc=%.out): %.out: %.so

PREDIFF= \
    '1,/^core-lib·version/d' \
    '/^stashing·as/d' \
    '/>spice/s/·$$//' \
    '/>verilog/s/·$$//' \
    '/^>>>>>$$/d' \
    '/^@/d' \
    '/^default·plugins:/d' \
    '/:·already.installed,·replacing/d'

# huh?
# $(TESTS:%.sh=%.log): Makefile

check_custom: $(CUSTOM_TESTS:%.sh=%.log)
check_bsim: $(BSIM6_TESTS:%.sp=%.log)

$(BSIM6_TESTS:%=%.out): $(PWD)/modelcard.nmos $(PWD)/modelcard.pmos
# $(BSIM6_TESTS:%.sp=%.log): BSIM6.1.1.so

$(SH_TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(AM_TESTS_ENVIRONMENT)$(GNUCAP) \
		$(GNUCAP_PLUGINS:%=-a ../src/.libs/%) $(GNUCAP_ARGS) $(LT) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

%.sp.out: %.sp
	$(AM_TESTS_ENVIRONMENT)$(GNUCAP) -a ../src/.libs/d_adms_i.so -a ../src/.libs/d_adms_vbr.so -b $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $@

PLUGIN_ARGS = -a c_debug.so -a d_adms_i.so -a d_adms_vbr.so

%.gc.out: %.gc
	$(AM_TESTS_ENVIRONMENT)$(GNUCAP) $(PLUGIN_ARGS) $(LT) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $@

# ugh, also creates some hxx files!
%.cc %.h: %.va $(ADMS_GC_XML) $(top_srcdir)/src/implicit.xml
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

abs.log abs.out: abs.gc cleantmp abs.gc
instanceparam.log instanceparam.out: cleantmp
inline.log inline.out: cleantmp

cleantmp:
	rm -rf .admstmp

if VPATH_BUILD
VPATH_CLEANFILES = modelcard.*
endif

CLEANFILES = \
    $(CUSTOM_TESTS) \
	 $(GCUF_TESTS) \
    $(VPATH_CLEANFILES) \
    .*.adms *.out \
    $(GC_TESTS) \
    $(INLINE_TESTS) \
    *.hxx .*.xml \
    core \
    rc_uic.sh \
    ddc_rc.sh \
    ddc.2.sh \
    *.vams

CLEANFILES+= $(check_LTLIBRARIES:%.la=%.h) \
             $(check_LTLIBRARIES:%.la=%.cc)
CLEANFILES+= .*.adms *.out *.hxx .*.xml \
    callfunctions.sh \
    inline.sh opamp.sh sensor_chip.sh

$(CUSTOM_TESTS:%.gc=%.gc.out): .P
$(CUSTOM_TESTS:%.gc=%.log): .P
$(INLINE_TESTS:%.gc=%.gc.out): .P
$(INLINE_TESTS:%.gc=%.log): .P
$(GC_TESTS:%.gc=%.gc.out): .P
$(GC_TESTS:%.gc=%.log): .P
$(GCUF_TESTS:%.gc=%.gc.out): .P
$(GCUF_TESTS:%.gc=%.log): .P
$(SH_TESTS:%.sh=%.log): .P
$(SH_TESTS:%.sh=%.out): .P
$(BSIM6_TESTS:%=%.out): .P
# $(BSIM6_TESTS:%.sp=%.log): .P
$(BSIM6_TESTS:%.gc=%.log): .P

.P:
	@:

.PHONY: %.reg clean-examples .P

include ../adms.mk

@am__include@ @top_srcdir@/gnuplot.mk
