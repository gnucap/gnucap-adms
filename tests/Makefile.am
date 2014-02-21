# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
# GNUCAP = LD_LIBRARY_PATH=../src/.libs
GNUCAP = $(shell echo gnucap | sed -e '$(transform)')

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = \
    if.la if2.la \
    vs.la cs.la cap.la \
    coil.la ccvs.la vccs.la \
    analogfunction.la \
    int.la pid.la vcvs.la \
    diff.la \
    mul.la \
    subcoil.la

AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

if VERILOG_EXTERN
EXT_TESTS =
endif

AM_TESTS_ENVIRONMENT = \
    export PATH='../src:$(PATH)'\
           LD_LIBRARY_PATH='../src/.libs'\
           GNUCAP='$(GNUCAP)'\
           REDIRECT='exec 2>&9'\
           GNUCAP_CPPFLAGS='-I$(abs_top_builddir)/src -I$(abs_top_srcdir)/src'\
           srcdir='$(srcdir)';
AM_TESTS_FD_REDIRECT = 9>&2
TEST_EXTENSIONS = .sh
TESTS = $(EXT_TESTS) $(GC_TESTS) $(GCUF_TESTS)
GC_TESTS = if.sh if2.sh coil.sh vs.sh cs.sh cap.sh \
    analogfunction.sh \
    int.sh \
    pid.sh \
	 diff.sh \
	 mul.sh \
    ccvs.sh \
    vccs.sh \
    vcvs.sh \
    subcoil.sh

# will only work with -uf
GCUF_TESTS = rc_uic.sh ddc_rc.sh ddc.2.sh

EXTRA_DIST = $(TESTS:%.sh=%.ref) $(GC_TESTS:%.sh=%.gc) rc_uic.gc ddc_rc.gc
FILTER = cat

$(GCUF_TESTS) $(GC_TESTS) $(EXT_TESTS): %.sh: %.gc Makefile %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "GNUCAP=\$${GNUCAP-$(GNUCAP)}" >> $@
	@echo "eval \$$REDIRECT" >> $@
	@echo "\$$GNUCAP \$$GNUCAP_ARGS $(LT) $< | sed $(subst ·, ,$(PREDIFF:%=-e %)) \\" >> $@
	@echo "$(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out" >> $@
	@echo "$(FILTER) $*.out < $(lastword $+) | diff -rup $(lastword $+) - || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

vccs.log ccvs.log cap.log: %.log: %.so
coil.log: coil.so
vs.log: vs.so
cs.log: cs.so
if.log: if.so
if2.log: if2.so

$(GC_TESTS:%.sh=%.log): %.log: %.so
$(GC_TESTS:%.sh=%.out): %.out: %.so

PREDIFF= \
    '1,/^\#·See·the·file/d'
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

$(TESTS:%.sh=%.log): Makefile

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../plugins/.libs/%) $(GNUCAP_ARGS) $< | \
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

%.reg: %.test %.ref
	diff $+

%.test: %.gc
	$(GNUCAP) $< | grep '^[a-zA-Z0-9]*=' > $@

check: $(CHECK:%=%.reg)

clean-local: clean-examples

clean-examples:
	rm -f *.out *.so *.cc *.hxx *.h

CLEANFILES = .*.adms *.out $(GC_TESTS) *.hxx .*.xml \
    core \
    rc_uic.sh \
    ddc_rc.sh
CLEANFILES+= $(check_LTLIBRARIES:%.la=%.h) \
             $(check_LTLIBRARIES:%.la=%.cc)

$(TESTS:%.sh=%.log): .P
$(TESTS:%.sh=%.out): .P

.P:
	@:

.PHONY: %.reg clean-examples .P

@am__include@ @top_srcdir@/gnuplot.mk
