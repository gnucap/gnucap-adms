# this file is part of gnucap-adms
# (c) 2011-2013 Felix Salfelder
# license: GPLv3+

GNUCAP_ADMS = ../src/gnucap-adms
GNUCAP = LD_LIBRARY_PATH=../src/.libs $(shell echo gnucap | sed -e '$(transform)')

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:
check_LTLIBRARIES = if.la if2.la # bug1.la
AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libe_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

%.so: %.la
	[ -f $@ ] || $(LN_S) .libs/$@ .

if VERILOG_EXTERN
EXT_TESTS =
endif

TEST_EXTENSIONS = .sh
AM_TESTS_FD_REDIRECT = 9>&2
GC_TESTS = opamp.sh sensor_chip.sh
TESTS = $(EXT_TESTS) $(GC_TESTS)
GC_TESTS = if.sh if2.sh

EXTRA_DIST = $(TESTS:%.sh=%.ref) $(GC_TESTS:%.sh=%.gc)

$(GC_TESTS) $(EXT_TESTS): %.sh: %.ref
	@echo "#!/bin/sh" > $@
	@echo "set -e" >> $@
	@echo "exec 2>&9" >> $@
	@echo "diff -rup $< $*.out || exit 1" >> $@
	@echo "rm $*.out" >> $@
	@chmod +x $@

$(TESTS): %.sh: %.out

if.out: if.so
if2.out: if2.so

PREDIFF= \
    '1,/^\#·See·the·file/d'
    '/^stashing·as/d' \
    '/^@/d' \
    '/:·already.installed,·replacing/d'

$(TESTS:%.sh=%.out): %.out: %.gc
	$(AM_V_at)$(GNUCAP) $(GNUCAP_PLUGINS:%=-a ../plugins/.libs/%) $(GNUCAP_ARGS) $< | \
	    sed $(subst ·, ,$(PREDIFF:%=-e %)) \
	    $(if $(POSTDIFF),| sed) $(subst ·, ,$(POSTDIFF:%=-e %)) > $*.out

i=0 1 2
ADMS_ES=$(i:%=-e ../src/gnucap_%.xml)
ADMS_GC_XML=../src/gnucap_0.xml

# ugh, also creates some hxx files!
%.cc %.h: %.va $(ADMS_GC_XML)
	admsXml @ADMSXML_FLAGS@ -I$(top_srcdir)/src $< $(ADMS_ES) -o $(notdir $*)

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

CLEANFILES = .*.adms *.out $(GC_TESTS) *.hxx .*.xml
CLEANFILES+= $(check_LTLIBRARIES:%.la=%.h) \
             $(check_LTLIBRARIES:%.la=%.cc)

.PHONY: %.reg clean-examples
