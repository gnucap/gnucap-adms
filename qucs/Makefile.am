include ../adms.mk

# BUG: unhardcode, commit installDATADIR patch to qucs-core
# or pick from Make2.
QUCS_PREFIX = ${HOME}/.local

# fix the serious ones first.
AM_CXXFLAGS = -Wno-unused-variable -Wno-unused-but-set-variable

QUCS_VERILOG = \
	andor4x2.la \
	andor4x3.la \
	andor4x4.la \
	binarytogrey4bit.la \
	comp_1bit.la \
	comp_2bit.la \
	comp_4bit.la \
	dff_SR.la \
	DLS_1ton.la \
	DLS_nto1.la \
	dmux2to4.la \
	dmux3to8.la \
	dmux4to16.la \
	fa1b.la \
	fa2b.la \
	gatedDlatch.la \
	greytobinary4bit.la \
	ha1b.la \
	hicumL0V1p12.la \
	hicumL0V1p2g.la \
	hicumL2V2p22.la \
	hicumL0V1p2.la \
	hicumL0V1p3.la \
	hicumL2V2p11.la \
	hpribin4bit.la \
	jkff_SR.la \
	log_amp.la \
	logic_0.la \
	logic_1.la \
	MESFET.la \
	mod_amp.la \
	mux2to1.la \
	mux4to1.la \
	mux8to1.la \
	nigbt.la \
	pad2bit.la \
	pad3bit.la \
	pad4bit.la \
	photodiode.la \
	phototransistor.la \
	potentiometer.la \
	tff_SR.la \
	EKV26MOS.la \
	vcresistor.la

QUCS_NOTYET=

#guesstopology issue
QUCS_NOTYET+= \
	bsim3v34pMOS.la \
	bsim4v30nMOS.la

# has an external temperature node...
QUCS_NOTYET+= fbh_hbt-2_2a.la

# scoping problem. model variable assignment in eval (Tdev)
QUCS_NOTYET+= hicumL2V2p31n.la \
	hicumL2V2p24.la \
	hicumL2V2p23.la

#adms segfaults
QUCS_NOTYET+= \
	bsim4v30pMOS.la \
	bsim3v34nMOS.la


try: $(QUCS_NOTYET)

QUCS_SHAREDDIR = $(QUCS_PREFIX)/share
QUCS_DATADIR = $(QUCS_SHAREDDIR)/qucs-core/verilog/

VPATH+= $(QUCS_DATADIR)

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:

AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

noinst_LTLIBRARIES = $(QUCS_VERILOG)
EXTRA_LTLIBRARIES = $(QUCS_NOTYET)

CLEANFILES=*.hxx *.vams *.cc *.h
