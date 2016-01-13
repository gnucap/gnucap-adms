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
	vcresistor.la \
	bsim3v34nMOS.la \
	bsim3v34pMOS.la

QUCS_NOTYET=

# has an external temperature node...
QUCS_NOTYET+= fbh_hbt-2_2a.la

# scoping problem. model variable assignment in eval (Tdev)
QUCS_NOTYET+= hicumL2V2p31n.la \
	hicumL2V2p24.la \
	hicumL2V2p23.la

#adms segfaults
QUCS_NOTYET+= \
	bsim4v30pMOS.la \
	bsim4v30nMOS.la


try: $(QUCS_NOTYET)

QUCS_SHAREDDIR = $(QUCS_PREFIX)/share
QUCS_DATADIR = $(QUCS_SHAREDDIR)/qucs-core/verilog/

VPATH+= $(QUCS_DATADIR)

AM_DEFAULT_SOURCE_EXT = .va
.va.cc:

# stated in adms.mk. does not work?
%.cc %.h: %.va $(top_srcdir)/src/implicit.xml

AM_CPPFLAGS = @GC_CPPFLAGS@ -I$(top_srcdir)/src
LIBS = ../src/libgnucap_adms.la
AM_LDFLAGS = -module -avoid-version -rpath /dev/null

noinst_LTLIBRARIES = $(QUCS_VERILOG)

nodist_andor4x2_la_SOURCES = andor4x2.va
nodist_andor4x3_la_SOURCES = andor4x3.va
nodist_andor4x4_la_SOURCES = andor4x4.va
nodist_binarytogrey4bit_la_SOURCES = binarytogrey4bit.va
nodist_comp_1bit_la_SOURCES = comp_1bit.va
nodist_comp_2bit_la_SOURCES = comp_2bit.va
nodist_comp_4bit_la_SOURCES = comp_4bit.va
nodist_dff_SR_la_SOURCES = dff_SR.va
nodist_DLS_1ton_la_SOURCES = DLS_1ton.va
nodist_DLS_nto1_la_SOURCES = DLS_nto1.va
nodist_dmux2to4_la_SOURCES = dmux2to4.va
nodist_dmux3to8_la_SOURCES = dmux3to8.va
nodist_dmux4to16_la_SOURCES = dmux4to16.va
nodist_fa1b_la_SOURCES = fa1b.va
nodist_fa2b_la_SOURCES = fa2b.va
nodist_gatedDlatch_la_SOURCES = gatedDlatch.va
nodist_greytobinary4bit_la_SOURCES = greytobinary4bit.va
nodist_ha1b_la_SOURCES = ha1b.va
nodist_hicumL0V1p12_la_SOURCES = hicumL0V1p12.va
nodist_hicumL0V1p2g_la_SOURCES = hicumL0V1p2g.va
nodist_hicumL2V2p22_la_SOURCES = hicumL2V2p22.va
nodist_hicumL0V1p2_la_SOURCES = hicumL0V1p2.va
nodist_hicumL0V1p3_la_SOURCES = hicumL0V1p3.va
nodist_hicumL2V2p11_la_SOURCES = hicumL2V2p11.va
nodist_hpribin4bit_la_SOURCES = hpribin4bit.va
nodist_jkff_SR_la_SOURCES = jkff_SR.va
nodist_log_amp_la_SOURCES = log_amp.va
nodist_logic_0_la_SOURCES = logic_0.va
nodist_logic_1_la_SOURCES = logic_1.va
nodist_MESFET_la_SOURCES = MESFET.va
nodist_mod_amp_la_SOURCES = mod_amp.va
nodist_mux2to1_la_SOURCES = mux2to1.va
nodist_mux4to1_la_SOURCES = mux4to1.va
nodist_mux8to1_la_SOURCES = mux8to1.va
nodist_nigbt_la_SOURCES = nigbt.va
nodist_pad2bit_la_SOURCES = pad2bit.va
nodist_pad3bit_la_SOURCES = pad3bit.va
nodist_pad4bit_la_SOURCES = pad4bit.va
nodist_photodiode_la_SOURCES = photodiode.va
nodist_phototransistor_la_SOURCES = phototransistor.va
nodist_potentiometer_la_SOURCES = potentiometer.va
nodist_tff_SR_la_SOURCES = tff_SR.va
nodist_EKV26MOS_la_SOURCES = EKV26MOS.va
nodist_vcresistor_la_SOURCES = vcresistor.va
nodist_bsim3v34pMOS_la_SOURCES = bsim3v34pMOS.va
nodist_bsim3v34nMOS_la_SOURCES = bsim3v34nMOS.va

# EXTRA_LTLIBRARIES = $(QUCS_NOTYET)

CLEANFILES=*.hxx *.vams *.cc *.h
