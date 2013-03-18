# Makefile for gnuplot related targets
# (c) 2012 felix salfelder
#
STYLE=data lines
GNUPLOT=gnuplot
SHELL=/bin/bash


# get 4-5.1-3,5,8-9.3
# return {4..5}:{{1..3},5,{8..9}}:{3..3}
define rowlist_sh
	comma=""; \
	l=""; rbr=""; lbr=""; \
	IFS=,; b=$$rowi; a=($$b); unset IFS; \
	for i in $${a[*]}; do \
		IFS=-; j=($$i); unset IFS;\
		la=$${i[0]}; \
		last=`expr $${#j[*]} - 1`; \
		S=$$( seq -s, $${j[0]} $${j[$$last]} );\
		S={$${j[0]}..$${j[$$last]}};\
		l=$$l$$comma$$S; \
		rowlist_=$$lbr$$l$$rbr; \
		comma=","; lbr="{"; rbr="}"; \
	done; 
endef

define do_pscols2
$(shell echo processing $(1)>&2; \
	colon=;\
	pscols_tmp=$$( echo $(1) | sed -r 's/[^.]*((\.[0-9,-]+)+).*/\1/'); \
	echo cols_raw $$pscols_tmp>&2; \
	IFS=.; C=($$pscols_tmp); unset IFS;\
	echo cols_raw $${C[*]}>&2; \
	ret=;\
	for rowi in $${C[*]}; do \
		rowlist_=; $(rowlist_sh) \
		ret=$$ret$$colon$$rowlist_; \
		colon=:; \
	done; echo $$ret; )
endef

sourcepattern='s/\(\(^\.\)*\)\.[0-9]\+\(\.[0-9,-]\+\)\+/\1.out/'

#FIXME: eps
.SECONDEXPANSION:
%.ps: $$(shell echo % | sed -e $(sourcepattern) )
	CH="";\
	head -n1 $< | grep '^[0-9\.Ee \t]*$$' >/dev/null; \
	[[ $$? = 0 ]] || CH="set key autotitle columnhead;"; \
	[[ $* =~ [^.]*\.[0-9]\.[0-9]\.[0-9] ]] && S=s; \
	IFS=,; spl=$*; a=($$spl); unset IFS; \
	# $*; \
	pscols_=$$(eval echo "\'$<\'\ using\ $(call do_pscols2,$*)\,"); \
	pscols_=$${pscols_%,} ; \
	echo -e "set style ${STYLE}; \n \
	$$CH \n \
	set datafile commentschars \"\" \n \
  	set terminal postscript ${ORIENT} enhanced ${COLOR} dashed lw ${LW} 'Helvetica' 14; \n \
  	${ECMD}  set output '$@'; \n \
 	$${S}plot $$pscols_" | $(GNUPLOT)

