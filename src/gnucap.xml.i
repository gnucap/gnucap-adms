<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- *** THIS FILE IS LICENSED UNDER THE TERMS OF THE GNU PUBLIC LICENSE v3 OR ABOVE *** -->

<!DOCTYPE admst SYSTEM "admst.dtd">
<admst version="2.3.0" xmlns:admst="http://mot-adms.sourceforge.net/xml-files/admst">
<!-- version: not renamed, load: renamed fname -->
<admst:message test="[/dbg_xml='yes']" format="**enter scope1**\n"/>

<!-- *********************************************************************************** -->
<!-- ******************************** global settings ********************************** -->
<!-- *********************************************************************************** -->

<!-- output filename is 1st argument in -o flag in command line-->
<admst:variable name="filename" select="%(/module/name)_common"/>
<admst:variable name="filename_state" select="0"/>
<admst:for-each select="/argv">
  <admst:variable name="arg" select="%(.)"/>

  <admst:if test="[position(.)=2]">
    <admst:variable name="filename" select="$arg"/>
  </admst:if>
  <admst:if test="[$(filename_state)='1']">
    <admst:variable name="filename" select="$arg"/>
    <admst:variable name="filename_state" select="2"/>
  </admst:if>
  <admst:if test="[$(filename_state)='0' and $(arg)='-o']">
    <admst:variable name="filename_state" select="1"/>
  </admst:if>
</admst:for-each>

<!-- prefix for local variables -->
<admst:variable name="var_prefix" select="__"/>
<admst:variable name="instancevar" select="c->"/>
<admst:variable name="modelvar" select="m->"/>
<admst:variable name="node_prefix" select="n_"/>
<admst:variable name="cont_prefix" select="c_"/>
<admst:variable name="state_prefix" select="ddt_"/>
<admst:variable name="ddt_type" select="x_"/>

<!-- variables used by tree traverse (by guesstopology?) -->
<admst:variable name="localvariables"/>
<admst:variable name="ddt_cur" select="0"/>
<admst:variable name="ddts"/>

<!-- variables used by tree traverse for hxx -->
<admst:variable name="_DYNAMIC" select="yes"/>
<admst:variable name="_STATIC" select="yes"/>
<admst:variable name="_DERIVATE" select="yes"/>
<admst:variable name="_DERIVATE2" select="yes"/>
<admst:variable name="_DERIVATEFORDDX" select="yes"/>
<admst:variable name="_DDX" select="yes"/>
<admst:variable name="eval_mode"/> <!-- ac, tr, dc?, ... -->
<admst:variable name="tr_begin"/>

<admst:variable name="requiredderivateforddx" select="no"/>
<admst:variable name="SkipFVariable" select="n"/>
<admst:variable name="ddxinsidederivate" select="no"/>
<admst:variable name="SkipFVariable" select="n"/>
<admst:variable name="pprobe"/>
<admst:variable name="qprobe"/>
<admst:variable name="e"/>
<admst:variable name="ep"/>
<admst:variable name="eq"/>
<admst:variable name="epq"/>
<admst:variable name="keep_ic" select="no"/>
<admst:variable name="ddtp" select="1."/>
<admst:variable name="pddt" select="-1"/>
<admst:variable name="eval_kept" select="no"/>
<admst:variable name="globalanalogfunction"/>

<!-- result for tree traverse (ngspiceVersion vs .hxx) -->
<admst:variable name="str_res_version" select=""/>
<admst:variable name="str_res_hxx_ac" select=""/>
<admst:variable name="str_res_hxx_tr" select=""/>

<!-- noise -->
<admst:value-to select="/simulator/package_name" value="fucap"/>
<admst:value-to select="/simulator/package_tarname" value="fucap"/>
<admst:value-to select="/simulator/package_version" value="0.0.1"/>
<admst:value-to select="/simulator/package_string" value="fucap 0.0.1"/>
<admst:value-to select="/simulator/package_bugreport" value="felix@salfelder.org"/>
<admst:for-each select="/module[extern='no']">
  <admst:variable name="has_sckt" select="%(count(instance))"/>

  <admst:new datatype="list" arguments="fnoise"><admst:variable name="fnoise" select="%(.)"/></admst:new>
  <admst:new datatype="list" arguments="tnoise"><admst:variable name="tnoise" select="%(.)"/></admst:new>
  <admst:new datatype="list" arguments="wnoise"><admst:variable name="wnoise" select="%(.)"/></admst:new>
  <admst:for-each select="contribution">
    <admst:variable name="contribution" select="%(.)"/>
    <admst:variable name="dependency" select="%(math/dependency)"/>
    <admst:choose>
      <admst:when test="rhs/tree/adms[datatypename='function']/..[name='flicker_noise']">
        <admst:push into="$fnoise/item" select="$contribution" onduplicate="ignore"/>
      </admst:when>
      <admst:when test="[$dependency='constant']/rhs/tree/adms[datatypename='function']/..[name='white_noise']">
        <admst:push into="$tnoise/item" select="$contribution" onduplicate="ignore"/>
      </admst:when>
      <admst:when test="[$dependency!='constant']/rhs/tree/adms[datatypename='function']/..[name='white_noise']">
        <admst:push into="$wnoise/item" select="$contribution" onduplicate="ignore"/>
      </admst:when>
    </admst:choose>
  </admst:for-each>
</admst:for-each>

<!-- create attribute with module name -->
<admst:for-each select="/module[extern='no']">
  <admst:new datatype="attribute" arguments="ngspicename">
  <admst:push into="../attribute" select="." onduplicate="abort"/>
  <admst:value-of select="../name"/>
  <admst:value-to select="value" value="%s"/>
  </admst:new>
</admst:for-each>

<!-- list of pointers for all ic (initial conditions) -->
<!-- bugs: a) already has empty or may be null element, so count(list)=1 -->
<!--       b) admst:for-each on this list works, but gives many warnings -->
<admst:new datatype="list" arguments="function">
	<admst:variable name="ddt_ic_list" select="" />
</admst:new>

<!--          thats why this list used only to filter contributions list -->

<admst:new datatype="list" arguments="contribution">
  <admst:variable name="adms_ic_list" select="" />
</admst:new>

<!-- return x = "yes" for contribution that is in uic-list -->
<admst:template match="is_ic_contribution">
  <admst:variable name="local_ret" value="no"/>
  <admst:variable name="local_id" value="%(id(.))"/>
  <admst:for-each select="$(adms_ic_list)[id(.)=$local_id]">
    <admst:variable name="local_ret" value="yes"/>
  </admst:for-each>
  <admst:return test="[$local_ret='yes']" name="x" string="yes"/>
  <admst:return test="[$local_ret='no']" name="x" string="no"/>
</admst:template>

<!-- to select all contributions in list use: -->
<!-- <admst:for-each select="contribution[is_ic_contribution(.)/[name='x']/value='yes']"> -->

<!-- ic probe names for ic ddt-->
<admst:variable name="ddt_p"/>
<admst:variable name="ddt_n"/>
<admst:template match="ic_probe_names">
  <admst:if test="[../datatypename='contribution']">
    <admst:variable name="ddt_p" select="%(../lhs/branch/pnode)" />
    <admst:variable name="ddt_n" select="%(../lhs/branch/nnode)" />
  </admst:if>
  <admst:if test="[./arguments[position(.)=1]/datatypename='probe']">
    <admst:variable name="ddt_p" select="%(./arguments[position(.)=1]/branch/pnode)" />
    <admst:variable name="ddt_n" select="%(./arguments[position(.)=1]/branch/nnode)" />
  </admst:if>
</admst:template>

<!-- *********************************************************************************** -->
<!-- *********************************** my utils ************************************** -->
<!-- *********************************************************************************** -->

<admst:variable  name="myF_par1" select=""/>
<admst:variable  name="myF_par2" select=""/>
<admst:variable  name="myF_ret1" select=""/>


<admst:template match="utils_compare_nodes_position">
<admst:message test="[/dbg_xml='yes']" format="*utils_compare_nodes_position*\n"/>
	<admst:variable  name="myF_ret1" select="not_found"/>

	<admst:if test="[$(myF_par1)=$(myF_par2)]">
		<admst:variable  name="myF_ret1" select="equal"/>
	</admst:if>

	<admst:for-each select="node[grounded='no']">
		<admst:variable name="local_name" select="%(name)"/>
		<admst:if test="[$myF_ret1='not_found']">
			<admst:if test="[$myF_par1=$local_name]">
				<admst:variable  name="myF_ret1" select="smaller"/>
			</admst:if>
			<admst:if test="[$myF_par2=$local_name]">
				<admst:variable  name="myF_ret1" select="bigger"/>
			</admst:if>
		</admst:if>
	</admst:for-each>
</admst:template>


<!-- *********************************************************************************** -->
<!-- ************************************** Make File .h ******************************* -->
<!-- *********************************************************************************** -->

<admst:open file="$(filename).h">
#include &lt;stdio.h&gt;
#include  &quot;math.h&quot;

#include &quot;u_xprobe.h&quot;
#include &quot;d_subckt.h&quot;
#include &quot;e_storag.h&quot;
#include &quot;e_model.h&quot;
#include &quot;e_compon.h&quot;
#include &quot;e_adms.h&quot;
#include &quot;u_parameter.h&quot;

<admst:for-each select="/module[extern='no']">
  <admst:value-of select="attribute[name='ngspicename']/value"/>
  <admst:variable name="module" select="%s"/>
  <admst:variable name="DEV_NAME" select="DEV_%(upper-case(name))"/>
  <admst:variable name="MODEL_NAME" select="MODEL_%(upper-case(name))"/>
  <admst:variable name="COMMON_NAME" select="COMMON_%(upper-case(name))"/>
#ifndef $(module)_header__
#define $(module)_header__

// customization - gnucap specific constants
#define VALUE_NAME &quot;#&quot;
#define TAIL_SIZE 1
#define _DDX
#define _DERIVATEFORDDX

#define SPICE_LETTER &quot;\\0&quot;
#define DEVICE_TYPE &quot;dev_$(module)&quot;
#define MODEL_TYPE &quot;$(module)&quot;

#define CONSTCtoK (273.15) // also in gnucap/constans.h P_CELSIUS0

#define EXIT_IF_ISNAN(var) assert(is_number(var));
inline void _used(double)  {}


<admst:for-each select="assignment">
  <admst:apply-templates select="rhs/tree" match="ddtcollect"/>
</admst:for-each>
<admst:for-each select="contribution">
  <admst:apply-templates select="rhs/tree" match="ddtcollect"/>
</admst:for-each>

// counted ddts %(count($ddts))
<admst:variable name="_num_states" select="%(count($ddts))"/>
// there are nodes that work different... hmmm
// BUG: number of stamp nodes is topology dependent!
<admst:variable name="_num_stamp_nodes" select="%(count(node[grounded='no' and discipline/name!='degradational']))"/>

class $(DEV_NAME);
class $(MODEL_NAME) : public MODEL_CARD {
	private:
#ifdef HAVE_PARA_BASE
		static std::map&lt;string, PARA_BASE $(MODEL_NAME)::*&gt; param_dict;
		static std::map&lt;string, PARA_BASE $(MODEL_NAME)::*&gt; param_dict_low;
#endif
	protected:
		explicit $(MODEL_NAME)(const $(MODEL_NAME)&amp; p);	// for clone
	public:
		explicit $(MODEL_NAME)(const $(DEV_NAME)* p);	// for dispatcher
		~$(MODEL_NAME)();
	public: // override virtual
		MODEL_CARD* clone()const {return new $(MODEL_NAME)(*this);}
	public:
		void precalc_first();

	public: // type
		void set_dev_type(const std::string&amp; nt);
		std::string dev_type()const	{ return _key;}

	public: // parameters
		int param_count() const;
		bool param_is_printable(int i) const;
		std::string param_name(int i) const;
		std::string param_name(int i, int j)const;
		std::string param_value(int i) const;
		void set_param_by_index(int, std::string&amp;, int);
		void set_param_by_name(std::string Name, std::string Value);
		bool IsParamGiven(const char *name) const;

  <admst:for-each select="variable[input='yes' and scope!='global_instance']">
    <admst:choose>
      <admst:when test="[type='real']">
        <admst:text format="\t\tPARAMETER&lt;double&gt; $(var_prefix)%(.);"/>
      </admst:when>
      <admst:when test="[type='integer']">
        <admst:text format="\t\tPARAMETER&lt;int&gt; $(var_prefix)%(.);"/>
      </admst:when>
      <admst:otherwise>
        <admst:fatal format="type of verilog-parameter is not real or interger\n"/>
      </admst:otherwise>
		</admst:choose> // %(input) %(scope) %(type) 

	</admst:for-each>

	public: // global_model

  <admst:for-each select="variable[input='no' and scope='global_model']">
    <admst:choose>
    <admst:when test="[type='real']">
      <admst:text format="\t\tdouble $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:when test="[type='integer']">
      <admst:text format="\t\tint $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:text format="\n\t\tstd::string $(var_prefix)%(.);\n"/>
    </admst:otherwise>
    </admst:choose>
  </admst:for-each>


	public:
		std::string _key;
		std::string _level;
};
<!-- fixme: use probe name macro... -->

/* ========================================================= */
<!-- create list for initial conditions, for now. only very specific patterns are supported -->
<!-- maybe can be cleared with: <admst:variable select="adms_ic_list" value="" /> -->
<admst:for-each select="contribution">
  <!-- catch I(a,b) <+ ddt(V(a,b))  -->
  <admst:if test="[./rhs/tree/datatypename='function' and ./rhs/tree/name='ddt']">
    <admst:push into="$(ddt_ic_list)" select="./rhs/tree" onduplicate="ignore"/>
    <admst:message format="ic candidate1 %(./rhs/tree) %(adms/id(.))\n"/>
    <admst:if test="[./rhs/tree/name='ddt']">
      <admst:if test="[./rhs/tree/arguments[position(.)=1]/datatypename='probe']">
        <admst:if test="[./rhs/tree/arguments[position(.)=1]/branch/pnode=./lhs/branch/pnode and ./rhs/tree/arguments[position(.)=1]/branch/nnode=./lhs/branch/nnode]">
          <admst:push into="$(adms_ic_list)" select="." onduplicate="ignore"/>
          <admst:message format="ic candidate2 %(.) %(adms/id(.))\n"/>
        </admst:if>
      </admst:if>
    </admst:if>
  </admst:if>
  <!-- catch I(a,b) <+ c*ddt(V(a,b))  -->
  <admst:if test="[./rhs/tree/datatypename='mapply_binary']">
    <admst:if test="[./rhs/tree/name='multtime']">
      <admst:if test="[./rhs/tree/arg1/datatypename='variable' and ./rhs/tree/arg2/datatypename='function' ]">
        <admst:message format="datatypename %(./rhs/tree/arg1/datatypename) : %(./rhs/tree/arg2/datatypename)  expr %(./rhs/tree/name)  \n"/>
        <admst:if test="[./rhs/tree/arg2/name='ddt']">
          <admst:if test="[./rhs/tree/arg2/arguments[position(.)=1]/datatypename='probe']">

<!-- hmm, how to store probes??
<admst:new datatype="probe" arguments="foobar">
  <admst:variable name="probeA" />
 select=./lhs/branch/pnode
</admst:new>
-->
              <admst:push into="$(ddt_ic_list)" select="./rhs/tree/arg2" onduplicate="ignore"/>

            <admst:if test="[./rhs/tree/arg2/arguments[position(.)=1]/branch/pnode=./lhs/branch/pnode
                        and ./rhs/tree/arg2/arguments[position(.)=1]/branch/nnode=./lhs/branch/nnode]">
              <admst:push into="$(adms_ic_list)" select="." onduplicate="ignore"/>
              <admst:message format="ic candidate3 %(.) %(adms/id(.))\n"/>
            </admst:if>
          </admst:if>
        </admst:if>
      </admst:if>

    </admst:if>
  </admst:if>
</admst:for-each>

class $(DEV_NAME) : public ADMS_BASE
{
	private:
		std::string _modelname;
		node_t _nodes[ %(count(node[grounded='no'])) ];

	private:
		explicit $(DEV_NAME)(const $(DEV_NAME)&amp; p);
	public:
		explicit $(DEV_NAME)();
		~$(DEV_NAME)();

	protected: // override virtual
		char id_letter()const	{untested();return SPICE_LETTER[0];}
		bool print_type_in_spice()const {return true;}
		std::string value_name()const {return VALUE_NAME;}
		uint_t max_nodes()const {return %(count(node[grounded='no' and location='external']));}
		uint_t min_nodes()const {return 0;}
		uint_t matrix_nodes()const {return %(count(node[grounded='no']));}
		// uint_t matrix_nodes()const {return net_nodes() + int_nodes();}
		uint_t net_nodes()const {return _net_nodes;}
		uint_t int_nodes()const {return %(count(node[grounded='no' and location='internal']));}
		CARD* clone()const {return new $(DEV_NAME)(*this);}

		void precalc_first();
		void expand();
		void precalc_last();
		void internal_precalc();

	public: // tr
		void tr_iwant_matrix();
		void tr_begin();

		void tr_restore();
		void store_values();
		void dc_advance();
		void tr_advance();
		void tr_regress();
		bool tr_needs_eval()const;
		bool do_tr();
		void tr_load();
		void tr_queue_eval();
		void tr_eval();
		void tr_eval_kept();
		bool conv_check()const;
		TIME_PAIR tr_review();
		void tr_accept();
		void tr_unload();
		double tr_involts()const	{unreachable();return NOT_VALID;}
		double tr_involts_limited()const {unreachable();return NOT_VALID;}
		double tr_amps()const	{itested();return NOT_VALID;}
		double tr_probe_num(const std::string&amp;)const;

	public: // ac
		void ac_iwant_matrix();
		void ac_begin();
		void do_ac();
		void ac_eval(); // -&gt; COMMON?
		void ac_load();
		COMPLEX ac_involts()const {unreachable();return NOT_VALID;}
		COMPLEX ac_amps()const {unreachable();return NOT_VALID;}
		XPROBE  ac_probe_ext(const std::string&amp;)const {itested(); return XPROBE(NOT_VALID, mtNONE);}
		uint_t tail_size()const {return TAIL_SIZE;}

	public: // tt
		void tt_advance();
		void stress_apply();
		void tt_accept();
		void tt_begin();
		void tr_stress_last(); // tt_stress??

	public: // misc
		void keep_ic();
		bool has_memory(){return $_num_states;}

	public: // ports
		std::string port_name(uint_t i)const {
			assert(i &lt; max_nodes());
			return node_name(i);
		}
	protected: //nodes
		std::string node_name(uint_t i)const {
			static std::string names[] = {
<admst:join select="node[grounded='no' and location='external']" separator=",">
				"%(name)"
</admst:join>
<admst:for-each select="node[grounded='no' and location='internal']">
   	<admst:value-of select="name"/>
   	<admst:text format=",\n\t\t\t\t&quot;%s&quot;"/>
</admst:for-each>
			};
			itested();
			assert(i &lt; max_nodes()+int_nodes());
			return names[i];
		}

	private: // parameters
		void set_param_by_name(std::string Name, std::string Value);

	protected: // implementation
		void guesstopology();
		void acLoad();
		void pzLoad(COMPLEX s);
		enum ELoadModes { // FIXME: remove
			MODE = 0x3, 				// old 'mode'
			MODETRAN = 0x1,
			MODEAC = 0x2,
			MODEDC = 0x70,				// old 'modedc'
			MODEDCOP = 0x10,
			MODETRANOP = 0x20,
			MODEDCTRANCURVE = 0x40,
			INITF = 0x3f00,				// old 'initf' parameters
			MODEINITFLOAT = 0x100,
			MODEINITJCT = 0x200,
			MODEINITFIX = 0x400,
			MODEINITSMSIG = 0x800,
			MODEINITTRAN = 0x1000,
			MODEINITPRED = 0x2000,
			MODEUIC = 0x10000l		// old 'nosolv' paramater
		};
		int load(int mode);

	private: // subdevices
<admst:for-each select="instance">
		COMPONENT* _%(name);
</admst:for-each>

	private: // global_instance (ask?)

<admst:for-each select="variable[input='no' and scope='global_instance']">
  <admst:choose>
    <admst:when test="[type='real']"><admst:text format="\n\t\tdouble $(var_prefix)%(.);\n"/></admst:when>
    <admst:when test="[type='integer']"><admst:text format="\n\t\tint $(var_prefix)%(.);\n"/></admst:when>
    <admst:otherwise><admst:text format="\n\t\tstd::string $(var_prefix)%(.);\n"/></admst:otherwise>
  </admst:choose>
</admst:for-each>

	public: // global_model... hmmm

  <admst:for-each select="variable[input='yes' and scope='global_model']">
    <admst:choose>
    <admst:when test="[type='real']">
      <admst:text format="\t\tdouble $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:when test="[type='integer']">
      <admst:text format="\t\tint $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:text format="\n\t\tstd::string $(var_prefix)%(.);\n"/>
    </admst:otherwise>
    </admst:choose>
  </admst:for-each>


	private: // analog functions from va
	<admst:apply-templates select="/module" match="analogfunctionH"/>

	private: // node numbers
		// stick to gnucap node ordering (external first)
		enum __Node {
  <admst:join select="node[grounded='no' and location='external']" separator=",">
    <admst:text format="\n\t\t\t__%(name)=%(position(.)-1)"/>
  </admst:join>
  <admst:for-each select="node[grounded='no' and location='internal']">,
			__%(name) /* %(location) %(direction) %(discipline) */
  </admst:for-each>,
			_num_nodes
		};
		// map nodenames to itl slots
		// some disciplines dont have itls
		enum __CNode {
  <admst:join select="node[grounded='no' and discipline/name!='degradational' and location='external']" separator=",">
    <admst:text format="\n\t\t\t$(cont_prefix)%(name)=%(position(.)-1)"/>
  </admst:join>
  <admst:for-each select="node[grounded='no' and discipline/name!='degradational' and location='internal']">,
			$(cont_prefix)%(name) /* %(location) %(direction) %(discipline) */
  </admst:for-each>
		};
	protected: // nodes sorted according to discipline
		enum ENode {
			n_GND=-1,
  <admst:join select="node[discipline/name='electrical' and grounded='no']" separator=",">
    <admst:text format="\n\t\t\t$node_prefix%(name)=__%(name)"/>
  </admst:join>
		};
		enum TNode {
  <admst:for-each select="node[discipline/name='thermal' and grounded='no']" separator=",">
    <admst:text format="\n\t\t\t$node_prefix%(name)=__%(name),"/>
  </admst:for-each>
		};
  <admst:if test="[count(node[discipline/name='degradational' and grounded='no'])!=0]">
		enum ANode {
    <admst:for-each select="node[discipline/name='degradational' and grounded='no']" separator=",">
        <admst:text format="\n\t\t\t$node_prefix%(name)=__%(name),"/>
      </admst:for-each>
    <admst:text format="\n\t\t\t_num_adp_nodes=%(count(node[discipline/name='degradational' and grounded='no']))"/>
		};
  </admst:if>

	private: // nodes / matrix indexes

// raw value
#define NR(p) _NR($(node_prefix) ## p)
inline double _NR(int p){ return _n[p].v0(); } 

#define BR(p,n) _BR($(node_prefix) ## p, $(node_prefix) ## n)
inline double _BR(int p, int n) {
	if (n==-1)
		return _n[p].v0();
	return volts_limited(_n[p],_n[n]);
}

// discipline wrapped value
#define NP(p) _NP($(node_prefix) ## p)
inline voltage_t _NP(ENode p){
	return std::min(1000., std::max(-1000., _n[p].v0()) );
}
inline double _NP(TNode p){
	assert(is_number(_n[p].v0()));
	return _n[p].v0() + P_CELSIUS0;
}
<admst:if test="[count(node[discipline/name='degradational' and grounded='no'])!=0]">
inline double _NP(ANode p){
	trace1("NP(ANode)", _n[p].m_());
	assert(is_number(_n[p].tt()));
	return _n[p].tt();
}
</admst:if>
#define BP(p,n) _BP($(node_prefix) ## p, $(node_prefix) ## n)
inline voltage_t _BP(ENode p, ENode n) { return volts_limited(_n[p],_n[n]);}

#define _N(n) _n[$(node_prefix) ## n]

	private:
		double it0[$_num_stamp_nodes];
		double it1[$_num_stamp_nodes];
		double i_kept[$_num_stamp_nodes];

		enum EMatrixEntries {

    <admst:variable name="local_var" select="II"/>
    <admst:for-each select="jacobian">
      <admst:value-of select="column/name"/>
      <admst:value-of select="row/name"/>
      <admst:text format="m_%s_%s=%(position(.)-1)"/>

      <admst:variable name="local_var" select="$(local_var)I"/>
      <admst:if test="[position(.)!=count(.)]">
        <admst:text format=", "/>
        <admst:if test="[$(local_var)='IIIIIII']">
          <admst:variable name="local_var" select="I"/>
          <admst:text format="\n\t\t\t"/>
        </admst:if>
      </admst:if>
    </admst:for-each>
		};

	bool m_required[%(count(jacobian))];	// true if matrix entry is not constant zero for current set of parameters (s. guesstop)
	// todo: distinguish between static and dynamic (tr-) matrix entries, the static ones must be loaded only once

	double m_entries[%(count(jacobian))];		// also used for static-entries in ac-analysis
	double m_entries_old[%(count(jacobian))];	// also used for dynamic-entries in ac-analyis
	double m_entries_kept[%(count(jacobian))]; // keep_ic, uic etc.

	enum EStates {
    <admst:variable name="numddt" select="0"/>
    <admst:for-each select="assignment">
      <admst:apply-templates select="rhs/tree" match="ddtcollect"/>
    </admst:for-each>
    <admst:for-each select="contribution">
      <admst:apply-templates select="rhs/tree" match="ddtcollect"/>
    </admst:for-each>
    <admst:for-each select="$(ddts)">
		_ddt_%(adms/id(.)), // %(.)
    </admst:for-each>
		_num_states = %(count($ddts)) // $(_num_states)
	};

	//	unsigned num_ddt(){return $(_num_states);}

	double* _states[OPT::_keep_time_steps]; // for each step: q0,i0, q1,i1, q2,i2, etc.
	double* _states_q1;	// for step t1:  q0, q1, q2, q3, etc.
	double* _states_kept;

//even states: _y
// odd states: _i

	double _DDT(double qq, unsigned stateno){
		unsigned qcap = stateno*2;
		trace4("$(DEV_NAME)::_DDT", stateno, qq, _sim->_mode, _sim->_phase);
<admst:if test="[$(_num_states)!=0]">
		trace2(&quot;$(DEV_NAME)::_DDT&quot;, uic_now(stateno), _sim-&gt;uic_now());
		if(uic_now(stateno)){
			_states[0][qcap] = qq;
#ifdef HAVE_KEEPCOEFF
			return ( qq + _states[1][stateno*2] ) * OPT::keepcoeff;
#endif
			incomplete(); // this won't work
		}
</admst:if>
		_states[0][stateno*2] = qq;

		if (_sim->analysis_is_ac()) {
			incomplete(); // shouldnt call DDT at all!
			return 0;
		}
		assert(_dt != 0.0);
		METHOD method;
		if( (_time[1]!=0) &amp;&amp; (_method_a==mTRAP) ) method = mTRAP;
		else method = mEULER;

		std::valarray&lt;FPOLY1&gt; q(OPT::_keep_time_steps);
		std::valarray&lt;FPOLY1&gt; i(OPT::_keep_time_steps);

		unsigned k = OPT::_keep_time_steps;
		// stupid copy. fixme.
		// gnucap stores things more usable.

		for (unsigned ii = 0; ii &lt; k; ii++)
		{
			trace5(&quot;$(DEV_NAME)::_DDT&quot;, ii, _states[ii][qcap], _states[ii][qcap+1], _time[0], _time[1] );
			assert( _states[ii][qcap] == _states[ii][qcap] );
			assert( _states[ii][qcap+1] == _states[ii][qcap+1] );
			q[ii] = FPOLY1(NOT_VALID, _states[ii][qcap], 0.0);
			i[ii] = FPOLY1(NOT_VALID, _states[ii][qcap+1], 0.0);// not necessary for ii==0
		}

		trace0("calling differentiate");
		i[0] = differentiate(&amp;q[0], &amp;i[0], _time, method);
		trace6("differentiate", _time[0], method, q[0].f0, q[1].f0, i[0].f0, _time[0]-_time[1]);
		trace1("differentiate", ( q[0].f0 - q[1].f0 ) / ( _time[0]-_time[1]));
		assert(i[0].f0 == i[0].f0);
		assert(i[0].f1 == i[0].f1);

		_states[0][qcap+1] = i[0].f0;

		assert(is_number(i[0].f0));
		if( _sim->analysis_is_static() ) {
			assert(i[0].f0 == 0.0);
		}

		return i[0].f0;
	}

	// noise: something is calculated, but never used
	<admst:text select="$fnoise/item" format="\n\tdouble fpnoise%(index($fnoise/item,.)), fenoise%(index($fnoise/item,.));\n"/>
	<admst:text select="$tnoise/item" format="\n\tdouble tnoise%(index($tnoise/item,.));\n"/>
        <admst:text select="$wnoise/item" format="\n\tdouble wnoise%(index($wnoise/item,.));\n"/>

	public:
		// initial conditions
		<admst:for-each select="contribution[is_ic_contribution(.)/[name='x']/value='yes']">
		// PARAMETER&lt;double&gt; ic_%(./lhs/branch/pnode)_%(./lhs/branch/nnode);
		</admst:for-each>
<admst:if test="[$(_num_states)!=0]">
	private: // from storage
		double tr_c_to_g(double c, double g)const;
		bool uic_now(unsigned stateno=-1)const;
</admst:if>
}; // DEV
// ------------------------------------------------------------------------ //
// ------------------------------------------------------------------------ //
class $(COMMON_NAME) : public COMMON_COMPONENT
{
	friend class $(DEV_NAME);
  private:
//    static map&lt;string, PARA_BASE $(COMMON_NAME)::*&gt; param_dict;
	public:
		void expand(const COMPONENT*);
		void precalc_first(const CARD_LIST*);
		void precalc_last(const CARD_LIST*);
	explicit $(COMMON_NAME)(const $(COMMON_NAME)&amp; p);
	explicit $(COMMON_NAME)(int c=0);
	~$(COMMON_NAME)();
	bool operator==(const COMMON_COMPONENT&amp;)const;
	COMMON_COMPONENT* clone()const {return new $(COMMON_NAME)(*this);}
	std::string name()const {itested();return &quot;$(module)&quot;;}
	void set_param_by_name(std::string name, std::string value);
<!--
	//void     set_param_by_index(int, std::string&, int);
	//bool     param_is_printable(int)const;
	//std::string param_name(int)const;
	//std::string param_name(int,int)const;
	//std::string param_value(int)const;
	//int param_count()const {return (8 + COMMON_COMPONENT::param_count());}
<!-- -->
	public: // initial conditions
  <admst:for-each select="contribution[is_ic_contribution(.)/[name='x']/value='yes']">
		PARAMETER&lt;double&gt; ic_%(./lhs/branch/pnode)_%(./lhs/branch/nnode);
  </admst:for-each>

<admst:if test="[$(_num_states)!=0]">
	public: // storage related

		void set_states(double*x){
		unsigned i=0;
		<admst:for-each select="$(ddts)">
			trace1(&quot;$(COMMON_NAME)::set_states ic_state_%(adms/id(.))&quot;, x[i]);
			ic_state_%(adms/id(.)) = x[i++];
		</admst:for-each>
		}
	
		// PARAMETER&lt;vector&lt;PARAMETER&lt;double &gt; &gt; &gt; state_ic;
  <admst:for-each select="$(ddts)">
		PARAMETER&lt;double&gt; ic_state_%(adms/id(.)); // %(.)
  </admst:for-each>
</admst:if>

	public: // global_instance input

<admst:for-each select="variable[input='yes' and scope='global_instance']">
  <admst:choose>
    <admst:when test="[type='real']">
      <admst:text format="\t\tPARAMETER&lt;double&gt; $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:when test="[type='integer']">
      <admst:text format="\t\tPARAMETER&lt;int&gt; $(var_prefix)%(.);\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="type of verilog-parameter is not real or interger\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:for-each>

};

//inline double $(DEV_NAME)::BP2( $(DEV_NAME)::ENodes p,  $(DEV_NAME)::ENodes n  ) const{
//	assert(is_number( _n[p].v0() - _n[n].v0() ));
//	return  min( 100. , max (  _n[p].v0() - _n[n].v0(), -100. ));
//}
</admst:for-each>

#endif
// vim should work with any ts=sw

</admst:open>


<!-- *********************************************************************************** -->
<!-- *********************************** Make d_*.cc *********************************** -->
<!-- *********************************************************************************** -->

<admst:template match="ddtcollect">
  <admst:message test="[/dbg_xml='yes']" format="*ddtcollect*\n"/>
  <admst:if test="[(adms/datatypename)='function' and name='ddt']">
    <admst:variable name="numddt" select="$numddt+1"/>
    <admst:push into="$ddts" select="." onduplicate="ignore"/>
  </admst:if>

  <admst:if test="[datatypename='mapply_unary'
    or datatypename='mapply_binary' or datatypename='mapply_ternary']">
    <admst:apply-templates select="arg1" match="ddtcollect"/>
  </admst:if>
  <admst:if test="[datatypename='mapply_binary' or datatypename='mapply_ternary']">
    <admst:apply-templates select="arg2" match="ddtcollect"/>
  </admst:if>
  <admst:if test="[datatypename='mapply_ternary']">
    <admst:apply-templates select="arg3" match="ddtcollect"/>
  </admst:if>
  <!-- <admst:if test="[datatypename='variable']">
    <admst:message format="fo2o %(name)"/>
    <admst:apply-templates select="module/assignment[lhs=%(.)]" match="ddtcollect"/>
  </admst:if> -->
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="e0">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)e*\n"/>
  <admst:apply-templates select="." match="%(adms/datatypename)">
    %(/simulator/tmp)
  </admst:apply-templates>
</admst:template>
<!-- -------------------------------------------------- -->
<!-- -------------------------------------------------- -->
<admst:open file="$(filename).cc">
/* AUTOMATICALLY GENERATED BY ADMS-XML
 * brought to you by Andreas Froese,
 *                   Felix Salfelder
 * do not edit
 */
#include &lt;limits&gt;
#include &lt;map&gt;
#include &lt;boost/assign.hpp&gt;
#include &lt;boost/algorithm/string.hpp>

#include &quot;$(filename).h&quot;
#include &quot;e_adms.h&quot;

using namespace std;
<admst:for-each select="/module[extern='no']">

//analog function stuff

  <admst:value-of select="attribute[name='ngspicename']/value"/>
  <admst:variable name="module" select="%s"/>
  <admst:variable name="DEV_NAME" select="DEV_%(upper-case(name))"/>
  <admst:variable name="MODEL_NAME" select="MODEL_%(upper-case(name))"/>

static $(DEV_NAME) p0;
static DISPATCHER&lt;CARD&gt;::INSTALL d0(&amp;device_dispatcher, DEVICE_TYPE, &amp;p0);

static $(MODEL_NAME) p1(&amp;p0);
static DISPATCHER&lt;MODEL_CARD&gt;::INSTALL d1(&amp;model_dispatcher, MODEL_TYPE, &amp;p1);

static $(COMMON_NAME) Default_$(COMMON_NAME)(CC_STATIC);

// ********************************* $(MODEL_NAME) ********************************* //
$(MODEL_NAME)::$(MODEL_NAME)(const $(DEV_NAME)* p) :
    MODEL_CARD(p),
  <!-- output='yes' - no matter, parametertype=instance or module - no matter also -->
  <admst:for-each select="variable[input='yes' and scope='global_model']">
	  <!--    <admst:apply-templates select="default" match="expression:stringify:noprobe"/> -->
	 $(var_prefix)%(name)(%(default/tree)),
  </admst:for-each>
    _key(&quot;dev_$(module)&quot;),
    _level()
{
}
// -------------------------------------------------------------------------- //
$(MODEL_NAME)::$(MODEL_NAME)(const $(MODEL_NAME)&amp; p) :
    MODEL_CARD(p),
  <admst:for-each select="variable[input='yes' and scope='global_model']">
    $(var_prefix)%(name)(p.$(var_prefix)%(name)),
  </admst:for-each>
    _key(p._key),
    _level(p._level)
{
}
// -------------------------------------------------------------------------- //
$(MODEL_NAME)::~$(MODEL_NAME)()
{
	incomplete();
}
// -------------------------------------------------------------------------- //
int $(MODEL_NAME)::param_count() const
{
	return(%(count(variable[input='yes'])));
}
// -------------------------------------------------------------------------- //
bool $(MODEL_NAME)::param_is_printable(int) const
{
	return(true);
}
// -------------------------------------------------------------------------- //
std::string $(MODEL_NAME)::param_name(int i)const
{
	switch(i)
	{
    <admst:for-each select="variable[input='yes']">
		case %(position(.)-1): return(&quot;%(name)&quot;);
    </admst:for-each>
		default: return(&quot;&quot;);
	}
	return(&quot;&quot;);
}
// -------------------------------------------------------------------------- //
std::string $(MODEL_NAME)::param_name(int i, int j)const
{
	if( j == 0 ) return(param_name(i));
	return(&quot;&quot;);
}
// -------------------------------------------------------------------------- //
std::string $(MODEL_NAME)::param_value(int i) const
{
	switch(i)
	{
    <admst:for-each select="variable[input='yes' and scope='global_model']">
		case %(position(.)-1): return($(var_prefix)%(name).string());
    </admst:for-each>
		default: return(&quot;&quot;);
	}
	return(&quot;&quot;);
}
// -------------------------------------------------------------------------- //
void $(MODEL_NAME)::set_param_by_index(int i, std::string&amp; value, int offset)
{
	USE(value);
	switch(i)
	{
    <admst:for-each select="variable[input='yes' and scope!='global_instance']">
		case %(position(.)-1): $(var_prefix)%(name)=value; break; // %(scope)
    </admst:for-each>
		default: throw Exception_Too_Many(i, %(count(variable[input='yes'])), offset); break;
	}
}
// -------------------------------------------------------------------------- //
#ifdef HAVE_PARA_BASE
std::map&lt;string, PARA_BASE $(MODEL_NAME)::*&gt; $(MODEL_NAME)::param_dict
<admst:if test="[count(variable[input='yes' and scope!='global_instance'])!=0]">
  = boost::assign::map_list_of
<admst:join select="variable[input='yes' and scope!='global_instance']" separator="">
  (&quot;%(name)&quot;, (PARA_BASE $(MODEL_NAME)::*)  (&amp;$(MODEL_NAME)::$(var_prefix)%(name)))
</admst:join>
</admst:if>;

std::map&lt;string, PARA_BASE $(MODEL_NAME)::*&gt; $(MODEL_NAME)::param_dict_low
<admst:if test="[count(variable[input='yes' and scope!='global_instance'])!=0]">
  = boost::assign::map_list_of
<admst:join select="variable[input='yes' and scope!='global_instance']" separator="">
    (toLower("%(name)"), (PARA_BASE $(MODEL_NAME)::*)
        (&amp;$(MODEL_NAME)::$(var_prefix)%(name)))
</admst:join>
</admst:if>;
#endif
// -------------------------------------------------------------------------- //
void $(MODEL_NAME)::set_param_by_name(string Name, string Value)
{
#ifdef HAVE_PARA_BASE
	PARA_BASE $(MODEL_NAME)::* x = (OPT::case_insensitive)?
		(param_dict_low[toLower(Name)]) : (param_dict[Name]);
	trace3("$(MODEL_NAME)::set_param_by_name", Name, OPT::case_insensitive, x);
	if(x) {
		PARA_BASE* p = &(this->*x);
		*p = Value;
		return;
	}
	throw Exception_No_Match(Name);
#else
	// fall back to linear search
	MODEL_CARD::set_param_by_name(Name,Value);
#endif
}
// -------------------------------------------------------------------------- //
void $(COMMON_NAME)::precalc_first(const CARD_LIST* par_scope)
{
	USE(par_scope);
  <admst:for-each select="variable[input='yes' and scope='global_instance']">
    <admst:apply-templates select="default/tree" match="subexpression:stringify:noprobe"/>
    <admst:choose>
      <admst:when test="[type='real']">
        <admst:text format="\te_val(&amp;($(var_prefix)%(name)), (double)%s, par_scope);\n"/>
      </admst:when>
      <admst:when test="[type='integer']">
        <admst:text format="\te_val(&amp;($(var_prefix)%(name)), (int)%s, par_scope);\n"/>
      </admst:when>
    </admst:choose>
  </admst:for-each>
  <admst:for-each select="$(ddts)">
	ic_state_%(adms/id(.)).e_val(NOT_INPUT,par_scope);
  </admst:for-each>
}
// -------------------------------------------------------------------------- //
void $(MODEL_NAME)::precalc_first()
{
	MODEL_CARD::precalc_first();
// -------------------------------------------------------------------------- //
	const CARD_LIST* par_scope = scope();
	assert(par_scope);
	MODEL_CARD::precalc_first();
  <admst:for-each select="variable[input='yes' and scope!='global_instance']">
    <admst:apply-templates select="default/tree" match="subexpression:stringify:noprobe"/>
    <admst:choose>
      <admst:when test="[type='real']">
        <admst:text format="\te_val(&amp;(this->$(var_prefix)%(name)), (double)%s, par_scope);\n"/>
      </admst:when>
      <admst:when test="[type='integer']">
        <admst:text format="\te_val(&amp;(this-&gt;$(var_prefix)%(name)), (int)%s, par_scope);\n"/>
      </admst:when>
    </admst:choose>
  </admst:for-each>

	// @(initial_model) block
	$(MODEL_NAME) *m = this; USE(m);
   <admst:apply-templates select="analog" match="analog:initial_model" required="yes"/>
}
// -------------------------------------------------------------------------- //
void $(MODEL_NAME)::set_dev_type(const std::string&amp; new_type)
{
	trace1("$(MODEL_NAME)::set_dev_type", new_type);

	_key = new_type;
	if (OPT::case_insensitive)
		notstd::to_lower(&amp;_key);
}
// ********************************* $(DEV_NAME) ********************************* //
$(DEV_NAME)::$(DEV_NAME)()
  :ADMS_BASE(),
   _nodes(),
	_states_q1(NULL)
<admst:for-each select="instance">,
    _%(name)(0)
</admst:for-each>
{
	std::fill_n(m_entries, %(count(jacobian)), 0);
	std::fill_n(m_entries_old, %(count(jacobian)), 0);

	attach_common(&amp;Default_$(COMMON_NAME));

	_n = _nodes;
	uint_t i = 0;
<admst:for-each select="node[grounded='no']" separator=";">
	// %(name) %(discipline) %(location)
  <admst:text test="[discipline/name='degradational']" format="\n\t_n[i].set_adp();"/>
	assert(!(_n[i++].n_()));
</admst:for-each>

	for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii){
		// not here!
		_states[ii] = NULL;
	}
}
// -------------------------------------------------------------------------- //
$(DEV_NAME)::$(DEV_NAME)(const $(DEV_NAME)&amp; p) :
    ADMS_BASE(p),
    _nodes(),
    _states_q1(NULL)
<admst:for-each select="instance">,
    _%(name)(0)
</admst:for-each>
{
	std::fill_n(m_entries, %(count(jacobian)), 0);
	std::fill_n(m_entries_old, %(count(jacobian)), 0);

	for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii){
		// not here!
		_states[ii] = NULL;
	}

	_n = _nodes;
	for (uint_t ii = 0; ii &lt; int_nodes() + max_nodes(); ++ii){
		_n[ii] = p._n[ii];
		trace3("$(DEV_NAME)::$(DEV_NAME)", ii, _n[ii].is_adp(), p._n[ii].is_adp());
	}

	for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii){
		_states[ii] = NULL;
	}
}
// -------------------------------------------------------------------------- //
<admst:if test="[$(_num_states)!=0]">
bool $(DEV_NAME)::uic_now(unsigned stateno)const{
	if(!_sim-&gt;uic_now()) return 0;
	const $(COMMON_NAME)* c = static_cast&lt;const $(COMMON_NAME)*&gt;(common());

	incomplete();
	switch(stateno){
<admst:for-each select="$(ddts)">
		case _ddt_%(adms/id(.)):
			return c-&gt;ic_state_%(adms/id(.)).has_good_value();
</admst:for-each>
		default:
			assert(0); return(0);
	}
}
</admst:if>
// -------------------------------------------------------------------------- //
$(DEV_NAME)::~$(DEV_NAME)()
{
<admst:if test="[$(_num_states)!=0]">
	if (_states[0]) {
		trace1(&quot;cleaning up states&quot;, long_label());
		for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii) {
			trace1(&quot;cleaning up states&quot;, ii);
			assert(_states[ii]);
			delete [] _states[ii];
			_states[ii] = 0;
		}

		assert(_states_q1);
		delete[] _states_q1;

	} else {
    //BUG! OPT::_keep_time_steps could have changed....
    for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii)
      assert(!_states[ii]);
	assert(!_states_q1);
	}
</admst:if>
}
// ---------------- all analog functions from verilog-module here: ---------- //

<admst:apply-templates select="/module" match="analogfunctionC"/>
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::set_param_by_name(std::string Name, std::string value)
{
	COMPONENT::set_param_by_name(Name,value);
}
// -------------------------------------------------------------------------- //
void $(COMMON_NAME)::expand(const COMPONENT* d)
{
	trace2(&quot;$(COMMON_NAME)::expand&quot;, modelname(), OPT::_keep_time_steps);

	COMMON_COMPONENT::expand(d);
	attach_model(d);
	const $(MODEL_NAME)* m = dynamic_cast&lt;const $(MODEL_NAME)*&gt;(model());
	if (!m) {
		throw Exception_Model_Type_Mismatch(d-&gt;long_label(), modelname(), &quot;$(module)&quot;);
	}
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::expand()
{
	trace2("$(DEV_NAME)::expand()", long_label(), common()->modelname());
	ADMS_BASE::expand();
	
	const $(COMMON_NAME)* c = static_cast&lt;const $(COMMON_NAME)*&gt;(common());
	assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m);
	trace2(&quot;$(DEV_NAME)::expand&quot;, c-&gt;modelname() ,OPT::_keep_time_steps);

	if (!subckt()) {
		new_subckt();
	}else{
	}

  // set m_required
  guesstopology();
  // allocate state vectors

<admst:if test="[$(_num_states)!=0]">
  if (!_states[0]) {
    // initialise in tr_begin() ?
    for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii)  // &gt;
    {
      assert(!_states[ii]);
      _states[ii] = new double[_num_states*2];
    }
    //assert(!_states_1);
    //_states_1 = new double[_num_states*2];

	assert(!_states_q1);
	_states_q1 = new double[_num_states];


  } else {
    for (int ii=0; ii&lt;OPT::_keep_time_steps; ++ii) //&gt;
            assert(_states[ii]);
//    assert(_states_1);
	assert(_states_q1);
  }
</admst:if>

	std::fill_n(it1, $_num_stamp_nodes, 0);

	// fix up internal nodes
	if (_sim->is_first_expand()) {
		precalc_first();
		precalc_last();
		trace4("$(DEV_NAME)::expand nodes", net_nodes(), max_nodes(), matrix_nodes(), int_nodes());
<admst:for-each select="node[grounded='no']" separator=";">
		// %(name) %(discipline) %(location)
</admst:for-each>

		for(uint_t ii = net_nodes(); ii &lt; matrix_nodes(); ++ii) {
#ifdef HAVE_ADP_NODE
			trace2("$(DEV_NAME)::expand nodes", ii, _n[ii].is_adp());
			if(!_n[ii].is_adp()){
#endif
				if( !_n[ii].n_() ) {
					if( ii &lt; max_nodes() ) {
						trace2("$(DEV_NAME)::expand node between _net_nodes and max_nodes(): set to ground", ii, node_name(ii));
						_n[ii].set_to_ground(this);
					} else {
						trace2("$(DEV_NAME)::expand node between max_nodes() and matrix_nodes(): alloc new", ii, node_name(ii));
						_n[ii].new_model_node(node_name(ii), this);
					}
				}
#ifdef HAVE_ADP_NODE
			} else if ( !_n[ii].a_()){
				trace0("$(DEV_NAME)::expand node, registering new adpnode");
				// no external adp nodes, so no unconnected ones.
				untested();
				_n[ii].new_model_adp_node(node_name(ii), this);
				assert(_n[ii].t_()!=INVALID_NODE);
			} else {
				trace1("$(DEV_NAME)::expand node,not registering new adpnode", _n[ii].a_()->long_label() );

			}
#endif
		}

		assert(max_nodes() != 0);
		for (uint_t ii = 0; ii &lt; matrix_nodes(); ++ii)
		{
			trace2("$(DEV_NAME)::expand matrix node check", ii, _n[ii].is_adp());
			assert(_n[ii].t_()!=INVALID_NODE);
#ifdef HAVE_ADP_NODE
			assert(_n[ii].n_() || _n[ii].is_adp());
#else
			assert(_n[ii].n_());
#endif
		}
<admst:for-each select="instance">
		if (!_%(name)) {
			trace1("$(DEV_NAME)::expand, instanciating", hp(_%(name)));
			// FIXME: recycle u_lang or move to e_adms.cc
			const CARD* p = LANGUAGE::find_proto(&quot;%(module)&quot;, this);
			if(!p) throw Exception("cannot find '%(module)'. Load module first?");
			_%(name) = dynamic_cast&lt;COMPONENT*&gt;(p-&gt;clone_instance());
			trace2("$(DEV_NAME)::expand", long_label(), _%(name)->has_common());
			if(_%(name)->has_common()){
				trace2("$(DEV_NAME)::expand", long_label(),_%(name)->mutable_common()->model());

				if(dynamic_cast&lt;const MODEL_CARD*&gt;(p)){
					untested();
					trace2("$(DEV_NAME)::expand attaching submodel", _%(name)->common()->modelname(), _%(name)->mutable_common()->model());
					_%(name)->set_dev_type("%(module)");
					_%(name)->mutable_common()->attach_model(_%(name));

				}
			}

			assert(_%(name));
			subckt()-&gt;push_front(_%(name));
		}
		{
			node_t nodes[] = {
  <admst:join select="terminal" separator=",">
				_n[n_%(nodefrominstantiator) /*%(nodefrommodule)*/]
  </admst:join> };
			trace2("expand", long_label(),_%(name)->mutable_common()->model());
			_%(name)->set_parameters("%(name)", this, _%(name)->mutable_common(), 0/*value*/, 0, NULL, %(count(terminal)), nodes);
			trace2("expand", long_label(),_%(name)->mutable_common()->model());

			PARAMETER&lt;double&gt;tmp; USE(tmp);
  <admst:for-each select="parameterset">
			tmp = <admst:apply-templates select="value/tree" match="e0"/>;
			_%(../name)-&gt;set_param_by_name(&quot;%(parameter)&quot;, tmp.string());
  </admst:for-each>
		}
</admst:for-each>
	}else{
		untested();
	}
  <admst:if test="[$has_sckt!=0]">
	assert(subckt());
	trace0("$(DEV_NAME)::expand: expanding netlist");
	subckt()->expand();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::precalc_first()
{
	ADMS_BASE::precalc_first();
  <admst:if test="[$has_sckt!=0]">
	if(subckt()) subckt()-&gt;precalc_first();
  </admst:if>
// 	&lt;admst:for-each select=&quot;variable[input='yes']&quot;&gt;
// 	 %(name).e_val(0,scope());
// 	&lt;/admst:for-each&gt;
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::precalc_last()
{
	ADMS_BASE::precalc_last();
  <admst:if test="[$has_sckt!=0]">
	if(subckt()) subckt()-&gt;precalc_last();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(COMMON_NAME)::precalc_last(const CARD_LIST* scope)
{
	COMMON_COMPONENT::precalc_last(scope);
  <admst:for-each select="contribution[is_ic_contribution(.)/[name='x']/value='yes']">
	ic_%(./lhs/branch/pnode)_%(./lhs/branch/nnode).e_val(0,scope);
  </admst:for-each>

  <admst:apply-templates select="." match="model"/>
#if 0 // hmm wrong place
  <admst:apply-templates select="analog" match="analog:initial_instance" required="yes"/>
#endif
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::internal_precalc()
{
  trace1(&quot; $(DEV_NAME)::internal_precalc()&quot;, long_label());

#define _STATIC
#define _DYNAMIC
  <admst:apply-templates select="analog" match="analog:initial_instance" required="yes"/>
#undef _STATIC
#undef _DYNAMIC

	set_converged();
	untested();
}
// -------------------------------------------------------------------------- //
#define tr_call_iwants(m_1, m_2) _sim->_aa.iwant(_n[m_1].m_(), _n[m_2].m_()); _sim->_lu.iwant(_n[m_1].m_(), _n[m_2].m_())
// -------------------------------------------------------------------------- //

void $(DEV_NAME)::tr_iwant_matrix()
{
	trace0("$(DEV_NAME)::tr_iwant_matrix()");
	assert(is_device());
  <admst:if test="[$has_sckt!=0]">
	if (subckt()) {
      subckt()->tr_iwant_matrix();
	}else{ // untested();
	}
  </admst:if>
  //	ADMS_BASE::tr_iwant_matrix();
  <admst:for-each select="jacobian">
    <admst:variable  name="myF_par1" select="%(row/name)"/>
    <admst:variable  name="myF_par2" select="%(column/name)"/>
    <admst:apply-templates select=".." match="utils_compare_nodes_position"/>
    <!-- if bigger check if 'smaller' exists -->
    <admst:if test="[$myF_ret1='bigger']">
      <admst:variable  name="local_ifelse" select="else"/>
      <admst:if test="(../jacobian[row/name=$myF_par2 and column/name=$myF_par1])">
        <admst:variable  name="local_ifelse" select="if"/>
      </admst:if>
      <admst:if test="[$local_ifelse='else']">
        <admst:variable  name="myF_ret1" select="smaller"/>
        // COUNTER DOESNT EXIST
      </admst:if>
    </admst:if>
    <!-- TODO if there are only node1:node1 jacobian (without any other combinations),
         it will not call iwant matrix
    -->
    <admst:if test="[$myF_ret1='smaller']">
      <!-- only if 1st node has smaller position than 2nd -->
      <admst:text format="\n\tif( m_required[m_$(myF_par1)_$(myF_par2)]"/>
      <admst:text format=" || m_required[m_$(myF_par2)_$(myF_par1)]"
        test="../jacobian[row/name=$myF_par2 and column/name=$myF_par1]"/>
      <admst:text format=" || m_required[m_$(myF_par1)_$(myF_par1)]"
        test="../jacobian[row/name=$myF_par1 and column/name=$myF_par1]"/>
      <admst:text format=" || m_required[m_$(myF_par2)_$(myF_par2)]"
        test="../jacobian[row/name=$myF_par2 and column/name=$myF_par2]"/>)
	  //tr_call_iwants(_N($(myF_par1)).m_(), _N($(myF_par2)).m_());
#ifdef HAVE_ADP_NODE
		if (!_n[$(node_prefix)$(myF_par1)].is_adp() && !_n[$(node_prefix)$(myF_par2)].is_adp())
#endif
		{
			trace2("tr_iwants $(myF_par1) $(myF_par2)", _n[$(node_prefix)$(myF_par1)].is_adp(), _n[$(node_prefix)$(myF_par2)].is_adp() );
			tr_call_iwants($(node_prefix)$(myF_par1), $(node_prefix)$(myF_par2));
		}
    </admst:if>
  </admst:for-each>
	//tr_iwant_matrix_extended();
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::ac_iwant_matrix()
{
	trace0("$(DEV_NAME)::ac_iwant_matrix()");
  <admst:for-each select="jacobian">
    <admst:variable  name="myF_par1" select="%(row/name)"/>
    <admst:variable  name="myF_par2" select="%(column/name)"/>
    <admst:apply-templates select=".." match="utils_compare_nodes_position"/>
    <!-- if bigger check if 'smaller' exists -->
    <admst:if test="[$myF_ret1='bigger']">
      <admst:variable  name="local_ifelse" select="else"/>
      <admst:if test="(../jacobian[row/name=$myF_par2 and column/name=$myF_par1])">
        <admst:variable  name="local_ifelse" select="if"/>
      </admst:if>
      <admst:if test="[$local_ifelse='else']">
        <admst:variable  name="myF_ret1" select="smaller"/>
        // COUNTER DOESNT EXIST
      </admst:if>
    </admst:if>
    <!-- TODO if only equal exists (without any other combinations), it will not call iwant matrix -->
    <admst:if test="[$myF_ret1='smaller']">
      <!-- only if 1st node has smaller position than 2nd -->
      <admst:text format="\n\tif( m_required[m_$(myF_par1)_$(myF_par2)]"/>
      <admst:text format=" || m_required[m_$(myF_par2)_$(myF_par1)]"
        test="../jacobian[row/name=$myF_par2 and column/name=$myF_par1]"/>
      <admst:text format=" || m_required[m_$(myF_par1)_$(myF_par1)]"
        test="../jacobian[row/name=$myF_par1 and column/name=$myF_par1]"/>
      <admst:text format=" || m_required[m_$(myF_par2)_$(myF_par2)]"
        test="../jacobian[row/name=$myF_par2 and column/name=$myF_par2]"/>)
#ifdef HAVE_ADP_NODE
		if (!_n[$(node_prefix)$(myF_par1)].is_adp() && !_n[$(node_prefix)$(myF_par2)].is_adp())
#endif
		_sim->_acx.iwant(_N($(myF_par1)).m_(), _N($(myF_par2)).m_());
    </admst:if>
  </admst:for-each>
	//ac_iwant_matrix_extended();
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_begin()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->tr_begin();
  </admst:if>
	const $(COMMON_NAME)* c = prechecked_cast&lt;const $(COMMON_NAME)*&gt;(common());
	assert(c); USE(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	assert(m);
	USE(m);
	
	trace0(&quot;$(DEV_NAME)::tr_begin&quot;);
	ADMS_BASE::tr_begin();

	// initial step

	<admst:apply-templates select="analog" match="analog:initial_step"/>
	// /initial step

	assert($_num_states == _num_states);
	for(int ii=0; ii&lt;OPT::_keep_time_steps; ii++){
		std::fill_n(_states[ii], ($_num_states)*2, 0);
	}
	std::fill_n(_states_q1, ($_num_states), 0);

<admst:if test="[$_num_states!=0]">
	unsigned k = 0;
  <admst:for-each select="$(ddts)">
	if(c-&gt;ic_state_%(adms/id(.)).has_good_value()){
		trace1(&quot;$(DEV_NAME)::tr_begin&quot;, c-&gt;ic_state_%(adms/id(.)));
		_states[1][2*k++] = c-&gt;ic_state_%(adms/id(.));
	}
  </admst:for-each>
</admst:if>

	// like d_vs.
	// &lt;admst:text format=&quot;\\n\\tfor(int i=0; i&lt;$_num_stamp_nodes; i++)&quot;/&gt;
	for(int i=0; i&lt;$_num_stamp_nodes; i++)
		it0[i] = it1[i]=0;

        //     in storag:
        // _m1 = _m0 = CPOLY1(0., 0., 0.);
	std::fill_n(m_entries, %(count(jacobian)), 0);
	std::fill_n(m_entries_old, %(count(jacobian)), 0);

	internal_precalc();

	// load static ??
        //&lt;admst:text format=&quot;\\n\\tfor(int i=0; i&lt;$_num_stamp_nodes; i++)&quot;/&gt;
        //	tr_load_source_point(_n[i], &amp;(it0[i]), &amp;(it1[i]));
        //

#ifdef load_in_begin

	<admst:for-each select="jacobian[static='yes']">
	<admst:value-of select="column/name"/>
	<admst:value-of select="row/name"/>
	<admst:value-of select="column/name"/>
	<admst:value-of select="row/name"/>
	<admst:value-of select="column/name"/>
	<admst:value-of select="row/name"/>
	<admst:value-of select="column/name"/>
	<admst:value-of select="row/name"/>
	if( m_required[m_%s_%s] )
		tr_load_point(_N(%s), _N(%s), &amp;(m_entries[m_%s_%s]), &amp;(m_entries_old[m_%s_%s]));
  </admst:for-each>
#endif

        q_load();

}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::dc_advance()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;dc_advance();
  </admst:if>
	trace0(&quot;$(DEV_NAME)::dc_advance&quot;);

	// ELEMENT:: does nothing

	ADMS_BASE::dc_advance();

<admst:if test="[$_num_states!=0]">
	// storage-like
	for (int t = 1;  t &lt; OPT::_keep_time_steps;  ++t) {
		//  notstd::copy_n(_states[0], 2*(unsigned)_num_states, _states[i]);
		//_i[t] = _i[0]; &lt;= what STORAGE does

		for( unsigned i=0; i&lt;_num_states; ++i){
			_states[t][2*i+1] = _states[0][2*i+1];
		}
	}
</admst:if>
	// ??
	//		for( unsigned i=0; i&lt;_num_states; ++i){
	//			_states_1[2*i+1]=_states[0][2*i+1];
	//		}
	//

}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tt_accept()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->tt_accept();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::stress_apply()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->stress_apply();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tt_advance()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->tt_advance();
  </admst:if>
	trace3("$(DEV_NAME)::tt_advance", short_label(), _sim->_time0, _sim->_dt0);
	if (_time[0] > _sim->_time0) {untested();
		for (int i=0; i&lt;OPT::_keep_time_steps-1; ++i) {
			_time[i] = _time[i+1];
			//_y[i] = _y[i+1];
			for( unsigned j=0; j&lt;_num_states; ++j){
				_states[i][2*j+1] = _states[i+1][2*j+1];
			}
		}
		_time[OPT::_keep_time_steps-1] = 0.;
		// _y[OPT::_keep_time_steps-1] = FPOLY1(0., 0., 0.);
	}

	for (int i=OPT::_keep_time_steps-1; i>=0; --i) {
		// FIXME: copy all timesteps to 0
		_time[i] = 0.;
		//    assert(_time[i] &lt; _time[i-1] || _time[i] == 0.);
	}
	_dt=0;
}
/*--------------------------------------------------------------------------*/
#ifdef HAVE_TT
void $(DEV_NAME)::tt_begin()
{
	trace1("$(DEV_NAME)::tt_begin", _sim->_Time0);
	ADMS_BASE::tt_begin();
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;tt_begin();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_stress_last()
{
	trace1("$(DEV_NAME)::tr_stress_last", _sim->_Time0);
	ADMS_BASE::tr_stress_last();
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;tr_stress_last();
  </admst:if>
}
#endif
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_restore()
{
	ADMS_BASE::tr_restore(); internal_precalc();
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;tr_restore();
  </admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_advance()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->tr_advance();
  </admst:if>
  trace3(&quot;DEV_LAMP::tr_advance&quot;, _sim-&gt;_time0, _time[0], _time[1]);
  ADMS_BASE::tr_advance(); // does time.

  // shift states, clear states[0]
  double* t = _states[OPT::_keep_time_steps-1];
  for (int ii = OPT::_keep_time_steps-1;  ii &gt; 0;  --ii){
    _states[ii] = _states[ii-1];
  }
  _states[0] = t;

  //initial guess.
  notstd::copy_n(_states[1], 2*(unsigned) _num_states, _states[0]);


// storage would do:
//for (int i=OPT::_keep_time_steps-1; i&gt;0; --i) {
//_i[i] = _i[i-1];
// }
//  element does:
//for (int i=OPT::_keep_time_steps-1; i&gt;0; --i) {
//  assert(_time[i] &lt; _time[i-1] || _time[i] == 0.);
//  _time[i] = _time[i-1];
//  _y[i] = _y[i-1];
//}
// _time[0] = _sim-&gt;_time0;
//
// _dt = _time[0] - _time[1];
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_regress()
{
	ADMS_BASE::tr_regress();
  <admst:if test="[$has_sckt!=0]">
	assert(subckt()); subckt()->tr_regress();
  </admst:if>
}
// -------------------------------------------------------------------------- //
bool $(DEV_NAME)::tr_needs_eval()const
{
	trace3(&quot;$(DEV_NAME)::tr_needs_eval&quot;, converged(), is_q_for_eval(), _sim-&gt;is_advance_iteration());
  <admst:if test="[$has_sckt!=0]">
	if( subckt()->tr_needs_eval()){
		return true;
	}
  </admst:if>
	if (is_q_for_eval()) return false;
	if (!converged()) return true;
	if (_sim-&gt;is_advance_iteration()) return true;
	if (_time[1] == 0){
		trace0(&quot;tr_needs_eval: //BUG// needed?? for ngspice jfet, but not for spice3f5 jfet&quot;);
		return true;
	}

	// check the node voltages, reference to ground
<admst:for-each select="node[grounded='no' and discipline/name!='degradational']">
	// probably simpler: if (!conchk(_sim->_n[ii]->vt1(), _n[ii]->v0(), 0, OPT::reltol*OPT::bypasstol))
	if (!conchk(_sim->_vt1[_n[$node_prefix%(name)].m_()], _n[$node_prefix%(name)].v0(), 0, OPT::reltol*OPT::bypasstol)) {
		return true;
	}
  <admst:variable name="name1" select="%(name)"/>
  <admst:variable name="upper" select="no"/>
  <admst:for-each select="../node[grounded='no' and discipline/name!='degradational']">
    <admst:if test="[$upper='yes']">
	// relative check %(name) $(name1) FIXME: iterate over branches?!
	if (!conchk((_sim->_vt1[_n[$node_prefix%(name)].m_()] - _sim->_vt1[_n[$node_prefix$(name1)].m_()]),
	(_n[$node_prefix%(name)].v0() - _n[$node_prefix$(name1)].v0()), 0, OPT::reltol*OPT::bypasstol)) {
		return true;
	}
    </admst:if>
    <admst:variable name="upper" select="yes" test="[name=$(name1)]"/>
  </admst:for-each>
</admst:for-each>

	trace0(&quot;$(DEV_NAME)::tr_needs_eval: no&quot;);
	return false;
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_eval_kept()
{
#if 0
<!-- <admst:apply-templates select="analog/code[datatypename='block']/item" match="tr_eval_kept" required="yes"/> -->
#endif // 0
}
// -------------------------------------------------------------------------- //
// MODEINITFLOAT = normal iteration
// MODEINITPRED  = 1st iter at a new time point
// MODEINITTRAN  = 1st iter at 1st time pt after initial DC
// MODEINITFIX   = like FLOAT, but honor options like &quot;off&quot;
// MODEINITJCT   = initial guess
// MODEINITSMSIG = like FLOAT, but setup for small signal, don't load arrays
// -------------------------------------------------------------------------- //
bool $(DEV_NAME)::do_tr()
{
	trace1(&quot;$(DEV_NAME)::do_tr&quot;, long_label());
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;do_tr();
  </admst:if>

	assert(OPT::_keep_time_steps &lt;= 8);
	assert(_dt == NOT_VALID || conchk(_dt, _time[0] - _time[1]));

	if (_sim-&gt;analysis_is_tran_dynamic())
	{
		trace1(&quot;$(DEV_NAME)::do_tr tran dynamic&quot;, long_label());
		tr_eval(); // if needed??
		q_load();

	} else {
		tr_eval(); // if needed??
		q_load();
	}

	// convergence check -- gnucap method
  <admst:if test="[$has_sckt!=0]">
	set_converged( subckt()-&gt;do_tr());
  </admst:if>
  <admst:if test="[$has_sckt=0]">
	set_converged(true);
  </admst:if>

<admst:if test="[not( count($ddts)=0) ]">
	for (unsigned ii = 0; ii &lt; _num_states; ++ii) { // similar e_cap.cc do_tr();
		set_converged(converged() &amp;&amp; conchk(_states[0][ii*2], _states_q1[ii]));	// was _states_1 originally
		//trace3(&quot;&quot;, ii, _states_1[ii], _states[0][ii]);
		//_states_1[ii] = _states[0][ii];
	}
</admst:if>

	for(int i=0; converged() &amp;&amp; i&lt;$_num_stamp_nodes; i++)
	{
		assert(is_number(it0[i]));
		bool conv= conchk(it0[i], it1[i]);
		if (!conv) {
			trace6(&quot;tr_eval: i not converged &quot;, i , it0[i], it1[i], _sim-&gt;iteration_number(), it0[i]-it1[i], conv );
			set_converged(false);
		}
	}


  <admst:text format="\n\tfor(int i=0; converged() && i<%(count(jacobian)); i++)"/>
	if (m_entries[i]){
		bool conv = conchk(m_entries[i], m_entries_old[i]);
		if (!conv) {
			trace4(&quot;tr_eval: m not converged&quot;, i , m_entries[i], m_entries_old[i], _sim-&gt;iteration_number() );
			set_converged(conv);
		}
	}

	bool needs_load = !converged();
	for (int i = 0; !needs_load &amp;&amp; i &lt; $_num_stamp_nodes; i++) {
		if (m_required[i]){
			needs_load = !conchk(it0[i], it1[i], 0, OPT::reltol*OPT::loadtol);
		}
	}

  for(int i=0; !needs_load &amp;&amp; i&lt;%(count(jacobian)); i++)
    if (m_required[i]){
      needs_load = !conchk(m_entries[i], m_entries_old[i], 0, OPT::reltol * OPT::loadtol);
    }

  store_values();

                // if (needs_load)
  trace2(&quot;$(DEV_NAME)::do_tr done &quot;, converged(), needs_load );
  return converged();
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::keep_ic()
{
#ifdef HAVE_KEEPCOEFF
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;keep_ic();
  </admst:if>
#endif

	//double x =  dn_diff(_n[OUT1].v0(),_n[OUT2].v0() );
	$(COMMON_NAME)* c = prechecked_cast&lt; $(COMMON_NAME)*&gt;(mutable_common());
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);

<admst:if test="[not($_num_states=0)]">
	tr_eval();
	store_values();
	double s[_num_states];
	for (unsigned i=0; i&lt;_num_states; i++){
		trace2(&quot;$(DEV_NAME)::keep_ic&quot;, long_label(), _states[0][i*2] );
		s[i]=_states[0][i*2];
	}
	c-&gt;set_states(s);
</admst:if>
	trace1(&quot;$(DEV_NAME)::keep_ic&quot;, long_label() );
	notstd::copy_n(it0, $_num_stamp_nodes , i_kept);
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::store_values()
{
  trace1(&quot;$(DEV_NAME)::store_values&quot;, long_label());
  // notstd::copy_n(_states[0],(unsigned) 2*_num_states, _states_1);

<admst:if test="[($_num_states!='0')]">
	for( unsigned i=0; i&lt;_num_states; ++i){
		_states_q1[i]=_states[0][2*i];
	}
</admst:if>

  // ELEMENT does y[0] -&gt; y_1
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_queue_eval()
{
	ADMS_BASE::tr_queue_eval();
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;tr_queue_eval();
	</admst:if>
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_load()
{
	trace5(&quot;$(DEV_NAME)::tr_load()&quot;, long_label(), _sim-&gt;is_inc_mode(),
	                   OPT::traceload, _loaditer, _sim-&gt;iteration_tag());

  <admst:if test="[$has_sckt!=0]">
	if( !OPT::traceload || !_sim-&gt;is_inc_mode()) {
		subckt()-&gt;tr_load();
	} else {
		// subdevices have queued themselves if necessary.
	}
	</admst:if>
	trace0(&quot;$(DEV_NAME)::tr_load&quot;);
          /// BUG: only dynamic!
	for(int i=0; i&lt;$_num_stamp_nodes; i++){
		trace3(&quot;$(DEV_NAME)::tr_load&quot;, i,  (it0[i]),  (it1[i]));
		tr_load_source_point(_n[i], &(it0[i]), &amp;(it1[i]));

		trace3("$(DEV_NAME)::tr_load done&quot;, &quot;i&quot;,     (it0[i]),  (it1[i]));
	}

          /// BUG: only load dynamic nodes! ?

  <admst:for-each select="jacobian">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
    if( m_required[m_%s_%s] ){
      trace4("$(DEV_NAME)::tr_load&quot;, &quot;%s&quot;, &quot;%s&quot;, m_entries[m_%s_%s], m_entries_old[m_%s_%s]);
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
      tr_load_point(_N(%s), _N(%s), &amp;(m_entries[m_%s_%s]), &amp;(m_entries_old[m_%s_%s]));
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
      trace4("$(DEV_NAME)::tr_load done &quot;, &quot;%s&quot;, &quot;%s&quot;, m_entries[m_%s_%s], m_entries_old[m_%s_%s]);

    }
  </admst:for-each>

  <admst:for-each select="jacobian">
	  <admst:value-of select="column/name"/>
	  <admst:value-of select="row/name"/>
	  <admst:value-of select="column/name"/>
	  <admst:value-of select="row/name"/>
	  <admst:value-of select="column/name"/>
	  <admst:value-of select="row/name"/>
	  if( !m_required[m_%s_%s] &amp;&amp; (m_entries[m_%s_%s]!=0) ){
		  error(bWARNING, "*** m_%s_%s\\n");
	  }

  </admst:for-each>

  _loaditer = _sim->iteration_tag();
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_unload()
{
	untested();
	incomplete();
	_sim-&gt;mark_inc_mode_bad();
	tr_load();
}
// -------------------------------------------------------------------------- //
TIME_PAIR $(DEV_NAME)::tr_review()
{
	assert(_dt == NOT_VALID || conchk(_dt, _time[0] - _time[1]));

	double timestep = NEVER;

    //----- (DEVtrunc)
	<admst:for-each select="$(ddts)">
	{
		std::valarray&lt;FPOLY1&gt; q(OPT::_keep_time_steps);

		for (int ii = 0; ii &lt; OPT::_keep_time_steps; ++ii)
		{
		assert(_states[ii]);
		q[ii].x  = NOT_VALID;
		q[ii].f0 = _states[ii][ _ddt_%(adms/id(.)) ];
		q[ii].f1 = NOT_VALID;
		}

		timestep = std::min(tr_review_trunc_error(&amp;q[0]), timestep);
		trace3(&quot;$(DEV_NAME)::tr_review _ddt_%(adms/id(.))&quot;, timestep,
		   _states[0][ _ddt_%(adms/id(.)) ],
		   _states[0][ _ddt_%(adms/id(.)) + 1 ]
	  );
	}
	</admst:for-each>
    //-----

	_time_by._error_estimate = tr_review_check_and_convert(timestep);
	_time_by._event = NEVER;

  <admst:if test="[$has_sckt!=0]">
	assert(subckt());
	_time_by.min( subckt()-&gt;tr_review() );
  </admst:if>

	return  _time_by;
}
// -------------------------------------------------------------------------- //
<admst:template match="model">
	const $(COMMON_NAME)* c = prechecked_cast&lt;const $(COMMON_NAME)*&gt;(this);
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);
</admst:template>
// -------------------------------------------------------------------------- //
<admst:template match="commonmodel">
	const $(COMMON_NAME)* c = prechecked_cast&lt;const $(COMMON_NAME)*&gt;(common());
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);
</admst:template>
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_accept()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()->tr_accept();
  </admst:if>
  <admst:apply-templates select="." match="commonmodel"/>

  trace0(&quot;$(DEV_NAME)::tr_accept&quot;);
  if (_sim-&gt;analysis_is_dcop() || _sim-&gt;analysis_is_ac()) {
    assert(_dt == NOT_VALID || conchk(_dt, _time[0] - _time[1]));
    // int mode;
    // mode = MODEINITSMSIG;
    // huh??    load(mode);
  } else { itested(); }

	// hmmm
	if(_sim-&gt;uic_now()){
		std::fill_n(it0, (int)$_num_stamp_nodes, 0);
#define NOJAC
#ifdef HAVE_KEEPCOEFF
# define _DDT( qq, stateno) (( _states[0][stateno*2] - _states[1][stateno*2] ) * OPT::keepcoeff)
#endif
<admst:for-each select="$(ddts)">
		double keep_%(adms/id(.))=1;
</admst:for-each>
#		include &quot;$(filename)_tr.hxx&quot;
#undef NOJAC
#undef _DDT
	}

	
}
// -------------------------------------------------------------------------- //
double $(DEV_NAME)::tr_probe_num(const std::string&amp; x)const
{
	trace3("$(DEV_NAME)::tr_probe_num", _dt, _time[0], _time[1]);
	assert(_dt == NOT_VALID || conchk(_dt, _time[0] - _time[1]));
<admst:if test="[($_num_states!='0')]">
	static std::string state_names[_num_states];
</admst:if>
	static std::string contr_names[$_num_stamp_nodes];

  <admst:apply-templates select="." match="commonmodel"/>

<admst:if test="[($_num_states!='0')]">
	// init states &amp; node names
	for(int i=0; i&lt;_num_states; i++)
	{
		std::stringstream ss;

		ss &lt;&lt; i;
		state_names[i] = &quot;state&quot; + ss.str();
	}
</admst:if>
	for(int i=0; i&lt;$_num_stamp_nodes; i++)
	{
		std::stringstream ss;
		ss &lt;&lt; i;
		contr_names[i] = &quot;cont&quot; + ss.str();
	}

<admst:if test="[($_num_states!='0')]">
	// states probes
	for (int ii=0; ii&lt;_num_states; ++ii)
		if (Umatch(x, state_names[ii] + ' '))
			return _states[0][ii];
</admst:if>
	// node probes
	for(int i=0; i&lt;$_num_stamp_nodes; i++)
		if( Umatch(x, contr_names[i] + ' ') )
			return it0[i];

	// verilog-variables ('global_...' eq. with 'ask' field)
  <admst:for-each select="variable[input='no' and scope='global_instance']">
	if( Umatch(x, "%(.) ") ) return $(var_prefix)%(.);
  </admst:for-each>
  <admst:for-each select="variable[input='no' and scope='global_model']">
	if( Umatch(x, "%(.) ") ) return m->$(var_prefix)%(.);
  </admst:for-each>

  return ADMS_BASE::tr_probe_num(x);
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::ac_begin()
{
	ADMS_BASE::ac_begin();
	internal_precalc();
	untested();
	// tr_accept();
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::do_ac()
{
  assert(OPT::_keep_time_steps &lt;= 8);
  assert(_dt == NOT_VALID || conchk(_dt, _time[0] - _time[1]));

  std::fill_n(it0, (int)$_num_stamp_nodes, 0);
  std::fill_n(it1, (int)$_num_stamp_nodes, 0);

  std::fill_n(m_entries, %(count(jacobian)), 0);
  std::fill_n(m_entries_old, %(count(jacobian)), 0);

  //pzLoad(_sim-&gt;_jomega);

  ac_eval();

  #if 0
  <admst:apply-templates select="analog/code[datatypename='block']" match="ac_load" required="yes"/>
  #endif
}
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::ac_load()
{
  <admst:if test="[$has_sckt!=0]">
	subckt()-&gt;ac_load();
  </admst:if>
	for(int i=0; i&lt;$_num_stamp_nodes; i++)
		ac_load_source_point(_n[i], COMPLEX(it0[i], it1[i]));

  <admst:for-each select="jacobian">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
	if( m_required[m_%s_%s] )
		ac_load_point(_N(%s), _N(%s), COMPLEX(m_entries[m_%s_%s], m_entries_old[m_%s_%s] * _sim-&gt;_jomega.imag()));
  </admst:for-each>
}

// -------------------------------------------------------------------------- //
// this is quite a hack.
// wouldnt be necessary if branches were subdevices ...
void $(DEV_NAME)::guesstopology ()
{
	const COMMON_COMPONENT* c = prechecked_cast&lt;const COMMON_COMPONENT*&gt;(common());
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);

	for(int i=0; i&lt;%(count(jacobian)); i++)
        	m_required[i] = false;

	#include &quot;$(filename)_top.hxx&quot;

  for(int i=0; i&lt;%(count(jacobian)); i++){
    trace1(&quot;guess m&quot;, m_required[i]);
  }
}

// -------------------------------------------------------------------------- //
<!--
#if 0
void $(DEV_NAME)::acLoad()
{
<admst:for-each select="jacobian[static='yes']">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:text format="\n\tif(m_required[m_%s_%s]) m_entries[m_%s_%s] += m_ac_static[m_%s_%s]; // static"/>
</admst:for-each>

<admst:text format="\n"/>
<admst:for-each select="jacobian[dynamic='yes']">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:text format="\n\tif(m_required[m_%s_%s]) m_entries_old[m_%s_%s] += _sim->_jomega.imag() * m_ac_dynamic[m_%s_%s]; // dynamic"/>
</admst:for-each>
}

// -------------------------------------------------------------------------- //
void $(DEV_NAME)::pzLoad(COMPLEX s)
{
<admst:for-each select="jacobian[static='yes']">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:text format="\n\tif(m_required[m_%s_%s]) m_entries[m_%s_%s] += m_ac_static[m_%s_%s] * s.real(); // static"/>
</admst:for-each>

<admst:text format="\n"/>
<admst:for-each select="jacobian[dynamic='yes']">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:text format="\n\tif(m_required[m_%s_%s]) m_entries[m_%s_%s] += m_ac_dynamic[m_%s_%s] * s.real(); //dynamic"/>
</admst:for-each>

<admst:text format="\n"/>
<admst:for-each select="jacobian[dynamic='yes']">
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:value-of select="column/name"/>
  <admst:value-of select="row/name"/>
  <admst:text format="\n\tif(m_required[m_%s_%s]) m_entries_old[m_%s_%s] += m_ac_dynamic[m_%s_%s] * s.imag(); //dynamic"/>
</admst:for-each>
}
#endif
<-->

// -------------------------------------------------------------------------- //
bool $(DEV_NAME)::conv_check()const
{
<admst:if test="[($_num_states!='0')]">
	for (unsigned ii = 0; ii &lt; _num_states; ++ii) {
		if (! conchk(_states[0][2*ii], _states_q1[ii]))
			return false;
	}
</admst:if>

	for(unsigned i=0; converged() &amp;&amp; i&lt;$_num_stamp_nodes; i++) {
		bool conv = conchk(it0[i], it1[i]);
		if (!conv) return false;
	}

   // FIXME: use only used entries.
	for(unsigned i=0; converged() &amp;&amp; i&lt;%(count(jacobian)); i++) {
		bool conv = conchk(m_entries[i], m_entries_old[i]);
		if (!conv) {
			trace4(&quot;conv check: not converged&quot;, i, m_entries[i], m_entries_old[i], _sim-&gt;iteration_number() );
			return false;
		}
	}
	return true;
}

// -------------------------------------------------------------------------- //
// write G and C into matrix and matrix_old
void $(DEV_NAME)::ac_eval()
{
	trace0(&quot;$(DEV_NAME)::ac_eval&quot;);
	std::fill_n(it0, (int)$_num_stamp_nodes, 0);
	std::fill_n(m_entries, %(count(jacobian)), 0);
	std::fill_n(m_entries_old, %(count(jacobian)), 0);

	const $(COMMON_NAME)* c = prechecked_cast&lt;const $(COMMON_NAME)*&gt;(common());
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);

	double CKTag0 = 0.; USE(CKTag0); // BUG (necessary?)
	double ddtone = 1.; USE(ddtone); // not necessary in the end?

<admst:for-each select="$(ddts)">
	double keep_%(adms/id(.))=1;
</admst:for-each>
	{
	// fixme: remerge _tr.hxx?
		#include &quot;$(filename)_ac.hxx&quot;
	}
}
// -------------------------------------------------------------------------- //
<admst:if test="[$(_num_states)!=0]">
// ripped from storage...
double $(DEV_NAME)::tr_c_to_g(double c, double g)const
{
//if(_sim-&gt;uic_now() &amp;&amp; have_ic()){
//		return 1; // ./OPT::shortckt;
//	}

	if (_sim-&gt;analysis_is_static()) {
		assert(_time[0] == 0.);
		return 0.;
	}else if (_sim-&gt;analysis_is_restore()) {itested();
		assert(_time[0] &gt; 0);
		return g;
		// no change, fake
	}else{
		assert(_sim-&gt;analysis_is_tran_dynamic());
		METHOD method;
		if (_time[1] == 0) {
		  method = mEULER; // Bogus current in previous step.  Force Euler.
		}else{
		  method = _method_a;
		}
		g = c / _dt;
		switch (method) {
			case mTRAPGEAR: incomplete();
			case mGEAR:	 g *= 3./2.;	break;
			case mTRAPEULER: incomplete();
			case mEULER: /* g *= 1 */	break;
			case mTRAP:	 g *= 2;	break;
		}
		return g;
	}
}
</admst:if>
// -------------------------------------------------------------------------- //
void $(DEV_NAME)::tr_eval()
{
	trace1("$(DEV_NAME)::tr_eval", _sim->analysis_is_static());
	std::fill_n(it0, $_num_stamp_nodes, 0);
	std::fill_n(m_entries, %(count(jacobian)), 0);

<admst:apply-templates select="analog/code[datatypename='block']/variable" match="generate_c:declare_variable"/>

	const $(COMMON_NAME)* c = prechecked_cast&lt;const $(COMMON_NAME)*&gt;(common());
	USE(c); assert(c);
	const $(MODEL_NAME)* m = prechecked_cast&lt;const $(MODEL_NAME)*&gt;(c-&gt;model());
	USE(m); assert(m);

<admst:if test="[$(_num_states)!=0]">
	double CKTag0 = tr_c_to_g(1., 1.);

	//if(_sim-&gt;uic_now() &amp;&amp; have_ic()){

	//	ddtone = 1./OPT::shortckt;
	//}else{
	//	ddtone = 1.;
	//}
<admst:for-each select="$(ddts)">

	double keep_%(adms/id(.))=1;
	if(uic_now(_ddt_%(adms/id(.)))){
		trace1("$(DEV_NAME)::tr_eval keeping", c->ic_state_%(adms/id(.)));
#ifdef HAVE_KEEPCOEFF
		keep_%(adms/id(.)) = OPT::keepcoeff;
#endif
	} else {
		trace3(&quot; $(DEV_NAME)::tr_eval not keeping &quot;, c-&gt;ic_state_%(adms/id(.)),
		    _sim-&gt;uic_now(), uic_now(_ddt_%(adms/id(.))) );
	}
</admst:for-each>

</admst:if>
	
	{
		#include &quot;$(filename)_tr.hxx&quot;
	}

	// this is what cap does.
	set_converged(conv_check());
}

// -------------------------------------------------------------------------- //
$(COMMON_NAME)::$(COMMON_NAME)(const $(COMMON_NAME)&amp; p) :
    COMMON_COMPONENT(p)
<admst:for-each select="$(ddts)">,
    ic_state_%(adms/id(.))(p.ic_state_%(adms/id(.)))
</admst:for-each>
<admst:for-each select="variable[input='yes' and scope='global_instance']">,
    $(var_prefix)%(.)(p.$(var_prefix)%(.))
</admst:for-each>
{
}
// -------------------------------------------------------------------------- //
$(COMMON_NAME)::$(COMMON_NAME)(int c) :
    COMMON_COMPONENT(c)
<admst:for-each select="$(ddts)">,
    ic_state_%(adms/id(.))(NOT_INPUT)
</admst:for-each>
<admst:for-each select="variable[input='yes' and scope='global_instance']">,
    $(var_prefix)%(.)(NOT_INPUT) /* default?? */
</admst:for-each>
{
}
// -------------------------------------------------------------------------- //
$(COMMON_NAME)::~$(COMMON_NAME)()
{
}
// -------------------------------------------------------------------------- //
bool $(COMMON_NAME)::operator==(const COMMON_COMPONENT&amp; x)const
{
	const $(COMMON_NAME)* p = dynamic_cast&lt;const $(COMMON_NAME)*&gt;(&amp;x);
	return (p
<admst:for-each select="$(ddts)">
	    &amp;&amp; ic_state_%(adms/id(.)) == p-&gt;ic_state_%(adms/id(.))
</admst:for-each>
	    &amp;&amp; COMMON_COMPONENT::operator==(x));
}
// -------------------------------------------------------------------------- //
#if 0
map&lt;string, PARA_BASE $(COMMON_NAME)::*&gt; $(COMMON_NAME)::param_dict
<admst:if test="[$(_num_states)!=0]">
  = boost::assign::map_list_of
	<admst:join select="$(ddts)" separator="">
		(&quot;_ddt_%(adms/id(.))&quot;, (PARA_BASE $(COMMON_NAME)::*) 
		(&amp;$(COMMON_NAME)::_ddt_%(adms/id(.)))) // %(.)
	</admst:join>
</admst:if>
;A
#endif
// -------------------------------------------------------------------------- //
void $(COMMON_NAME)::set_param_by_name(std::string name, std::string value)
{
	USE(value);

<admst:for-each select="$ddt_ic_list[.!='']">
	// %(.), %(adms/id(.))
	<admst:apply-templates select="." match="ic_probe_names"/>

	if( Umatch(name,&quot;ic_$(ddt_p)_$(ddt_n)") ) {
		ic_state_%(adms/id(.)) = value; return;
	}
	if( Umatch(name,&quot;ic_$(ddt_n)_$(ddt_p)") ) {
		ic_state_%(adms/id(.)) = "-("+value+")"; return;
	}

</admst:for-each>


#if 0
	PARA_BASE $(COMMON_NAME)::* x = (param_dict[Name]);

  if(x) {
    PARA_BASE* p = &amp;(this-&gt;*x);
    *p = Value;
		return;
	}
#endif
	throw Exception_No_Match(name);
}
// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //

  </admst:for-each>
</admst:open> <!-- cc -->

<!-- *********************************************************************************** -->
<!-- *********************************** constant text ********************************* -->
<!-- *********************************************************************************** -->

<!-- ----------------------------------------------------------------- -->
<admst:template match="c:myfunc:math_h">
<admst:message test="[/dbg_xml='yes']" format="*c:myfunc:math_h*\n"/>
inline double doh(double x)            {return  log(1+exp(x)); }
inline double d_doh(double x)          {return  exp(x)/(exp(x) + 1); }
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="Printhxxwrapperdefines">
<admst:message test="[/dbg_xml='yes']" format="*Printhxxwrapperdefines*\n"/>
// ******************** load definitions ***************************
#define _eval_static_residual1(p,v)\\
	it0[p ## Node]-=v;
#define _keep_static_residual2(p,n,v)\\
	kept[p ## Node]-=v;\\
        kept[n ## Node]+=v;

        <!-- still in use -->
#define _eval_static_residual2(p,n,v)\\
	it0[$node_prefix ##p ]-=v;\\
	it0[$node_prefix ##n]+=v;


</admst:template>

<!-- *********************************************************************************** -->
<!--******************************** tree traverse - datatypename ********************** -->
<!-- *********************************************************************************** -->

<admst:template match="expression:stringify:noprobe">
  <admst:apply-templates select="tree" match="subexpression:differentiate"/>
  <admst:value-of select="/simulator/tmp"/>
</admst:template>
<admst:template match="subexpression:stringify:noprobe">
  <admst:apply-templates select="." match="subexpression:differentiate"/>
  <admst:value-of select="/simulator/tmp"/>
</admst:template>
<admst:template match="subexpression:differentiate">
  <admst:value-of select="./adms/datatypename"/>
  <admst:apply-templates select="." match="%s"/>
  <admst:if test="/simulator/probe">
  // actually differntiating: %(.) wrt %(/simulator/probe)

    <admst:choose>
      <admst:when test="adms[datatypename='probe']">
        <admst:choose>
          <admst:when test="[.=/simulator/probe]">
            <admst:value-to select="/simulator/ddx" value="1.0"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="adms[datatypename='variable']">
        <admst:value-of select="probe"/>
        <admst:if-not-inside select="/simulator/probe" list="%p">
          <admst:value-to select="/simulator/ddx" value="0.0"/>
        </admst:if-not-inside>
        <admst:value-of select="probe"/>
        <admst:if-inside select="/simulator/probe" list="%p">
          <admst:choose>
            <admst:when test="[insource='yes']">
              <admst:value-of select="/simulator/probe/branch/nnode/name"/>
              <admst:value-of select="/simulator/probe/branch/pnode/name"/>
              <admst:value-of select="/simulator/probe/nature/access"/>
              <admst:value-of select="name"/>
              <admst:value-to select="/simulator/ddx" value="%s_%s%s_%s"/>
            </admst:when>
            <admst:otherwise>
              <admst:value-to select="/simulator/ddx" value="0.0"/>
            </admst:otherwise>
          </admst:choose>
        </admst:if-inside>
      </admst:when>
      <admst:when test="adms[ datatypename='number' or datatypename='variable']">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:when>
    </admst:choose>
  </admst:if>
</admst:template>


<!-- analog//callfunctions -->
<admst:template match="callfunction">
  <admst:choose>
    <admst:when test="function[name='\$strobe']">
      <admst:text format="_strobe("/>
    </admst:when>
    <admst:when test="function[name='\$warning']">
      <admst:text format="_warning("/>
    </admst:when>
    <admst:when test="function[name='\$error']">
      <admst:text format="_error("/>
    </admst:when>
    <admst:when test="function[name='\$finish']">
      <admst:text format="_finish("/>
    </admst:when>
    <admst:when test="function[name='\$stop']">
      <admst:text format="_stop("/>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="/simulator/tmp"/>
      <admst:value-of select="function/name"/>
      <admst:error format="function not supported: %s(%s)\n"/>
    </admst:otherwise>
  </admst:choose>
  <admst:reset select="/simulator/tmp"/>
  <admst:join select="function/arguments" separator=",">
    <admst:apply-templates select="." match="expression:stringify:noprobe"/>%s</admst:join>
  <admst:text format=");\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="mapply_unary">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)mapply_unary*\n"/>
  <admst:choose>
    <admst:when test="[name='plus']">
      <admst:choose>
        <admst:when test="arg1/math[value=0.0]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(+%s)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='minus']">
      <admst:choose>
        <admst:when test="arg1/math[value=0.0]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(-%s)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='not']">
      <admst:choose>
        <admst:when test="arg1/math[value=0.0]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(!%s)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='bw_not']">
      <admst:choose>
        <admst:when test="arg1/math[value=0.0]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(~%s)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="name"/>
      <admst:error format="%s: function not handled\n"/>
    </admst:otherwise>
  </admst:choose>
  <admst:if test="/simulator/probe">
    <admst:choose>
      <admst:when test="/simulator[tmp='0.0']">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:when>
      <admst:otherwise>
        <admst:choose>
          <admst:when test="[name='plus']">
            <admst:choose>
              <admst:when test="/simulator[ddx!='0.0']">
                <admst:value-of select="/simulator/ddx"/>
                <admst:value-to select="/simulator/ddx" value="(+%s)"/>
              </admst:when>
            </admst:choose>
          </admst:when>
          <admst:when test="[name='minus']">
            <admst:choose>
              <admst:when test="/simulator[ddx!='0.0']">
                <admst:value-of select="/simulator/ddx"/>
                <admst:value-to select="/simulator/ddx" value="(-%s)"/>
              </admst:when>
            </admst:choose>
          </admst:when>
          <admst:when test="[name='not']">
            <admst:choose>
              <admst:when test="/simulator[ddx!='0.0']">
                <admst:value-of select="/simulator/ddx"/>
                <admst:value-to select="/simulator/ddx" value="(!%s)"/>
              </admst:when>
            </admst:choose>
          </admst:when>
          <admst:when test="[name='bw_not']">
            <admst:choose>
              <admst:when test="/simulator[ddx!='0.0']">
                <admst:value-of select="/simulator/ddx"/>
                <admst:value-to select="/simulator/ddx" value="(~%s)"/>
              </admst:when>
            </admst:choose>
          </admst:when>
        </admst:choose>
      </admst:otherwise>
    </admst:choose>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="mapply_binary">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)mapply_binary*\n"/>
  <admst:choose>
    <admst:when test="[name='bw_equr']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s^~%s)"/>
    </admst:when>
    <admst:when test="[name='bw_equl']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s~^%s)"/>
    </admst:when>
    <admst:when test="[name='bw_xor']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s^%s)"/>
    </admst:when>
    <admst:when test="[name='bw_or']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s|%s)"/>
    </admst:when>
    <admst:when test="[name='bw_and']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&amp;%s)"/>
    </admst:when>
    <admst:when test="[name='multmod']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s%%%s)"/>
    </admst:when>
    <admst:when test="[name='or']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s||%s)"/>
    </admst:when>
    <admst:when test="[name='and']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&amp;&amp;%s)"/>
    </admst:when>
    <admst:when test="[name='equ']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s==%s)"/>
    </admst:when>
    <admst:when test="[name='notequ']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s!=%s)"/>
    </admst:when>
    <admst:when test="[name='lt']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&lt;%s)"/>
    </admst:when>
    <admst:when test="[name='lt_equ']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&lt;=%s)"/>
    </admst:when>
    <admst:when test="[name='gt']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&gt;%s)"/>
    </admst:when>
    <admst:when test="[name='gt_equ']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&gt;=%s)"/>
    </admst:when>
    <admst:when test="[name='shiftr']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&gt;&gt;%s)"/>
    </admst:when>
    <admst:when test="[name='shiftl']">
      <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
      <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
      <admst:value-to select="/simulator/tmp" value="(%s&lt;&lt;%s)"/>
    </admst:when>
    <admst:when test="[name='addp']">
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)and(arg2/math/value=0.0)]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="arg1/math[value=0.0]">
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(+%s)"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
            <admst:variable name="dx" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="arg2/math[value=0.0]">
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:variable name="x" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
          </admst:if>
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:variable name="y" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
          </admst:if>
          <admst:value-to select="/simulator/tmp" value="($x+$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='addm']">
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)and(arg2/math/value=0.0)]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="arg1/math[value=0.0]">
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="(-%s)"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
            <admst:variable name="dx" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="arg2/math[value=0.0]">
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:value-to select="/simulator/tmp" value="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:variable name="x" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
          </admst:if>
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:variable name="y" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
          </admst:if>
          <admst:value-to select="/simulator/tmp" value="($x-$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='multtime']">
      <admst:variable name="x" select="0.0"/>
      <admst:variable name="y" select="0.0"/>
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)or(arg2/math/value=0.0)]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="[(arg1/math/value=1.0)and(arg2/math/value=1.0)]">
          <admst:value-to select="/simulator/tmp" value="1.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:variable name="x" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
          </admst:if>
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
          </admst:if>
          <admst:variable name="y" select="%s"/>
          <admst:value-to select="/simulator/tmp" value="($x*$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='multdiv']">
      <admst:variable name="x" select="0.0"/>
      <admst:variable name="y" select="0.0"/>
      <admst:choose>
        <admst:when test="arg1/math[value=0.0]">
          <admst:value-to select="/simulator/tmp" value="0.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:when test="[(arg1/math/value=1.0)and(arg2/math/value=1.0)]">
          <admst:value-to select="/simulator/tmp" value="1.0"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dx" select="0.0"/>
            <admst:variable name="dy" select="0.0"/>
          </admst:if>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
          <admst:variable name="x" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dx" select="%s"/>
          </admst:if>
          <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
          <admst:variable name="y" select="%s"/>
          <admst:if test="/simulator/probe">
            <admst:value-of select="/simulator/ddx"/>
            <admst:variable name="dy" select="%s"/>
          </admst:if>
          <admst:value-to select="/simulator/tmp" value="($x/$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="name"/>
      <admst:error format="%s: function not handled\n"/>
    </admst:otherwise>
  </admst:choose>

  <admst:if test="/simulator/probe">
    <admst:choose>
      <admst:when test="[name='addp']">
        <admst:choose>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:value-to select="/simulator/ddx" value="(+$dy)"/>
          </admst:when>
          <admst:when test="[$dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="$dx"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-to select="/simulator/ddx" value="($dx+$dy)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='addm']">
        <admst:choose>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:value-to select="/simulator/ddx" value="(-$dy)"/>
          </admst:when>
          <admst:when test="[$dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="$dx"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-to select="/simulator/ddx" value="($dx-$dy)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multtime']">
        <admst:choose>
          <admst:when test="[$x='0.0' and $y='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0' and $dy='1.0']">
            <admst:value-to select="/simulator/ddx" value="($x)"/>
          </admst:when>
          <admst:when test="[$dx='1.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="($y)"/>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:value-to select="/simulator/ddx" value="($x*$dy)"/>
          </admst:when>
          <admst:when test="[$dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="$dx*$y"/>
          </admst:when>
          <admst:when test="[$dx='1.0' and $dy='1.0']">
            <admst:value-to select="/simulator/ddx" value="($x+$y)"/>
          </admst:when>
          <admst:when test="[$dx='1.0']">
            <admst:value-to select="/simulator/ddx" value="($y+($dy*$x))"/>
          </admst:when>
          <admst:when test="[$dy='1.0']">
            <admst:value-to select="/simulator/ddx" value="($dx*$y)+$x"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:value-to select="/simulator/ddx" value="$dy"/>
          </admst:when>
          <admst:when test="[$y='1.0']">
            <admst:value-to select="/simulator/ddx" value="$dx"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-to select="/simulator/ddx" value="(($dx*$y)+($x*$dy))"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multdiv']">
        <admst:choose>
          <admst:when test="[$x='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:choose>
              <admst:when test="[$dy='1.0']">
                <admst:value-to select="/simulator/ddx" value="(-1/($y*$y))"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-to select="/simulator/ddx" value="(-$dy/($y*$y))"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:choose>
              <admst:when test="[$dy='1.0']">
                <admst:value-to select="/simulator/ddx" value="(-$x/($y*$y))"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-to select="/simulator/ddx" value="(-($x*$dy)/($y*$y))"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$dx='1.0']">
            <admst:choose>
              <admst:when test="[$dy='0.0']">
                <admst:value-to select="/simulator/ddx" value="(1./$y)"/>
              </admst:when>
              <admst:when test="[$dy='1.0']">
                <admst:value-to select="/simulator/ddx" value="(($y-$x)/($y*$y))"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-to select="/simulator/ddx" value="(($y-($x*$dy))/($y*$y))"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:otherwise>
            <admst:choose>
              <admst:when test="[$y='1.0']">
                <admst:value-to select="/simulator/ddx" value="$dx"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:value-to select="/simulator/ddx" value="$dx/$y"/>
              </admst:when>
              <admst:when test="[$dy='1.0']">
                <admst:value-to select="/simulator/ddx" value="(($dx*$y)-$x)/($y*$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-to select="/simulator/ddx" value="($dx*$y-$x*$dy)/($y*$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:otherwise>
        <admst:value-to select="/simulator/ddx" value=""/>
      </admst:otherwise>
    </admst:choose>
  </admst:if>

</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="mapply_ternary">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)mapply_ternary*\n"/>
  <admst:apply-templates select="arg3" match="subexpression:stringify:noprobe"/>
  <admst:value-of select="/simulator/ddx"/>
  <admst:variable name="dz" select="%s"/>
  <admst:apply-templates select="arg2" match="subexpression:stringify:noprobe"/>
  <admst:value-of select="/simulator/ddx"/>
  <admst:variable name="dy" select="%s"/>
  <admst:apply-templates select="arg1" match="subexpression:stringify:noprobe"/>
  <admst:variable name="x" select="%s"/>
  <admst:value-to select="/simulator/tmp" value="($x?%s:%s)"/>
  <admst:if test="/simulator/probe">
    <admst:value-to select="/simulator/ddx" value="($x?$dy:$dz)"/>
  </admst:if>
</admst:template>

<!-- *********************************************************************************** -->
<!-- ******************************* tree traverse ************************************* -->
<!-- *********************************************************************************** -->

<!-- ----------------------------------------------------------------- -->
<!-- analog/![initializeModel|initializeInstance|initial_model|initial_instance|initial_step|noise] -->
<!-- from ngspiceMODULEload.c.xml -->
<admst:template match="analog:evaluate">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:evaluate*\n"/>
  // evaluate template
  <admst:if test="code">
    <admst:assert test="code/adms[datatypename='block']" format="expecting datatypename=block\n"/>
    <admst:for-each select="code/item">
      <admst:if test="adms[datatypename='block']">
        <admst:if test="[name!='initial_model' and name!='initial_instance'
           and name!='initializeModel' and name!='initializeInstance'
           and name!='initial_step']">
          <admst:apply-templates select="." match="block:local:declaration"/>
        </admst:if>
      </admst:if>
      <admst:if test="adms[datatypename!='block']">
        <admst:apply-templates select="." match="block:local:declaration"/>
      </admst:if>
    </admst:for-each>
    <admst:apply-templates select="code" match="variable:declaration"/>
    <admst:for-each select="code/item">
      <admst:choose>
        <admst:when test="adms[datatypename!='block']">
          <admst:value-of select="./adms/datatypename"/>
          <admst:apply-templates select="." match="ngspiseMODULEload%s"/>
        </admst:when>
        <admst:otherwise>
          <admst:if test="[name!='initial_model' and name!='initial_instance'
            and name!='initializeModel' and name!='initializeInstance'
            and name!='initial_step']">
            <admst:apply-templates select="." match="block"/>
          </admst:if>
        </admst:otherwise>
      </admst:choose>
    </admst:for-each>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- -------- expression//function: mapping verilog-name == C-name of function ------- -->
<admst:template match="function:getname">
  <admst:message test="[/dbg_xml='yes']" format="*MODULE function:getname: %(.)*\n"/>
  <admst:choose>
    <admst:when test="[name='abs']"><admst:return name="function:getname" value="abs"/></admst:when>
    <admst:when test="[name='\$shrinkl']"><admst:return name="function:getname" value="shrinkl"/></admst:when>
    <admst:when test="[name='\$shrinka']"><admst:return name="function:getname" value="shrinka"/></admst:when>
    <admst:when test="[name='log']"><admst:return name="function:getname" value="log10"/></admst:when>
    <admst:when test="[name='ln']"><admst:return name="function:getname" value="logE"/></admst:when>
    <admst:when test="[name='limexp']"><admst:return name="function:getname" value="limexp"/></admst:when>
    <admst:when test="[name='\$limexp']"><admst:return name="function:getname" value="limexp"/></admst:when>
    <admst:when test="[name='\$model']"><admst:return name="function:getname" value="_modelname"/></admst:when>
    <admst:when test="[name='\$instance']">
      <admst:return name="function:getname" value="_instancename"/>
    </admst:when>
    <admst:when test="[name='\$temperature']">
      <admst:return name="function:getname" value="_ambient_temp"/>
    </admst:when>
    <admst:when test="[name='\$nominal_temperature']">
      <admst:return name="function:getname" value="_circuit_tnom"/>
    </admst:when>
    <admst:otherwise>
      <admst:message test="[/dbg_xml='yes']" format="otherwise (not impl?) %(name)\n"/>
      <admst:value-of select="name"/>
      <admst:return name="function:getname" value="%s"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<admst:template match="function:getname_push">
  <admst:message format="*function:getname_push %(.)*\n"/>
  <!--  <admst:value-of select="(%(functiongetname)(.)/[name='functiongetname'])"/> -->
  <admst:apply-templates select="." match="function:getname">
	  <admst:value-of select="returned('function:getname')/value"/>
  </admst:apply-templates>
</admst:template>

<!-- -------------------------------------------------------------- -->
<admst:template match="function">
  <admst:message format="* check: old function template?\n"/>
  <admst:choose>
    <admst:when test="[name='ddt']">
      <admst:for-each select="arguments[position(.)=1]">
        <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
        <admst:value-to select="/simulator/tmp" value="/*tmp*/ %s"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="[name='\$given']">
      <admst:for-each select="arguments[position(.)=1]">
        <admst:if test="adms[datatypename!='variable']">
          <admst:error format="$given: argument is not a variable\n"/>
        </admst:if>
        <admst:if test="[input='no']">
          <admst:value-of select="name"/>
          <admst:error format="$given(%s): argument is not a parameter\n"/>
        </admst:if>
        <admst:choose>
          <admst:when test="[parametertype='model']">
            <admst:value-to select="/simulator/tmp" value="m->$(var_prefix)%(name).has_good_value()"/>
          </admst:when>
          <admst:when test="[parametertype='instance']">
            <admst:value-to select="/simulator/tmp" value="/*instance*/$(var_prefix)%(name).has_good_value()"/>
          </admst:when>
          <admst:otherwise>
            <admst:error format="$given(%s): should not be reached\n"/>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$model']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="%s"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$instance']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="%s"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$temperature']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="%s"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$nominal_temperature']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="%s"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$vt']">
      <admst:choose>
       <admst:when test="arguments">
        <admst:choose>
         <admst:when test="arguments[count(.)=1]">
          <admst:apply-templates select="." match="function:assert:onearg"/>
          <admst:for-each select="arguments[position(.)=1]">
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
            <admst:value-to select="/simulator/tmp" value="_vt(%s)"/>
          </admst:for-each>
         </admst:when>
         <admst:otherwise>
           <admst:error format="$vt(...): too many args"/>
         </admst:otherwise>
        </admst:choose>
       </admst:when>
       <admst:otherwise>
         <admst:apply-templates select="." match="function:assert:noarg"/>
         <admst:value-to select="/simulator/tmp" value="_vt_nom"/>
       </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='\$scale']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:value-to select="/simulator/tmp" value="_scale"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$abstime']">
      <admst:apply-templates select="." match="function:assert:noarg"/>
      <admst:value-to select="/simulator/tmp" value="_abstime"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='\$options']">
      <admst:for-each select="arguments[position(.)=1]">
        <admst:if test="adms[datatypename!='string']">
          <admst:error format="$given: argument is not a string\n"/>
        </admst:if>
        <admst:choose>
          <admst:when test="[value='OPTm_hier']">
            <admst:value-to select="/simulator/tmp" value="_circuit_m_hier"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-of select="value"/>
            <admst:fatal format="$options(%s): bad argument []\n"/>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='ddx' or name='\$derivate']">
      <admst:for-each select="arguments">
        <admst:if test="[position(.)=2]">
          <admst:if test="adms[datatypename!='probe']">
            <admst:value-of select="../name"/>
            <admst:error format="%s: second argument is not a probe\n"/>
          </admst:if>
          <admst:value-of select="branch/nnode/name"/>
          <admst:value-of select="branch/pnode/name"/>
          <admst:value-of select="nature/access"/>
        </admst:if>
      </admst:for-each>
      <admst:for-each select="arguments">
        <admst:if test="[position(.)=1]">
          <admst:if test="adms[datatypename!='variable']">
            <admst:value-of select="../name"/>
            <admst:error format="%s: first argument is not a variable\n"/>
          </admst:if>
          <admst:value-of select="name"/>
        </admst:if>
      </admst:for-each>
      <admst:value-to select="/simulator/tmp" value="%s_%s%s_%s"/>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='floor']">
      <admst:apply-templates select="." match="function:assert:onearg"/>
      <admst:for-each select="arguments[position(.)=1]">
        <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
        <admst:value-to select="/simulator/tmp" value="floor(%s)"/>
      </admst:for-each>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='ceil']">
      <admst:apply-templates select="." match="function:assert:onearg"/>
      <admst:for-each select="arguments[position(.)=1]">
        <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
        <admst:value-to select="/simulator/tmp" value="ceil(%s)"/>
      </admst:for-each>
      <admst:if test="/simulator/probe">
        <admst:value-to select="/simulator/ddx" value="0.0"/>
      </admst:if>
    </admst:when>
    <admst:when test="[name='pow' or name='hypot' or name='min' or name='max']">
      <admst:value-of select="index(./subexpression/expression/function,.)"/>
      <admst:variable name="index" select="%s"/>
      <admst:if test="/simulator/probe">
        <admst:for-each select="arguments">
          <admst:choose>
            <admst:when test="[position(.)=1]">
              <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
              <admst:variable name="x" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dx" select="%s"/>
            </admst:when>
            <admst:when test="[position(.)=2]">
              <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
              <admst:variable name="y" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dy" select="%s"/>
            </admst:when>
            <admst:otherwise>
              <admst:count select="../arguments"/>
              <admst:value-of select="../name"/>
              <admst:error format="%s(...): two arguments expected - %s found(s) \n"/>
            </admst:otherwise>
          </admst:choose>
        </admst:for-each>
        <admst:choose>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:value-of select="name"/>
            <admst:value-to select="/simulator/ddx" value="(__dFy_%s_$index*$dy)"/>
          </admst:when>
          <admst:when test="[$dy='0.0']">
            <admst:value-of select="name"/>
            <admst:value-to select="/simulator/ddx" value="(__dFx_%s_$index*$dx)"/>
          </admst:when>
          <admst:otherwise>
            <admst:apply-templates select="." match="function:getname_push"/>
            <admst:apply-templates select="." match="function:getname_push"/>
            <admst:value-to select="/simulator/ddx" value="(__dFx_%s_$index*$dx+__dFy_%s_$index*$dy)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:if>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="__%s_$index"/>
    </admst:when>
    <admst:when test="[name='div']">
      <admst:value-of select="index(./subexpression/expression/function,.)"/>
      <admst:variable name="index" select="%s"/>
      <admst:if test="/simulator/probe">
        <admst:for-each select="arguments">
          <admst:choose>
            <admst:when test="[position(.)=1]">
              <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
              <admst:variable name="x" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dx" select="%s"/>
            </admst:when>
            <admst:when test="[position(.)=2]">
              <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
              <admst:variable name="y" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dy" select="%s"/>
            </admst:when>
            <admst:otherwise>
              <admst:count select="../arguments"/>
              <admst:value-of select="../name"/>
              <admst:error format="%s(...): two arguments expected - %s found(s) \n"/>
            </admst:otherwise>
          </admst:choose>
        </admst:for-each>
        <admst:choose>
          <admst:when test="[$dx='0.0' and $dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:when test="[$dx='0.0']">
            <admst:value-to select="/simulator/ddx" value="(__dFy_%(name)_$index*$dy)"/>
          </admst:when>
          <admst:when test="[$dy='0.0']">
            <admst:value-to select="/simulator/ddx" value="(__dFx_%(name)_$index*$dx)"/>
          </admst:when>
          <admst:otherwise>
            <admst:value-to select="/simulator/ddx" value="(__dFx_%(name)_$index*$dx+__dFy_%(name)_$index*$dy)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:if>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="__%s_$index"/>
    </admst:when>
    <admst:when test="[class='builtin']">
      <admst:value-of select="index(./subexpression/expression/function,.)"/>
      <admst:variable name="index" select="%s"/>
      <admst:if test="/simulator/probe">
        <admst:for-each select="arguments">
          <admst:choose>
            <admst:when test="[position(.)=1]">
              <admst:apply-templates select="." match="stringifry_expression"/>
              <admst:variable name="x" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dx" select="%s"/>
            </admst:when>
            <admst:otherwise>
              <admst:count select="../arguments"/>
              <admst:value-of select="../name"/>
              <admst:error format="%s(...): one argument expected - %s found(s) \n"/>
            </admst:otherwise>
          </admst:choose>
        </admst:for-each>
        <admst:choose>
          <admst:when test="[$dx='0.0']">
            <admst:value-to select="/simulator/ddx" value="0.0"/>
          </admst:when>
          <admst:otherwise>
            <admst:apply-templates select="." match="function:getname_push"/>
            <admst:value-to select="/simulator/ddx" value="$dx*__d_%s_$index"/>
          </admst:otherwise>
        </admst:choose>
      </admst:if>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:value-to select="/simulator/tmp" value="__%s_$index"/>
    </admst:when>
    <admst:otherwise>
      <admst:choose>
        <admst:when test="[name='\$simparam']">
           <admst:apply-templates select="." match="function:simparam"/>
        </admst:when>
        <admst:when test="[name='analysis']">
           <admst:apply-templates select="." match="function:analysis"/>
        </admst:when>
        <admst:otherwise>
          <admst:value-of select="name"/>
          <admst:variable name="function" select="%s"/>
          <admst:variable name="args" select=""/>
          <admst:for-each select="arguments">
            <admst:value-of select="position(.)"/>
            <admst:variable name="index" select="%s"/>
            <admst:if test="[not($args='')]">
              <admst:variable name="args" select="$args,"/>
            </admst:if>
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
            <admst:variable name="arg$index" select="%s"/>
            <admst:variable name="args" select="$args$(arg$index)"/>
          </admst:for-each>
          <admst:value-to select="/simulator/tmp" value="$function($args)"/>
          <admst:if test="/simulator/probe">
            <admst:variable name="dargs" select="$args"/>
            <admst:for-each select="arguments">
              <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
              <admst:variable name="x" select="%s"/>
              <admst:value-of select="/simulator/ddx"/>
              <admst:variable name="dargs" select="$dargs,%s"/>
            </admst:for-each>
            <admst:value-to select="/simulator/ddx" value="d_$function($dargs)"/>
          </admst:if>
        </admst:otherwise>
      </admst:choose>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="function:assert:noarg">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)function:assert:noarg*\n"/>
  <admst:if test="[not(nilled(arguments))]">
    <admst:value-of select="name"/>
    <admst:error format="%s: should not have arguments\n"/>
  </admst:if>
</admst:template>
<admst:template match="function:assert:onearg">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)function:assert:onearg*\n"/>
  <admst:if test="arguments[not(count(.)=1)]">
    <admst:value-of select="name"/>
    <admst:error format="%s: should have one argument exactly\n"/>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="function:analysis">
<admst:message test="[/dbg_xml='yes']" format="*(ngspvers)function:analysis*\n"/>
  <admst:value-of select="name"/>
  <admst:variable name="function" select="%s"/>
  <admst:variable name="args" select=""/>
  <admst:for-each select="arguments">
    <admst:value-of select="position(.)"/>
    <admst:variable name="index" select="%s"/>
    <admst:if test="[not($args='')]">
      <admst:variable name="args" select="$args,"/>
    </admst:if>
    <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
    <admst:variable name="arg$index" select="%s"/>
    <admst:variable name="args" select="$args$(arg$index)"/>
  </admst:for-each>
  <admst:choose>
    <admst:when test="[$arg1='&quot;noise&quot;']">
      <admst:value-to select="/simulator/tmp" value="0.0"/>
      <admst:error format="$function($args): replaced by 0.0\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:error format="$function($args) -- not implemented in ngspice interface\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="function:simparam">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)function:simparam*\n"/>
  <admst:value-of select="name"/>
  <admst:variable name="function" select="%s"/>
  <admst:variable name="args" select=""/>
  <admst:for-each select="arguments">
    <admst:value-of select="position(.)"/>
    <admst:variable name="index" select="%s"/>
    <admst:if test="[not($args='')]">
      <admst:variable name="args" select="$args,"/>
    </admst:if>
    <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
    <admst:variable name="arg$index" select="%s"/>
    <admst:variable name="args" select="$args$(arg$index)"/>
  </admst:for-each>
  <admst:choose>
    <admst:when test="[$arg1='&quot;gdev&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_gdev"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;gmin&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_gmin"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;imax&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_imax"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;imelt&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_imelt"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;iteration&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_iteration"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;scale&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_scale"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;shrink&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_shrink"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;simulatorSubversion&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_simulatorSubversion"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;simulatorVersion&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_simulatorVersion"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;sourceScaleFactor&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_sourceScaleFactor"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;tnom&quot;']">
      <admst:value-to select="/simulator/tmp" value="_circuit_tnom"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;checkjcap&quot;']">
      <admst:value-to select="/simulator/tmp" value="1.0"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;maxmosl&quot;']">
      <admst:value-to select="/simulator/tmp" value="1.0"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;maxmosw&quot;']">
      <admst:value-to select="/simulator/tmp" value="1.0"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;minmosl&quot;']">
      <admst:value-to select="/simulator/tmp" value="1.0e-12"/>
    </admst:when>
    <admst:when test="[$arg1='&quot;minmosw&quot;']">
      <admst:value-to select="/simulator/tmp" value="1.0e-12"/>
    </admst:when>
    <admst:otherwise>
      <admst:error format="$function($args) -- not implemented in ngspice interface\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="variable:declaration">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)variable:declaration*\n"/>
  <admst:for-each select="module/evaluation/variable">
    <admst:assert test="adms[datatypename='variable']" format="expecting datatypename=variable\n"/>
    <admst:if test="[scope='local']">
      <admst:if test="[static='no' and dynamic='yes']">#if defined(_DYNAMIC)\n</admst:if>
      <admst:value-of select="name"/>
      <admst:if test="[type='integer']">int %s;\n</admst:if>
      <admst:if test="[type='real']">double %s=0.0/0.0;\n</admst:if>
      <admst:if test="[type='string']">char* %s;\n</admst:if>
      <admst:if test="[insource='yes']">
        <admst:if test="probe">
          <admst:text format="#if defined(_DERIVATE) // BAR\n"/>
          <admst:for-each select="probe">
            <admst:value-of select="branch/nnode/name"/>
            <admst:value-of select="branch/pnode/name"/>
            <admst:value-of select="nature/access"/>
            <admst:value-of select="../name"/>
            <admst:text format="double %s_%s%s_%s=0.0;\n"/>
          </admst:for-each>
          <admst:text format="#endif /*_DERIVATE*/\n"/>
        </admst:if>
      </admst:if>
      <admst:if test="[static='no' and dynamic='yes']">#endif /*_DYNAMIC*/\n</admst:if>
    </admst:if>
    <admst:if test="[scope!='local']">
      <admst:if test="[insource='yes']">
        <admst:if test="probe">
          <!-- admst:text format="#if defined(_DERIVATE) // BLI\n"/ -->
          <admst:if test="[$_DERIVATE='yes']">
            <admst:for-each select="probe">
              <admst:value-of select="branch/nnode/name"/>
              <admst:value-of select="branch/pnode/name"/>
              <admst:value-of select="nature/access"/>
              <admst:value-of select="../name"/>
              <admst:text format="double %s_%s%s_%s=0.0;\n"/>
            </admst:for-each>
          <!--admst:text format="#endif /*_DERIVATE*/\n"/-->
          </admst:if>
        </admst:if>
      </admst:if>
    </admst:if>
  </admst:for-each>
  <admst:reset select="module/evaluation/variable"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- save all variables used for local declaration -->
<admst:template match="block:local:declaration">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)block:local:declaration*\n"/>
  <admst:text format="// block:local:declaration\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='assignment']">
      <admst:push into="module/evaluation/variable" select="lhs" onduplicate="ignore"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:for-each select="item">
        <admst:apply-templates select="." match="block:local:declaration" required="yes"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="then" match="block:local:declaration" required="yes"/>
      <admst:apply-templates select="else" match="block:local:declaration" required="yes"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="block:local:declaration" required="yes"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']">
    </admst:when>
    <admst:when test="adms[datatypename='nilled']">
    </admst:when>
    <admst:when test="adms[datatypename='callfunction']">
    </admst:when>
    <admst:when test="adms[datatypename='case']">
      <admst:error format="case statement: please implement me! (local declaration)\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='blockvariable']"/>
    <admst:otherwise>
      <admst:value-of select="admst(.)"/>
      <admst:value-of select="adms/datatypename"/>
      <admst:error format="'datatypename=%s': should not be reached %s\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//blockvariable -->
<admst:template match="blockvariable">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)blockvariable*\n"/>
  <admst:for-each select="variable">
    <admst:if test="[type='integer']">int $(var_prefix)%(name);\n</admst:if>
    <admst:if test="[type='real']">double $(var_prefix)%(name);\n</admst:if>
    <admst:if test="[type='string']">char* $(var_prefix)%(name);\n</admst:if>
    <admst:text test="[insource='yes']" select="probe" format="double %(../name)_%(nature/access)%(branch/pnode/name)_%(branch/nnode/name);\n"/>
  </admst:for-each>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//block -->
<admst:template match="block">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)block*\n"/>
  <admst:assert test="[name!='/']" format="expecting subblock\n"/>
  <admst:text format="{\n"/>
  <admst:for-each select="item">
    <admst:value-of select="./adms/datatypename"/>
    <admst:apply-templates select="." match="%s" required="yes"/>
  </admst:for-each>
  <admst:text format="}\n"/>
</admst:template>

<!-- from here is not guesstopology, but model_initialize, at least next 2 functions  -->
<!-- ----------------------------------------------------------------- -->
<!-- analog/[initializeModel|initializeInstance|initial_model|initial_instance|initial_step|noise] -->
<admst:template match="block:initial">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)block:initial*\n"/>
  <admst:assert test="adms[datatypename='block']" format="expecting datatypename=block\n"/>
  <admst:apply-templates select="." match="block:local:declaration"/>
  <admst:apply-templates select="." match="variable:declaration"/>
  <admst:apply-templates select="." match="block" required="yes"/>
</admst:template>
<admst:template match="analog:initial_instance">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:initial_instance*\n"/>
  <admst:if test="code">
    <admst:if test="code/adms[datatypename='block']">
      <admst:for-each select="code/item">
        <admst:if test="adms[datatypename='block']">
         <admst:apply-templates select="[name='initial_instance' or name='initializeInstance']"
           match="block:initial"/>
        </admst:if>
      </admst:for-each>
    </admst:if>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="analog:initial_model">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:initial_model*\n"/>
  <admst:if test="code">
    <admst:if test="code/adms[datatypename='block']">
      <admst:for-each select="code/item">
        <admst:if test="adms[datatypename='block']">
        <admst:apply-templates select="[name='initial_model' or name='initializeModel']"
          match="block:initial"/>
        </admst:if>
      </admst:for-each>
    </admst:if>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="analog:initial_step">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:initial_step*\n"/>
  <admst:if test="code">
    <admst:if test="code/adms[datatypename='block']">
      <admst:for-each select="code/item">
        <admst:if test="adms[datatypename='block']">
          <admst:apply-templates select="[name='initial_step']" match="block:initial"/>
        </admst:if>
      </admst:for-each>
    </admst:if>
  </admst:if>
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:initial_step<-*\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="analog:noise">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)analog:noise*\n"/>
  <admst:if test="code">
    <admst:if test="code/adms[datatypename='block']">
      <admst:for-each select="code/item">
        <admst:if test="adms[datatypename='block']">
          <admst:apply-templates select="[name='noise']" match="block:initial"/>
        </admst:if>
      </admst:for-each>
    </admst:if>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//function: local assignment handling -->
<admst:template match="function:assignment">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)function:assignment %(.)*\n"/>
  <admst:for-each select="function[class='builtin']">
    <admst:choose>
      <admst:when test="arguments[count(.)=1]">
        <admst:value-of select="position(.)-1"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:text format="_%s(__%s_%s,"/>
        <admst:join select="arguments" separator=",">
          <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
          <admst:text format="(%s)"/>
        </admst:join>
        <admst:text format=")\n"/>
        <admst:value-of select="position(.)-1"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:text format="EXIT_IF_ISNAN(__%s_%s)\n"/>
      </admst:when>
      <admst:when test="arguments[count(.)=2]">
        <admst:value-of select="position(.)-1"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:text format="_%s"/>
        <admst:text test="[name='div']" format="0"/>
        <admst:text format="(__%s_%s,"/>
        <admst:join select="arguments" separator=",">
          <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
          <admst:text format="%s"/>
        </admst:join>
        <admst:text format=")\n"/>
        <admst:value-of select="position(.)-1"/>
        <admst:apply-templates select="." match="function:getname_push"/>
        <admst:text format="EXIT_IF_ISNAN(__%s_%s)\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:value-of select="name"/>
        <admst:error format="%s: function not handled\n"/>
      </admst:otherwise>
    </admst:choose>
  </admst:for-each>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="ddx:function:computation:reference">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)ddx:function:computation\n"/>
  <admst:if test="[datatypename='contribution' or lhs/insource='yes']">
	  <!-- <admst:if test="lhs[insource='yes']"> -->
    <admst:if test="rhs[hasVoltageDependentFunction='yes']">
      <admst:text format="#if defined(_DERIVATE)\n"/>
      <admst:for-each select="rhs/function">
        <admst:if test="arguments[count(.)=1]">
          <admst:for-each select="arguments[position(.)=1]">
            <admst:if test="math[dependency!='constant']">
              <admst:value-of select="../position(.)-1"/>
              <admst:apply-templates select=".." match="function:getname_push"/>
              <admst:text format="double __d_%s_%s=0.0;\n"/>
            </admst:if>
          </admst:for-each>
        </admst:if>
        <admst:if test="arguments[count(.)=2]">
          <admst:for-each select="arguments">
            <admst:if test="[position(.)=1]">
              <admst:if test="[(../name='div') or (math/dependency!='constant')]">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select="." match="function:getname_push"/>
                <admst:text format="double __dFx_%s_%s=0.0;\n"/>
              </admst:if>
            </admst:if>
            <admst:if test="[position(.)=2]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="double __dFy_%s_%s=0.0;\n"/>
              </admst:if>
            </admst:if>
          </admst:for-each>
        </admst:if>
      </admst:for-each>
      <admst:text format="#endif /* _DERIVATE */\n"/>
      <admst:text format="#if defined(_DERIVATE)\n"/>
      <admst:for-each select="rhs/function">
        <admst:if test="arguments[count(.)=1]">
          <admst:for-each select="arguments[position(.)=1]">
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
            <admst:choose>
              <admst:when test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_d_%s(__%s_%s,__d_%s_%s,(%s))\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__%s_%s)\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__d_%s_%s)\n"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_%s(__%s_%s,(%s))\n"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:if>
        <admst:if test="arguments[count(.)=2]">
          <admst:value-of select="./position(.)-1"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:text format="_%s(__%s_%s,"/>
          <admst:text test="[name='div']" format="__dFx_%(name)_%(position(.)-1),"/>
          <admst:join select="arguments" separator=",">
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
            <admst:text format="%s"/>
          </admst:join>
          <admst:text format=")\n"/>
          <admst:for-each select="arguments">
            <admst:if test="[position(.)=1]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_dx_%s(__dFx_%s_%s,__%s_%s,"/>
                <admst:join select="../arguments" separator=",">
                  <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
                  <admst:text format="%s"/>
                </admst:join>
                <admst:text format=")\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__dFx_%s_%s)\n"/>
              </admst:if>
            </admst:if>
            <admst:if test="[position(.)=2]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_dy_%s(__dFy_%s_%s,"/>
                <admst:text test="[../name='div']" format="__dFx_%(../name)_%(../position(.)-1),"/>
                <admst:text format="__%s_%s,"/>
                <admst:join select="../arguments" separator=",">
                  <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
                  <admst:text format="%s"/>
                </admst:join>
                <admst:text format=")\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__dFy_%s_%s)\n"/>
              </admst:if>
            </admst:if>
          </admst:for-each>
          <admst:value-of select="position(.)-1"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:text format="EXIT_IF_ISNAN(__%s_%s)\n"/>
        </admst:if>
      </admst:for-each>
      <admst:text format="#else\n"/>
    </admst:if>
  </admst:if>
  <admst:apply-templates select="rhs" match="function:assignment"/>
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)ddx:function:computationrest\n"/>
  <admst:if test="[datatypename='contribution' or lhs/insource='yes']">
  <!--<admst:if test="[lhs/insource='yes']"> -->
    <admst:if test="rhs[hasVoltageDependentFunction='yes']">
      <admst:text format="#endif\n"/>
    </admst:if>
  </admst:if>
</admst:template>

<!-- analog//function: ddx handling -->
<admst:template match="ddx:function:computation">
  <admst:variable name="doit__" select="no"/>
  <admst:if test="[datatypename='contribution' or lhs/insource='yes']">
  <!--<admst:if test="[lhs/insource='yes']"> -->
    <admst:if test="rhs[hasVoltageDependentFunction='yes']">
      <admst:variable name="doit__" select="yes"/>
    </admst:if>
  </admst:if>
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)ddx:function:computation
    ...  . %(.)*
    ...  l %(lhs)*
    ...  r %(rhs)*
    ...  t %(lhs/type)*
    ...  $(doit__) \n "/>

  <!-- old hack...
  <admst:if test="lhs[insource!='no']">
    <admst:if test="rhs[hasVoltageDependentFunction='yes']">
      <admst:variable name="doit__" select="yes"/>
    </admst:if>
  </admst:if>
  -->

  <admst:if test="[$doit__='yes']">
    <admst:if test="[$_DERIVATE='yes']">
      /* derivate stuff */

      <admst:for-each select="rhs/function">
        <admst:if test="arguments[count(.)=1]">
          <admst:for-each select="arguments[position(.)=1]">
            <admst:if test="math[dependency!='constant']">
              <admst:value-of select="../position(.)-1"/>
              <admst:apply-templates select=".." match="function:getname_push"/>
              <admst:text format="double __d_%s_%s = 0.0; // 1st arg of one: %(.)\n"/>
            </admst:if>
          </admst:for-each>
        </admst:if>
        <admst:if test="arguments[count(.)=2]">
          <admst:for-each select="arguments">
            <admst:if test="[position(.)=1]">
              <admst:if test="[(../name='div') or (math/dependency!='constant')]">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="double __dFx_%s_%s = 0.0; // first arg of two: %(.)\n"/>
              </admst:if>
            </admst:if>
            <admst:if test="[position(.)=2]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="double __dFy_%s_%s=0.0; //case 3?\n"/>
              </admst:if>
            </admst:if>
          </admst:for-each>
        </admst:if>
      </admst:for-each>
    </admst:if> <!-- _DERIVATE -->
	 <admst:if test="[$_DERIVATE!='yes']">/*no _DERIVATE*/</admst:if>
	 <admst:if test="[$_DERIVATE='yes']"> // derivate... %(.)

      <admst:for-each select="rhs/function">
        <admst:if test="arguments[count(.)=1]">
          <admst:for-each select="arguments[position(.)=1]">
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
            <admst:choose>
              <admst:when test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_d_%s(__%s_%s,__d_%s_%s,(%s))\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__%s_%s);\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__d_%s_%s);\n"/>
              </admst:when>
              <admst:otherwise>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_%s(__%s_%s,(%s))\n"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:if>
        <admst:if test="arguments[count(.)=2]">
          <admst:value-of select="./position(.)-1"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:text format="_%s(__%s_%s,"/>
          <admst:text test="[name='div']" format="__dFx_%(name)_%(position(.)-1),"/>
          <admst:join select="arguments" separator=",">
            <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
				<admst:text format="%s /**/"/>
          </admst:join>
			 <admst:text format=") // 3592\n"/>
          <admst:for-each select="arguments">
            <admst:if test="[position(.)=1]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_dx_%s(__dFx_%s_%s,__%s_%s,"/>
                <admst:join select="../arguments" separator=",">
                  <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
                  <admst:text format="%s"/>
                </admst:join>
					 <admst:text format=") // 3606\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__dFx_%s_%s);\n"/>
              </admst:if>
            </admst:if>
            <admst:if test="[position(.)=2]">
              <admst:if test="math[dependency!='constant']">
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="_dy_%s(__dFy_%s_%s,"/>
                <admst:text test="[../name='div']" format="__dFx_%(../name)_%(../position(.)-1),"/>
                <admst:text format="__%s_%s,"/>
                <admst:join select="../arguments" separator=",">
                  <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
                  <admst:text format="%s"/>
                </admst:join>
                <admst:text format=")\n"/>
                <admst:value-of select="../position(.)-1"/>
                <admst:apply-templates select=".." match="function:getname_push"/>
                <admst:text format="EXIT_IF_ISNAN(__dFy_%s_%s);\n"/>
              </admst:if>
            </admst:if>
          </admst:for-each>
          <admst:value-of select="position(.)-1"/>
          <admst:apply-templates select="." match="function:getname_push"/>
          <admst:text format="EXIT_IF_ISNAN(__%s_%s);\n"/>
        </admst:if>
      </admst:for-each>
    </admst:if> <!-- _DERIVATE -->
  </admst:if>  <!-- 2 conditions -->
  <admst:if test="[$_DERIVATE='no' or $doit__='no']">
    // function assignment %(.)

    <admst:apply-templates select="rhs" match="function:assignment"/>
  </admst:if>
  <admst:if test="[$_DERIVATE!='no' and $doit__!='no']">
    // function assignment %(.)

#if 0

	 <admst:apply-templates select="rhs" match="function:assignment"/>
#endif
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//assignment -->
<admst:template match="assignment">
// assignment ngspiceVersion.xml %(.)

  <admst:text test="rhs[not(nilled(function[class='builtin']))]" format="{\n"/>

  <admst:for-each select="rhs/function">
    <admst:value-of select="position(.)-1"/>
    <admst:apply-templates select="." match="function:getname_push"/>
    <admst:text format="double __%s_%s=0.0;\n"/>
  </admst:for-each>

  <admst:apply-templates select="." match="ddx:function:computation"/>

  <admst:if test="lhs[prototype/insource='yes']">
    <admst:if test="rhs/probe">

// iterating rhs probes

      <admst:for-each select="rhs/probe">
// probe %(.)

        <admst:value-of select="."/>
        <admst:value-to select="/simulator/probe" value="%p"/>
        <admst:apply-templates select="../tree" match="subexpression:differentiate"/>
        <admst:value-of select="/simulator/ddx"/>
        <admst:value-of select="branch/nnode/name"/>
        <admst:value-of select="branch/pnode/name"/>
        <admst:value-of select="nature/access"/>
        <admst:value-of select="../../lhs/name"/>
        <admst:text format="%s_%s%s_%s=%s;\n"/>
        <admst:value-of select="branch/nnode/name"/>
        <admst:value-of select="branch/pnode/name"/>
        <admst:value-of select="nature/access"/>
        <admst:value-of select="../../lhs/name"/>
        <admst:text format="EXIT_IF_ISNAN(%s_%s%s_%s); // 3515\n"/>
      </admst:for-each>

    </admst:if>
  </admst:if>

  <admst:apply-templates select="lhs" match="variable:lhs" required="yes"/>
  <admst:apply-templates select="rhs" match="expression:stringify:noprobe"/>
  <admst:text format="=%s;\n"/>
  <admst:text format="EXIT_IF_ISNAN("/>
  <admst:apply-templates select="lhs" match="variable:lhs" required="yes"/>
  <admst:text format="); // 3526\n"/>

  <admst:if test="lhs[insource='yes']">
    <admst:value-of select="rhs/probe"/>
    <admst:if-inside select="lhs/probe" list="%p">
      <admst:if test="lhs/probe">

        <admst:for-each select="lhs/probe">
          <admst:value-of select="../../rhs/probe"/>
          <admst:if-not-inside select="." list="%p">
            <admst:value-of select="branch/nnode/name"/>
            <admst:value-of select="branch/pnode/name"/>
            <admst:value-of select="nature/access"/>
            <admst:value-of select="../name"/>
            <admst:text format="%s_%s%s_%s=0.0;\n"/>
          </admst:if-not-inside>
        </admst:for-each>

      </admst:if>
    </admst:if-inside>
  </admst:if>

  <admst:text test="rhs[not(nilled(function[class='builtin']))]" format="}\n"/>

  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)assignment<-*\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- ---------- analog//tr_eval_kept ------------------------------- -->
<admst:template match="tr_eval_kept">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)tr_eval_kept*\n"/>
  <admst:variable name="eval_kept" select="yes"/>
  <admst:variable name="keep_ic" select="no"/>


  <admst:text format="// var decl\n"/>
  <admst:apply-templates select="module/evaluation/variable" match="variable:declaration"/>
  <admst:text format="// var decl done\n"/>


  <admst:reset select="module/evaluation/variable"/>
  <admst:for-each select="item">
    <admst:choose>
      <admst:when test="adms[datatypename!='block']">
        <admst:text format="// not a block\n"/>
        <admst:text format="// fetch states\n"/>
        <admst:apply-templates select="." match="%(adms/datatypename)"/>
        <admst:text format="// fetched states\n"/>
        <admst:text format="// not a block done\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:if test="[name!='initial_model' and name!='initial_instance']">
          <admst:text format="// not initial_\n"/>

          <admst:text format="{ //block \n"/>
          <admst:text select="item" format=" // %(adms/datatypename) \n"/>
          <admst:apply-templates select="item" match="%(adms/datatypename)"/>
          <admst:text format="} //blovk... \n"/>

          <admst:text format="// done block\n"/>
        </admst:if>
      </admst:otherwise>
    </admst:choose>
  </admst:for-each>
  <admst:variable name="eval_kept" select="no"/>
  <admst:variable name="keep_ic" select="no"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- ---------- analog//keep_ic ------------------------------- -->
<admst:template match="keep_ic">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)keep_ic*\n"/>
  <admst:message format=" //keep_ic have here: %(adms/datatypename) name %(name)\n"/>
  <admst:for-each select="item">
    <!-- <admst:message format=" //keep_ic looking around: have here: %(adms/datatypename) name %(name): %(.)\n"/>
    -->
  </admst:for-each>

  <admst:variable name="keep_ic" select="yes"/>
  <admst:apply-templates select="module/evaluation/variable" match="variable:declaration"/>
  <admst:reset select="module/evaluation/variable"/>
  <admst:for-each select="item">
    <admst:choose>
      <admst:when test="adms[datatypename!='block']">
        <admst:text format="// not a block\n"/>
        <admst:warning format="keep ic takes a block, not %(adms/datatypename) - %(.)\n"/>
        <!-- <admst:apply-templates select="." match="%(adms/datatypename)"/> -->
        <admst:text format="// not a block done\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:if test="[name!='initial_model' and name!='initial_instance']">
          <admst:text format="// not initial_\n"/>
          <!-- <admst:apply-templates select="." match="block"/> -->

          <admst:text format="{ //block \n"/>


          <admst:apply-templates select="item" match="%(adms/datatypename)"/>
          <admst:text format="} //blovk. \n"/>

          <admst:text format="// done block\n"/>
        </admst:if>
      </admst:otherwise>
    </admst:choose>
  </admst:for-each>
  <admst:variable name="keep_ic" select="no"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- ---------- analog//contribution ------------------------------- -->
<admst:template match="ac_load">
  <admst:message test="[/dbg_xml='yes']" format="*ac_load*\n"/>
  incomplete();
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- ---------- analog//contribution ------------------------------- -->

<admst:template match="contribution">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)contribution*\n"/>
  <admst:choose>
    <admst:when test="[whitenoise='no' and flickernoise='no']">
      <admst:apply-templates select="." match="contribution:nonoise" required="yes"/>
    </admst:when>
    <admst:otherwise>
      <admst:apply-templates select="." match="contribution:noise" required="yes"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- ----------------------  analog//conditional ------------------ -->
<admst:template match="conditional">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)conditional*\n"/>
  <admst:if test="if[dynamic='yes']">
    <admst:choose>
      <admst:when test="[nilled(else)]">
        <admst:text format="#ifdef _DYNAMIC /*&lt;dynamic_ifthen&gt;*/\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:text format="#ifdef _DYNAMIC /*&lt;dynamic_ifthenelse&gt;*/\n"/>
      </admst:otherwise>
    </admst:choose>
  </admst:if>
  <admst:if test="if[not(nilled(function[class='builtin']))]">
    <admst:text format="{\n"/>
    <admst:for-each select="if/function">
      <admst:value-of select="position(.)-1"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:text format="double __%s_%s=0.0;\n"/>
    </admst:for-each>
    <admst:apply-templates select="if" match="function:assignment"/>
  </admst:if>
  <admst:apply-templates select="if" match="expression:stringify:noprobe"/>
  <admst:text format="if\n(%s)\n"/>
  <admst:if test="then/adms[datatypename!='block']">
    <admst:text format="{\n"/>
  </admst:if>
  <admst:value-of select="then/adms/datatypename"/>
  <admst:apply-templates select="then" match="%s" required="yes"/>
  <admst:if test="then/adms[datatypename!='block']">
    <admst:text format="}\n"/>
  </admst:if>
  <admst:if test="else">
    <admst:text format="else\n"/>
    <admst:choose>
      <admst:when test="else/adms[datatypename='block']">
        <admst:value-of select="else/adms/datatypename"/>
        <admst:apply-templates select="else" match="%s" required="yes"/>
      </admst:when>
      <admst:otherwise>
        <admst:text format="{\n"/>
        <admst:value-of select="else/adms/datatypename"/>
        <admst:apply-templates select="else" match="%s" required="yes"/>
        <admst:text format="}\n"/>
      </admst:otherwise>
    </admst:choose>
  </admst:if>
  <admst:if test="if[not(nilled(function[class='builtin']))]">
    <admst:text format="}\n"/>
  </admst:if>
  <admst:if test="if[dynamic='yes']">
    <admst:choose>
      <admst:when test="[nilled(else)]">
        <admst:text format="#endif /*&lt;/dynamic_ifthen&gt;*/\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:text format="#endif /*&lt;/dynamic_ifthenelse&gt;*/\n"/>
      </admst:otherwise>
    </admst:choose>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//case -->
<admst:template match="case">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)case*\n"/>
  <admst:error format="case statement: please implement me! (inside block)\n"/>
  <admst:text format="/*CASE*/;\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//nilled -->
<admst:template match="nilled">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)nilled*\n"/>
  <admst:text format=";\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//whileloop -->
<admst:template match="whileloop">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)whileloop*\n"/>
  <admst:if test="while[dynamic='yes']">
    <admst:text format="#ifdef _DYNAMIC /*&lt;dynamic_while&gt;*/\n"/>
  </admst:if>
  <admst:if test="while[not(nilled(function[class='builtin']))]">
    <admst:text format="{\n"/>
    <admst:for-each select="while/function">
      <admst:value-of select="position(.)-1"/>
      <admst:apply-templates select="." match="function:getname_push"/>
      <admst:text format="double __%s_%s=0.0;\n"/>
    </admst:for-each>
    <admst:apply-templates select="while" match="function:assignment"/>
  </admst:if>
  <admst:apply-templates select="while" match="expression:stringify:noprobe"/>
  <admst:text format="while\n(%s)\n"/>
  <admst:if test="whileblock/adms[datatypename!='block']">
    <admst:text format="{\n"/>
  </admst:if>
  <admst:value-of select="whileblock/adms/datatypename"/>
  <admst:apply-templates select="whileblock" match="%s" required="yes"/>
  <admst:if test="whileblock/adms[datatypename!='block']">
    <admst:text format="}\n"/>
  </admst:if>
  <admst:if test="while[not(nilled(function[class='builtin']))]">
    <admst:text format="}\n"/>
  </admst:if>
  <admst:if test="while[dynamic='yes']">
    <admst:text format="#endif /*&lt;/dynamic_while&gt;*/\n"/>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- expression//probe -->
<admst:template match="probe">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)probe*\n"/>
  <admst:choose>
    <admst:when test="branch/nnode[grounded='no']">
      <admst:value-of select="branch/nnode/name"/>
      <admst:value-of select="branch/pnode/name"/>
      <admst:value-to select="/simulator/tmp" value="BP(%s,%s)"/>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="branch/pnode/name"/>
      <admst:value-to select="/simulator/tmp" value="NP(n_%s)"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- expression//node -->
<admst:template match="node">
<admst:message test="[/dbg_xml='yes']" format="*(ngspvers)node*\n"/>
  <admst:value-of select="name"/>
  <admst:error format="module node not expected here ... %s\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- expression//string -->
<admst:template match="string">
<admst:message test="[/dbg_xml='yes']" format="*(ngspvers)string*\n"/>
  <admst:value-of select="value"/>
  <admst:value-to select="/simulator/tmp" value="&quot;%s&quot;"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- expression//number -->
<admst:template match="number">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)number*\n"/>
  <admst:choose>
    <admst:when test="[scalingunit='1']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="%s"/>
    </admst:when>
    <admst:when test="[scalingunit='E']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+18)"/>
    </admst:when>
    <admst:when test="[scalingunit='P']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+15)"/>
    </admst:when>
    <admst:when test="[scalingunit='T']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+12)"/>
    </admst:when>
    <admst:when test="[scalingunit='G']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+9)"/>
    </admst:when>
    <admst:when test="[scalingunit='M']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+6)"/>
    </admst:when>
    <admst:when test="[scalingunit='k']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+3)"/>
    </admst:when>
    <admst:when test="[scalingunit='h']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+2)"/>
    </admst:when>
    <admst:when test="[scalingunit='D']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e+1)"/>
    </admst:when>
    <admst:when test="[scalingunit='d']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-1)"/>
    </admst:when>
    <admst:when test="[scalingunit='c']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-2)"/>
    </admst:when>
    <admst:when test="[scalingunit='m']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-3)"/>
    </admst:when>
    <admst:when test="[scalingunit='u']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-6)"/>
    </admst:when>
    <admst:when test="[scalingunit='n']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-9)"/>
    </admst:when>
    <admst:when test="[scalingunit='A']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-10)"/>
    </admst:when>
    <admst:when test="[scalingunit='p']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-12)"/>
    </admst:when>
    <admst:when test="[scalingunit='f']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-15)"/>
    </admst:when>
    <admst:when test="[scalingunit='a']">
      <admst:value-of select="value"/>
      <admst:value-to select="/simulator/tmp" value="(%s*1.0e-18)"/>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="scalingunit"/>
      <admst:error format="scaling unit not supported: %s\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- analog//contribution[noise] -->
<admst:template match="contribution:noise">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)contribution:noise*\n"/>
  <admst:fatal format="template not in use: remved\n"/>
  <admst:if test="[flickernoise='yes']">
    <admst:text format="ngspice_flickernoise(%(lhs/branch/pnode/name),%(lhs/branch/nnode/name)"/>
    <admst:for-each select="rhs/tree/arguments">
      <admst:apply-templates select="." match="%(datatypename)"/>
      <admst:value-of select="/simulator/tmp"/>
      <admst:text format=",%s"/>
    </admst:for-each>
    <admst:text test="[count(rhs/tree/arguments)=2]" format=",NULL"/>
  </admst:if>
  <admst:if test="[whitenoise='yes']">
    <admst:text format="ngspice_whitenoise(%(lhs/branch/pnode/name),%(lhs/branch/nnode/name)"/>
    <admst:for-each select="rhs/tree/arguments">
      <admst:apply-templates select="." match="%(datatypename)"/>
      <admst:value-of select="/simulator/tmp"/>
      <admst:text format=",%s"/>
    </admst:for-each>
    <admst:text test="[count(rhs/tree/arguments)=1]" format=",NULL"/>
  </admst:if>
  <admst:text format=")\n"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- variable:rhs -->
<admst:template match="variable">
  <admst:value-of select="name"/>
  <admst:if test="[parametertype='analogfunction']">
    <admst:value-to select="/simulator/tmp" value="%s"/>
  </admst:if>
  <admst:if test="[input='yes' and parametertype='model']">
    <admst:value-to select="/simulator/tmp" value="m-&gt;$(var_prefix)%s"/>
  </admst:if>
  <admst:if test="[input='yes' and parametertype='instance']">
    <admst:value-to select="/simulator/tmp" value="/*instance*/&gt;$(var_prefix)%s"/>
  </admst:if>
  <admst:if test="[input='no' and scope='global_model']">
    <admst:value-to select="/simulator/tmp" value="m-&gt;$(var_prefix)%s"/>
  </admst:if>
  <admst:if test="[input='no' and scope='global_instance']">
    <admst:value-to select="/simulator/tmp" value="/*gli*/$(var_prefix)%s"/>
  </admst:if>
  <admst:if test="[parametertype!='analogfunction' and scope='local']">
    <admst:value-to select="/simulator/tmp" value="/*else*/%s"/>
  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- variable:lhs -->
<admst:template match="variable:lhs">
  <admst:message test="[/dbg_xml='yes']" format="*(ngspvers)variable:lhs*\n"/>
  <admst:text test="[input='yes' and parametertype='model']" format="m->$(var_prefix)%(name)"/>
  <admst:text test="[input='yes' and parametertype='instance']" format="/*i*/$(var_prefix)%(name)"/>
  <admst:text test="[input='no' and scope='global_model']" format="/*glm*/m->$(var_prefix)%(name)"/>
  <admst:text test="[input='no' and scope='global_instance']" format="/*gli*/$(var_prefix)%(name)"/>
  <admst:text test="[scope='local']" format="%(name)"/>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- handle analog//callfunctions -->
<admst:template match="callfunction">
  <admst:message format="callfunction: untested\n"/>
  <admst:choose>
    <admst:when test="function[name='\$strobe']">
      <admst:text format="fprintf(stdout"/>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="function/name"/>
#warning &quot;%s not supported by this interface&quot;
// FIXME(
    </admst:otherwise>
  </admst:choose>
  <admst:for-each select="function/arguments">
    <admst:apply-templates select="./tree" match="expression:stringify:noprobe"/>
    <admst:text format=",%s"/>
  </admst:for-each>
  <admst:text format=");\n"/>
  <admst:choose>
    <admst:when test="function[name='\$strobe']">
      <admst:text format="fprintf(stdout,&quot;\\n&quot;);\n"/>
    </admst:when>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="generate_c:declare_variable">
  <admst:message test="[/dbg_xml='yes']" format="*generate_c:declare_variable*\n"/>
  <admst:if test="[type='integer']">\tint %(name);</admst:if>
  <admst:if test="[type='real']">\tdouble %(name) = std::numeric_limits&lt;double&gt;::quiet_NaN();</admst:if>
  <admst:if test="[type='string']">\tchar* %(name);</admst:if>
      _used(%(name));

</admst:template>

<!-- *********************************************************************************** -->
<!-- *********************************** analog function ******************************* -->
<!-- *********************************************************************************** -->

<admst:template match="afunction:getname"> <!--same as function:getname?-->
  <admst:message test="[/dbg_xml='yes']" format="*(af)afunction:getname*\n"/>
  <admst:choose>
    <admst:when test="[name='abs']"><admst:return name="afunction:getname" string="abs"/></admst:when>
    <admst:when test="[name='\$shrinkl']"><admst:return name="afunction:getname" string="shrinkl"/></admst:when>
    <admst:when test="[name='\$shrinka']"><admst:return name="afunction:getname" string="shrinka"/></admst:when>
    <admst:when test="[name='log']"><admst:return name="afunction:getname" string="log10"/></admst:when>
    <admst:when test="[name='ln']"><admst:return name="afunction:getname" string="logE"/></admst:when>
    <admst:when test="[name='limexp']"><admst:return name="afunction:getname" string="limexp"/></admst:when>
    <admst:when test="[name='\$limexp']"><admst:return name="afunction:getname" string="limexp"/></admst:when>
    <admst:when test="[name='\$vt']"><admst:return name="afunction:getname" string="vt"/></admst:when>
    <admst:when test="[name='\$model']"><admst:return name="afunction:getname" string="_modelname"/></admst:when>
    <admst:when test="[name='\$instance']"><admst:return name="afunction:getname" string="_instancename"/></admst:when>
    <admst:when test="[name='\$temperature']">
      <admst:return name="afunction:getname" string="_ambient_temp"/>
    </admst:when>
    <admst:when test="[name='\$nominal_temperature']"><admst:return name="afunction:getname" string="_circuit_tnom"/></admst:when>
    <admst:otherwise><admst:return name="afunction:getname" string="%(name)"/></admst:otherwise>
  </admst:choose>
</admst:template>

<!-- -------------------------------------------------------------------- -->
<admst:template match="af:print:expression">
  <admst:message test="[/dbg_xml='yes']" format="*(af)af:print:expression*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='expression']">
      <admst:apply-templates select="tree" match="af:print:expression">
        <admst:variable name="expression" select="%(returned('x')/value)"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:variable name="dx_%(name)" select="%(returned('dx.%(name)')/value)"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:return name="x" string="$expression"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:return name="dx.%(name)" string="$(dx_%(name))"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='probe']">
      <admst:fatal format="probe not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='variable']">
      <admst:variable name="variable" select="%(name)"/>
      <admst:return name="x" string="$variable"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:variable name="ddx" select="%(name)"/>
        <admst:choose>
          <admst:when test="[$variable='$ddx']">
            <admst:return name="dx.$ddx" string="1.0"/>
          </admst:when>
          <admst:when test="../..[input='yes']">
            <admst:return name="dx.$ddx" string="0.0"/>
          </admst:when>
          <admst:otherwise>
            <admst:return name="dx.$ddx" string="$(variable)_$ddx"/>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>

    </admst:when>
    <admst:when test="adms[datatypename='mapply_unary']">
      <admst:if test="[name='plus']">
        <admst:variable name="op" select="+"/>
      </admst:if>
      <admst:if test="[name='minus']">
        <admst:variable name="op" select="-"/>
      </admst:if>
      <admst:if test="[name='not']">
        <admst:variable name="op" select="!"/>
      </admst:if>
      <admst:if test="[name='bw_not']">
        <admst:variable name="op" select="~"/>
      </admst:if>
      <admst:apply-templates select="arg1" match="af:print:expression"><admst:variable name="arg1" select="%(returned('x')/value)"/></admst:apply-templates>
      <admst:return name="x" string="($op$arg1)"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:return name="dx.%(name)" string="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='mapply_binary']">
      <admst:apply-templates select="arg1" match="af:print:expression">
        <admst:variable name="x" select="%(returned('x')/value)"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:variable name="dx_%(name)" select="%(returned('dx.%(name)')/value)"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:apply-templates select="arg2" match="af:print:expression">
        <admst:variable name="y" select="%(returned('x')/value)"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:variable name="dy_%(name)" select="%(returned('dx.%(name)')/value)"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:choose>
        <admst:when test="[name='addp']">
          <admst:choose>
            <admst:when test="[$x='0.0' and $y='0.0']">
              <admst:return name="x" string="0.0"/>
            </admst:when>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" string="(+$y)"/>
            </admst:when>
            <admst:when test="[$y='0.0']">
              <admst:return name="x" string="$x"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" string="($x+$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="df" select="%(name)"/>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:when>
              <admst:when test="[$y='0.0']">
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" string="(+$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" string="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" string="($dx+$dy)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='addm']">
          <admst:choose>
            <admst:when test="[$x='0.0' and $y='0.0']">
              <admst:return name="x" string="0.0"/>
            </admst:when>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" string="(-$y)"/>
            </admst:when>
            <admst:when test="[$y='0.0']">
              <admst:return name="x" string="$x"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" string="($x-$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="df" select="%(name)"/>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:when>
              <admst:when test="[$y='0.0']">
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" string="(-$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" string="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" string="($dx-$dy)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='multtime']">
          <admst:choose>
            <admst:when test="[$x='0.0' or $y='0.0']">
              <admst:return name="x" string="0.0"/>
            </admst:when>
            <admst:when test="[$x='1.0' and $y='1.0']">
              <admst:return name="x" string="1.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" string="($x*$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="df" select="%(name)"/>
            <admst:choose>
              <admst:when test="[$x='0.0' or $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0' and $y='1.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='1.0']">
                <admst:return name="dx.$df" string="($x)"/>
              </admst:when>
              <admst:when test="[$dx='1.0' and $dy='0.0']">
                <admst:return name="dx.$df" string="($y)"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" string="($x*$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" string="$dx*$y"/>
              </admst:when>
              <admst:when test="[$dx='1.0' and $dy='1.0']">
                <admst:return name="dx.$df" string="($x+$y)"/>
              </admst:when>
              <admst:when test="[$dx='1.0']">
                <admst:return name="dx.$df" string="($y+($dy*$x))"/>
              </admst:when>
              <admst:when test="[$dy='1.0']">
                <admst:return name="dx.$df" string="($dx*$y)+$x"/>
              </admst:when>
              <admst:when test="[$x='1.0']">
                <admst:return name="dx.$df" string="$dy"/>
              </admst:when>
              <admst:when test="[$y='1.0']">
                <admst:return name="dx.$df" string="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" string="(($dx*$y)+($x*$dy))"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='multdiv']">
          <admst:choose>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" string="0.0"/>
            </admst:when>
            <admst:when test="[$x='1.0' and $y='1.0']">
              <admst:return name="x" string="1.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" string="($x/$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="df" select="%(name)"/>
            <admst:choose>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0' and $y='1.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$(df))"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$x='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" string="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0']">
                <admst:choose>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" string="(-1/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" string="(-$dy/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:choose>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" string="(-$x/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" string="(-($x*$dy)/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:when test="[$dx='1.0']">
                <admst:choose>
                  <admst:when test="[$dy='0.0']">
                    <admst:return name="dx.$df" string="(1/$y)"/>
                  </admst:when>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" string="(($y-$x)/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" string="(($y-($x*$dy))/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:otherwise>
                <admst:choose>
                  <admst:when test="[$y='1.0']">
                    <admst:return name="dx.$df" string="$dx"/>
                  </admst:when>
                  <admst:when test="[$dy='0.0']">
                    <admst:return name="dx.$df" string="$dx/$y"/>
                  </admst:when>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" string="(($dx*$y)-$x)/($y*$y)"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" string="($dx*$y-$x*$dy)/($y*$y)"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:otherwise>
          <admst:choose>
            <admst:when test="[name='bw_equr']">
              <admst:return name="x" string="($x^~$y)"/>
            </admst:when>
            <admst:when test="[name='bw_equl']">
              <admst:return name="x" string="($x~^$y)"/>
            </admst:when>
            <admst:when test="[name='bw_xor']">
              <admst:return name="x" string="($x^$y)"/>
            </admst:when>
            <admst:when test="[name='bw_or']">
              <admst:return name="x" string="($x|$y)"/>
            </admst:when>
            <admst:when test="[name='bw_and']">
              <admst:return name="x" string="($x&amp;$y)"/>
            </admst:when>
            <admst:when test="[name='or']">
              <admst:return name="x" string="($x||$y)"/>
            </admst:when>
            <admst:when test="[name='and']">
              <admst:return name="x" string="($x&amp;&amp;$y)"/>
            </admst:when>
            <admst:when test="[name='equ']">
              <admst:return name="x" string="($x==$y)"/>
            </admst:when>
            <admst:when test="[name='multmod']">
              <admst:return name="x" string="((int)$x%%(int)$y)"/>
            </admst:when>
            <admst:when test="[name='notequ']">
              <admst:return name="x" string="($x!=$y)"/>
            </admst:when>
            <admst:when test="[name='lt']">
              <admst:return name="x" string="($x&lt;$y)"/>
            </admst:when>
            <admst:when test="[name='lt_equ']">
              <admst:return name="x" string="($x&lt;=$y)"/>
            </admst:when>
            <admst:when test="[name='gt']">
              <admst:return name="x" string="($x&gt;$y)"/>
            </admst:when>
            <admst:when test="[name='gt_equ']">
              <admst:return name="x" string="($x&gt;=$y)"/>
            </admst:when>
            <admst:when test="[name='shiftr']">
              <admst:return name="x" string="($x&gt;&gt;$y)"/>
            </admst:when>
            <admst:when test="[name='shiftl']">
              <admst:return name="x" string="($x&lt;&lt;$y)"/>
            </admst:when>
            <admst:otherwise>
              <admst:error format="%(name): function not handled\n"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:return name="dx.%(name)" string="0.0"/>
          </admst:for-each>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="adms[datatypename='mapply_ternary']">
      <admst:apply-templates select="arg1" match="af:print:expression">
        <admst:variable name="x" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:if test="[name='conditional']">
        <admst:apply-templates select="arg2" match="af:print:expression">
          <admst:variable name="y" select="%(returned('x')/value)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="dy_%(name)" select="%(returned('dx.%(name)')/value)"/>
          </admst:for-each>
        </admst:apply-templates>
        <admst:apply-templates select="arg3" match="af:print:expression">
          <admst:variable name="z" select="%(returned('x')/value)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="dz_%(name)" select="%(returned('dx.%(name)')/value)"/>
          </admst:for-each>
        </admst:apply-templates>
        <admst:return name="x" string="($x?$y:$z)"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:variable name="df" select="%(name)"/>
          <admst:return name="dx.$df" string="($x?$(dy_$df):$(dz_$df))"/>
        </admst:for-each>
      </admst:if>
    </admst:when>

    <admst:when test="adms[datatypename='function']">
      <admst:apply-templates select="." match="afunction:getname">
        <admst:variable name="function" select="%(returned('afunction:getname')/value)"/>
      </admst:apply-templates>
      <admst:variable name="args" select=""/>
      <admst:for-each select="arguments">
        <admst:if test="[not($args='')]">
          <admst:variable name="args" select="$args,"/>
        </admst:if>
        <admst:apply-templates select="." match="af:print:expression">
          <admst:variable name="index" select="%(index(../arguments,.))"/>
          <admst:variable name="args" select="$args%(returned('x')/value)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="arg$(index)_%(name)" select="%(returned('dx.%(name)')/value)"/>
          </admst:for-each>
        </admst:apply-templates>
      </admst:for-each>
      <admst:choose>
        <admst:when test="[ name='cos' or name='sin' or name='tan' or name='cosh' or name='sinh' or name='tanh' or name='acos' or name='asin'
                            or name='atan' or name='ln' or name='log' or name='exp' or name='sqrt' or name='abs' or name='limexp'
                            or name='pow' or name='hypot' or name='min' or name='max' ]">
          <admst:return name="x" string="_$function($args)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="name" path="name"/>
            <admst:variable name="ret" select=""/>
            <admst:for-each select="../../arguments">
              <admst:if test="[not($ret='')]">
                <admst:variable name="ret" select="$ret+"/>
              </admst:if>
              <admst:variable name="index" select="%(index(../arguments,.))"/>
              <admst:variable name="ret" select="$(ret)_d$(index)_$function($args)*($(arg$(index)_$name))"/>
            </admst:for-each>
            <admst:return name="dx.$name" string="$ret"/>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='ceil' or name='floor']">
          <admst:return name="x" string="$function($args)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="name" path="name"/>
            <admst:return name="dx.$name" string="0.0"/>
          </admst:for-each>
        </admst:when>
        <admst:otherwise>
          <admst:return name="x" string="$(module)_$function($args)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:variable name="name" path="name"/>
            <admst:variable name="darg" select=""/>
            <admst:for-each select="../../arguments">
              <admst:variable name="index" select="%(index(../arguments,.))"/>
              <admst:variable name="darg" select="$darg,($(arg$(index)_$name))"/>
            </admst:for-each>
            <admst:return name="dx.$name" string="$(module)_d_$function($args$darg)"/>
          </admst:for-each>
        </admst:otherwise>
      </admst:choose>
    </admst:when>

    <admst:when test="adms[datatypename='string']">
      <admst:return name="x" string="&quot;%(value)&quot;"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:return name="dx.%(name)" string="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='number']">
      <admst:choose>
        <admst:when test="[scalingunit='1']">
          <admst:return name="x" string="%(value)"/>
        </admst:when>
        <admst:when test="[scalingunit='E']">
          <admst:return name="x" string="(%(value)*1.0e+18)"/>
        </admst:when>
        <admst:when test="[scalingunit='P']">
          <admst:return name="x" string="(%(value)*1.0e+15)"/>
        </admst:when>
        <admst:when test="[scalingunit='T']">
          <admst:return name="x" string="(%(value)*1.0e+12)"/>
        </admst:when>
        <admst:when test="[scalingunit='G']">
          <admst:return name="x" string="(%(value)*1.0e+9)"/>
        </admst:when>
        <admst:when test="[scalingunit='M']">
          <admst:return name="x" string="(%(value)*1.0e+6)"/>
        </admst:when>
        <admst:when test="[scalingunit='k']">
          <admst:return name="x" string="(%(value)*1.0e+3)"/>
        </admst:when>
        <admst:when test="[scalingunit='h']">
          <admst:return name="x" string="(%(value)*1.0e+2)"/>
        </admst:when>
        <admst:when test="[scalingunit='D']">
          <admst:return name="x" string="(%(value)*1.0e+1)"/>
        </admst:when>
        <admst:when test="[scalingunit='d']">
          <admst:return name="x" string="(%(value)*1.0e-1)"/>
        </admst:when>
        <admst:when test="[scalingunit='c']">
          <admst:return name="x" string="(%(value)*1.0e-2)"/>
        </admst:when>
        <admst:when test="[scalingunit='m']">
          <admst:return name="x" string="(%(value)*1.0e-3)"/>
        </admst:when>
        <admst:when test="[scalingunit='u']">
          <admst:return name="x" string="(%(value)*1.0e-6)"/>
        </admst:when>
        <admst:when test="[scalingunit='n']">
          <admst:return name="x" string="(%(value)*1.0e-9)"/>
        </admst:when>
        <admst:when test="[scalingunit='A']">
          <admst:return name="x" string="(%(value)*1.0e-10)"/>
        </admst:when>
        <admst:when test="[scalingunit='p']">
          <admst:return name="x" string="(%(value)*1.0e-12)"/>
        </admst:when>
        <admst:when test="[scalingunit='f']">
          <admst:return name="x" string="(%(value)*1.0e-15)"/>
        </admst:when>
        <admst:when test="[scalingunit='a']">
          <admst:return name="x" string="(%(value)*1.0e-18)"/>
        </admst:when>
        <admst:otherwise>
          <admst:fatal format="scaling unit not supported: %(scalingunit)\n"/>
        </admst:otherwise>
      </admst:choose>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:return name="dx.%(name)" string="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="%(datatypename): not handled inside expression\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- -------------------------------------------------------------------- -->
<admst:template match="af:print">
  <admst:message test="[/dbg_xml='yes']" format="*(af)af:print*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='callfunction']">
      <admst:choose>
        <admst:when test="function[name='\$strobe']">
          <admst:variable name="outputfile" select="stdout"/>
        </admst:when>
      </admst:choose>
      <admst:variable name="args" select=""/>
      <admst:for-each select="function/arguments">
        <admst:apply-templates select="." match="af:print:expression">
          <admst:variable name="index" select="%(index(../arguments,.))"/>
          <admst:variable name="args" select="$args,%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="fprintf($outputfile$args); fprintf($outputfile,&quot;\\n&quot;);\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="af:print:expression">
        <admst:variable name="whileblock" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:apply-templates select="while" match="af:print">
        <admst:variable name="while" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:return name="x" string="while($whileblock)\n$while"/>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="if" match="af:print:expression">
        <admst:variable name="if" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:apply-templates select="then" match="af:print">
        <admst:variable name="then" select="%(returned('x')/value)"/>
        </admst:apply-templates>
      <admst:if test="else">
        <admst:apply-templates select="else" match="af:print">
          <admst:variable name="then" select="$(then)else\n%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:if>
      <admst:return name="x" string="if($if)\n$then"/>
    </admst:when>
    <admst:when test="adms[datatypename='case']">
      <admst:apply-templates select="case" match="af:print:expression">
        <admst:variable name="case" select="switch ((int)%(returned('x')/value)) {\n"/>
      </admst:apply-templates>
      <admst:for-each select="caseitem">
        <admst:variable name="condition" select=""/>
        <admst:for-each select="condition">
          <admst:variable name="condition" select="$condition case %(.):"/>
        </admst:for-each>
        <admst:variable name="case" select="$case $condition"/>
        <admst:if test="[defaultcase='yes']">
          <admst:variable name="case" select="$case default:"/>
        </admst:if>
        <admst:variable name="case" select="$case \n"/>
        <admst:apply-templates select="code" match="af:print">
          <admst:variable name="case" select="$case%(returned('x')/value) break;\n"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="$case }"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']">
      <admst:fatal format="contribution not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='assignment']">
      <admst:apply-templates select="rhs" match="af:print:expression">
        <admst:return name="x" string="%(../lhs/name)=%(returned('x')/value);\n"/>
      </admst:apply-templates>
    </admst:when>
    <admst:when test="adms[datatypename='nilled']">
      <admst:return name="x" string=";"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:variable name="block" select=""/>
      <admst:for-each select="item">
        <admst:apply-templates select="." match="af:print">
          <admst:variable name="block" select="$block%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="{$block}"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="%(datatypename): not handled inside blocks\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>


<admst:template match="v2c:converttype">
  <admst:choose>
    <admst:when test="[type='integer']">
      <admst:text format="int"/>
    </admst:when>
    <admst:when test="[type='real']">
      <admst:text format="double"/>
    </admst:when>
    <admst:when test="[type='string']">
      <admst:text format="char*"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="should not be reached\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<admst:template match="analogfunctionH">
  <admst:for-each select="/module/analogfunction">
    <admst:variable name="function" select="%(name)"/>
\t\t
    <admst:apply-templates select="." match="v2c:converttype"/> $function (
    <admst:join select="variable[input='yes']" separator=", ">
      <admst:apply-templates select="." match="v2c:converttype"/>
      <admst:text format=" %(name)"/>
    </admst:join>);
\t\t
    <admst:apply-templates select="." match="v2c:converttype"/> d_$(function) (
    <admst:join select="variable[input='yes']" separator=", ">
      <admst:apply-templates select="." match="v2c:converttype"/><admst:text format=" %(name)"/>
    </admst:join>, <admst:join select="variable[input='yes']" separator=", ">
      <admst:apply-templates select="." match="v2c:converttype"/><admst:text format=" d_%(name)"/>
    </admst:join>);
  </admst:for-each>
</admst:template>

<admst:template match="analogfunctionC">
  <admst:for-each select="/module/analogfunction">
    <admst:variable name="globalanalogfunction" select="%(.)"/>
    <admst:variable name="function" select="%(name)"/>
    <admst:apply-templates select="." match="v2c:converttype"/>
$(DEV_NAME)::$function (
    <admst:join select="variable[input='yes']" separator=", ">
      <admst:apply-templates select="." match="v2c:converttype"/><admst:text format=" %(name)"/>
    </admst:join>)
{
    <admst:text format="\tdouble $function; "/>
    <admst:for-each select="variable[input='no' and output='no']">
      <admst:apply-templates select="." match="v2c:converttype"/>
      <admst:text format=" %(name);\n"/>
    </admst:for-each>
    <admst:apply-templates select="tree" match="af:print">
      <admst:text format="%(returned('x')/value)"/>
    </admst:apply-templates>
return $function;
}
double $(DEV_NAME)::d_$(function) (
    <admst:join select="variable[input='yes']" separator=", ">

    <admst:apply-templates select="." match="v2c:converttype"/>
    <admst:text format=" %(name)"/>
    </admst:join>,
    <admst:join select="variable[input='yes']" separator=", ">
      <admst:apply-templates select="." match="v2c:converttype"/>
      <admst:text format=" d_%(name)"/>
    </admst:join>)
{
    <admst:text format="double $function"/>
    <admst:for-each select="$globalanalogfunction/variable[input='yes']">
      <admst:variable name="ddx" select="%(name)"/>
      <admst:text format="; double $(function)_$(ddx)"/>
    </admst:for-each>
    <admst:for-each select="variable[input='no' and output='no']">
      <admst:variable name="name" select="%(name)"/>
      <admst:text format="; "/>
      <admst:apply-templates select="." match="v2c:converttype"/>
      <admst:text format=" $(name)"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:variable name="ddx" select="%(name)"/>
        <admst:text format="; "/>
        <admst:apply-templates select="." match="v2c:converttype"/><admst:text format=" $(name)_$(ddx)"/>
      </admst:for-each>
    </admst:for-each>;
    <admst:apply-templates select="tree" match="af:print:derivate">
    <admst:text format="%(returned('x')/value)"/>
    </admst:apply-templates>
return 
    <admst:join select="$globalanalogfunction/variable[input='yes']" separator="+">
      <admst:text format="$(function)_%(name)*d_%(name)"/>
    </admst:join>;
}

  </admst:for-each>
</admst:template>

<admst:template match="af:print:derivate">
  <admst:message test="[/dbg_xml='yes']" format="*(af)af:print:derivate*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='callfunction']">
      <admst:choose>
        <admst:when test="function[name='\$strobe']">
          <admst:variable name="outputfile" select="stdout"/>
        </admst:when>
      </admst:choose>
      <admst:variable name="args" select=""/>
      <admst:for-each select="function/arguments">
        <admst:apply-templates select="." match="af:print:expression">
          <admst:variable name="index" select="%(index(../arguments,.))"/>
          <admst:variable name="args" select="$args,%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="fprintf($outputfile$args); fprintf($outputfile,&quot;\\n&quot;);\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="af:print:expression">
        <admst:variable name="whileblock" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:apply-templates select="while" match="af:print:derivate">
        <admst:variable name="while" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:return name="x" string="while($whileblock)\n$while"/>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="if" match="af:print:expression">
        <admst:variable name="if" select="%(returned('x')/value)"/>
      </admst:apply-templates>
      <admst:apply-templates select="then" match="af:print:derivate">
        <admst:variable name="then" select="%(returned('x')/value)"/>
        </admst:apply-templates>
      <admst:if test="else">
        <admst:apply-templates select="else" match="af:print:derivate">
          <admst:variable name="then" select="$(then)else\n%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:if>
      <admst:return name="x" string="if($if)\n$then"/>
    </admst:when>
    <admst:when test="adms[datatypename='case']">
      <admst:apply-templates select="case" match="af:print:expression">
        <admst:variable name="case" select="switch ((int)%(returned('x')/value)) {\n"/>
      </admst:apply-templates>
      <admst:for-each select="caseitem">
        <admst:variable name="condition" select=""/>
        <admst:for-each select="condition">
          <admst:variable name="condition" select="$condition case %(.):"/>
        </admst:for-each>
        <admst:variable name="case" select="$case $condition"/>
        <admst:if test="[defaultcase='yes']">
          <admst:variable name="case" select="$case default:"/>
        </admst:if>
        <admst:variable name="case" select="$case \n"/>
        <admst:apply-templates select="code" match="af:print:derivate">
          <admst:variable name="case" select="$case%(returned('x')/value) break;\n"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="$case }"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']">
      <admst:fatal format="contribution not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='assignment']">
      <admst:variable name="lhs" select="%(lhs/name)"/>
      <admst:apply-templates select="rhs" match="af:print:expression">
        <admst:variable name="rhs" select=""/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:variable name="rhs" select="$rhs$(lhs)_%(name)=%(returned('dx.%(name)')/value);\n"/>
        </admst:for-each>
        <admst:variable name="rhs" select="$rhs$lhs=%(returned('x')/value);\n"/>
      </admst:apply-templates>
      <admst:return name="x" string="{$rhs}\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='nilled']">
      <admst:return name="x" string=";"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:variable name="block" select=""/>
      <admst:for-each select="item">
        <admst:apply-templates select="." match="af:print:derivate">
          <admst:variable name="block" select="$block%(returned('x')/value)"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" string="{$block}"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="%(datatypename): not handled inside blocks\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- 0 -->
</admst>

<admst version="2.3.0">
<admst:message test="[/dbg_xml='yes']" format="**enter scope2 (guess)**\n"/>

<!-- *********************************************************************************** -->
<!-- *********************** tree traverse - guesstopology ***************************** -->
<!-- *********************************************************************************** -->

<admst:for-each select="/module[extern='no']">
  <admst:value-of select="attribute[name='ngspicename']/value"/>
  <admst:variable name="module" select="%s"/>
  <admst:variable name="DEV_NAME" select="DEV_%(upper-case(name))"/>
  <admst:variable name="MODEL_NAME" select="MODEL_%(upper-case(name))"/>

  <admst:open file="$(filename)_top.hxx">
    <admst:apply-templates select="analog/code[datatypename='block']/variable"
      match="generate_c:declare_variable"/>

    <admst:for-each select="analog/code[datatypename='block']/item">
      <admst:if test="[(datatypename!='block') or (datatypename='block'
              and name!='initial_model'
              and name!='initializeModel'
              and name!='initial_instance'
              and name!='initializeInstance')]">
        <admst:apply-templates select="." match="ngspiseMODULEguesstopology.cevaluatetopology"/>
      </admst:if>
    </admst:for-each>
  </admst:open>
</admst:for-each>

<!--admst:for-each select="analog/code[datatypename='block']/item">
<admst:if test="[(datatypename!='block') or (datatypename='block'
and name!='initial_model' and name!='initializeModel' and name!='initial_instance'
and name!='initializeInstance')]">
<admst:apply-templates select="." match="evaluate.localvariables"/>
</admst:if>
</admst:for-each>
<admst:for-each select="$localvariables">
<admst:apply-templates select="" match="generate_c:declare_variable"/>
</admst:for-each-->

<!-- ----------------------------------------------------------------- -->
<admst:template match="ngspiseMODULEguesstopology.cevaluatetopology">
  <admst:message test="[/dbg_xml='yes']"
     format="*(guess)ngspiseMODULEguesstopology.cevaluatetopology*\n"/>

  <admst:if test="[datatypename='block']">
	{ // adms-block

    <admst:apply-templates select="variable" match="generate_c:declare_variable"/>
  </admst:if>

  <admst:choose>
    <admst:when test="[datatypename='callfunction']"/>
    <admst:when test="[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
      <admst:apply-templates select="while" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
    </admst:when>
    <admst:when test="[datatypename='conditional']">
      <admst:if test="if[nilled(variable[OPdependent='yes'])]">
        <admst:apply-templates select="if" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
        <admst:choose>
          <admst:when test="if/math[dependency='constant']">
            <admst:apply-templates select="if" match="expression:stringify:noprobe"/>
            <admst:text format="if\n(%s)\n"/>
            <admst:text format="{\n"/>
              <!--admst:apply-templates select="then/variable" match="generate_c:declare_variable"/-->
              <admst:apply-templates select="then" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
            <admst:text format="}\n"/>
            <admst:if test="[exists(else)]">
              <admst:text format="else\n"/>
              <admst:text format="{\n"/>
              <admst:apply-templates select="else" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
              <admst:text format="}\n"/>
            </admst:if>
          </admst:when>
          <admst:otherwise>
            <admst:apply-templates select="then" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
            <admst:apply-templates select="else" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
          </admst:otherwise>
        </admst:choose>
      </admst:if>
    </admst:when>
    <admst:when test="[datatypename='contribution']">
      <admst:if test="lhs[discipline/potential=nature]">
        <admst:choose>
          <admst:when test="lhs/branch[grounded='no']">
            <admst:text test="lhs/branch/nnode[location='internal']"
              format="_nodes[$node_prefix%(lhs/branch/nnode/name)]
                      = _nodes[$node_prefix%(lhs/branch/pnode/name)]; // nnode collapsed \n"/>

            <admst:text test="lhs/branch/pnode[location='internal']"
              format="_nodes[$node_prefix%(lhs/branch/pnode/name)] = _nodes[$node_prefix%(lhs/branch/nnode/name)]; // pnode collapsed \n"/>

          </admst:when>
          <admst:otherwise>
            <admst:text format="$node_prefix%(lhs/branch/pnode/name) = 0; /* pnode collapsed to GND */ \n"/>
          </admst:otherwise>
        </admst:choose>
      </admst:if>
      <admst:variable name="contribution" select="%(.)"/>
      <admst:variable name="psource" select="%(lhs/branch/pnode)"/>
      <admst:variable name="nsource" select="%(lhs/branch/nnode)"/>
      <admst:for-each select="rhs/probe">
        <admst:variable name="pprobe" select="%(branch/pnode)"/>
        <admst:variable name="nprobe" select="%(branch/nnode)"/>
        <admst:choose>
          <admst:when test="$contribution[static='yes']"> <admst:text format="  static_"/> </admst:when>
          <admst:when test="$contribution[dynamic='yes']"> <admst:text format="  dynamic_"/> </admst:when>
          <admst:when test="$contribution[whitenoise='yes']"> <admst:text format="  whitenoise_"/> </admst:when>
          <admst:when test="$contribution[flickernoise='yes']"> <admst:text format="  flickernoise_"/> </admst:when>
        </admst:choose>
        <admst:choose>
          <admst:when test="[($nprobe/grounded='no')and($nsource/grounded='no')]">
            <admst:text format="jacobian4(%($psource/name),%($nsource/name),%($pprobe/name),%($nprobe/name))\n"/>
          </admst:when>
          <admst:when test="[($nprobe/grounded='no')and($nsource/grounded='yes')]">
            <admst:text format="jacobian2p(%($psource/name),%($pprobe/name),%($nprobe/name))\n"/>
          </admst:when>
          <admst:when test="[$nsource/grounded='no']">
            <admst:text format="jacobian2s(%($psource/name),%($nsource/name),%($pprobe/name))\n"/>
          </admst:when>
          <admst:when test="[$nsource/grounded='yes']">
            <admst:text format="jacobian1(%($psource/name),%($pprobe/name))\n"/>
          </admst:when>
        </admst:choose>
      </admst:for-each>
    </admst:when>
    <admst:when test="[datatypename='assignment']">
      <!-- ****** TODO: "fixed" here ********** -->
      <!--admst:if test="[(lhs/insource='yes') and (lhs/OPdependent='no')]">
        <admst:apply-templates select="lhs" match="variable:lhs"/>
          <admst:text format="/* 4998 */ ="/>
        <admst:apply-templates select="rhs" match="expression:stringify:noprobe"/>
        <admst:text format="%s;\n"/>
      </admst:if-->
    </admst:when>
    <admst:when test="[datatypename='block']">
      <admst:apply-templates select="item" match="ngspiseMODULEguesstopology.cevaluatetopology"/>
    </admst:when>
    <admst:when test="[datatypename='expression']"/>
    <admst:when test="[datatypename='probe']"/>
    <admst:when test="[datatypename='variable']"/>
    <admst:when test="[datatypename='mapply_unary']"/>
    <admst:when test="[datatypename='mapply_binary']"/>
    <admst:when test="[datatypename='mapply_ternary']"/>
    <admst:when test="[datatypename='function']"/>
    <admst:when test="[datatypename='number']"/>
    <admst:when test="[datatypename='string']"/>
    <admst:when test="[datatypename='nilled']"/>
    <admst:when test="[datatypename='blockvariable']"/>
    <admst:otherwise>
      <admst:fatal format="%(datatypename): adms element not implemented\n"/>
    </admst:otherwise>
  </admst:choose>

  <admst:if test="[datatypename='block']">
	} // adms-block

  </admst:if>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="function">
<admst:message test="[/dbg_xml='yes']" format="*(guess)function*\n"/>
  <admst:variable name="function" select="%(name)"/>
  <admst:variable name="args" select=""/>
  <admst:for-each select="arguments">
    <admst:variable test="[$args='']" name="args" select="$args,"/>
    <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
    <admst:variable name="args" select="$args%s"/>
  </admst:for-each>
  <admst:choose>
    <admst:when test="[name='\$simparam']">
      <admst:apply-templates select="." match="function:simparam"/>
    </admst:when>
    <admst:when test="[name='analysis']">
      <admst:apply-templates select="." match="function:analysis"/>
    </admst:when>
    <admst:when test="[name='\$given']">
      <admst:for-each select="arguments[position(.)=1]">
        <admst:if test="[datatypename!='variable']">
          <admst:error format="$given: argument is not a variable\n"/>
        </admst:if>
        <admst:if test="[input='no']">
          <admst:value-of select="name"/>
          <admst:error format="$given(%s): argument is not a parameter\n"/>
        </admst:if>
        <admst:choose>
          <admst:when test="[parametertype='model']">
            <admst:value-of select="name"/>
            <admst:value-to select="/simulator/tmp" value="m->$var_prefix)%s.has_good_value())"/>
          </admst:when>
          <admst:when test="[parametertype='instance']">
            <admst:value-of select="name"/>
            <admst:value-to select="/simulator/tmp" value="$(var_prefix)%s.has_good_value()"/>
          </admst:when>
          <admst:otherwise>
            <admst:error format="$given(%s): should not be reached\n"/>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="mycode" select=""/>
      <admst:if test="[exists(arguments)]">
        <admst:for-each select="arguments">
          <admst:apply-templates select="." match="subexpression:stringify:noprobe"/>
          <admst:choose>
            <admst:when test="[$mycode='']">
              <admst:variable name="mycode" select="%s"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="mycode" select="$mycode,%s"/>
            </admst:otherwise>
          </admst:choose>
        </admst:for-each>
        <admst:variable name="mycode" select="($mycode)"/>
      </admst:if>
      <admst:variable name="mycode" select="%(fname(.)/[name='fname']/value)$mycode"/>
      <admst:value-to select="/simulator/tmp" value="$mycode"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<admst:template match="evaluate.localvariables">
<admst:message test="[/dbg_xml='yes']" format="*(guess)evaluate.localvariables*\n"/>
  <admst:choose>
    <admst:when test="[datatypename='assignment']">
      <admst:if test="[(lhs/insource='yes') and (lhs/OPdependent='no')]">
        <admst:push select="lhs[scope='local']" into="$localvariables" onduplicate="ignore"/>
      </admst:if>
    </admst:when>
    <admst:when test="[datatypename='block']">
      <admst:apply-templates select="item" match="evaluate.localvariables"/>
    </admst:when>
    <admst:when test="[datatypename='conditional']">
      <admst:push select="if/variable[scope='local' and OPdependent='no']"
        into="$localvariables" onduplicate="ignore"/>
      <admst:apply-templates select="then" match="evaluate.localvariables"/>
      <admst:apply-templates select="else" match="evaluate.localvariables"/>
    </admst:when>
    <admst:when test="[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="evaluate.localvariables"/>
    </admst:when>
    <admst:when test="[datatypename='case']">
      <admst:apply-templates select="caseitem/code" match="evaluate.localvariables"/>
    </admst:when>
    <admst:when test="[datatypename='contribution']"/>
    <admst:when test="[datatypename='nilled']"/>
    <admst:when test="[datatypename='callfunction']"/>
    <admst:when test="[datatypename='blockvariable']"/>
    <admst:otherwise>
      <admst:error format="'%(datatypename): should not be reached\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!-- -------------------------same as function:getname?-------------------------------- -->
<admst:template match="fname">
  <admst:message test="[/dbg_xml='yes']" format="*(guess)fname*\n"/>
  <admst:choose>
    <admst:when test="[name='div']"><admst:return name="fname" value="_div1"/></admst:when>
    <admst:when test="[name='abs']"><admst:return name="fname" value="fabs"/></admst:when>
    <admst:when test="[name='\$shrinkl']"><admst:return name="fname" value="shrinkl"/></admst:when>
    <admst:when test="[name='\$shrinka']"><admst:return name="fname" value="shrinka"/></admst:when>
    <admst:when test="[name='log']"><admst:return name="fname" value="log10"/></admst:when>
    <admst:when test="[name='ln']"><admst:return name="fname" value="log"/></admst:when>
    <admst:when test="[name='limexp']"><admst:return name="fname" value="limexp"/></admst:when>
    <admst:when test="[name='\$limexp']"><admst:return name="fname" value="limexp"/></admst:when>
    <admst:when test="[name='\$model']"><admst:return name="fname" value="_modelname"/></admst:when>
    <admst:when test="[name='\$instance']"><admst:return name="fname" value="_instancename"/></admst:when>
    <admst:when test="[name='\$temperature']">
      <admst:return name="fname" value="_ambient_temp"/>
    </admst:when>
    <admst:otherwise><admst:return name="fname" value="%(name)"/></admst:otherwise>
  </admst:choose>
</admst:template>

<!-- 1 -->
</admst>

<admst version="2.3.0">
<admst:message test="[/dbg_xml='yes']" format="**enter scope3 (tr_eval)**\n"/>

<!-- *********************************************************************************** -->
<!-- *********************************** ngspiceMODULE.hxx.xml ************************* -->
<!-- *********************************************************************************** -->

<!-- ----------------------------------------------------------------- -->
<!-- ------------------ create hxx for tr ---------------------------- -->
<admst:variable name="eval_mode" select="tr"/>
<admst:for-each select="/module[extern='no']">
  <admst:value-of select="attribute[name='ngspicename']/value"/>
  <admst:variable name="module" select="%s"/>
  <admst:variable name="DEV_NAME" select="DEV_%(upper-case(name))"/>
  <admst:variable name="MODEL_NAME" select="MODEL_%(upper-case(name))"/>

  <admst:open file="$(filename)_tr.hxx">
    <admst:apply-templates select="analog/code" match="analog:evaluate" required="yes"/>
  </admst:open>
</admst:for-each>

<!-- ----------------------------------------------------------------- -->
<!-- ------------------ create hxx for ac ---------------------------- -->
<admst:variable name="eval_mode" select="ac"/>
<admst:for-each select="/module[extern='no']">
  <admst:value-of select="attribute[name='ngspicename']/value"/>
  <admst:variable name="module" select="%s"/>
  <admst:variable name="DEV_NAME" select="DEV_%(upper-case(name))"/>
  <admst:variable name="MODEL_NAME" select="MODEL_%(upper-case(name))"/>

  <admst:open file="$(filename)_ac.hxx">
    <admst:apply-templates select="analog/code" match="analog:evaluate" required="yes"/>
  </admst:open>
</admst:for-each>

<!-- ----------------------------------------------------------------- -->
<!-- -------------------- analog//contribution ----------------------- -->
<!-- never called ?? admst:template match="contribution:ddt">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)contribution:ddt*\n"/>
  <admst:choose>
    <admst:when test="[whitenoise='yes' or flickernoise='yes']">
        <admst:warning format=" not implemented.\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:apply-templates select="." match="contribution:nonoise:ddt"/>
    </admst:otherwise>
  </admst:choose>
</admst:template-->

<!-- ----------------------------------------------------------------- -->
<admst:template match="contribution:nonoise:ddt">
  <admst:error format="not called?\n"/>
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)contribution:nonoise:ddt*\n"/>

  <admst:if test="rhs[not(nilled(function[class='builtin']))]">
  <admst:text format="{ // 5206 \n"/>
    <admst:for-each select="rhs/function">
      <admst:value-of select="position(.)-1"/>
      <admst:apply-templates select="." match="function:getCid"/>
      <admst:text format="double __%s_%s=0.0;\n"/>
    </admst:for-each>
    <admst:apply-templates select="." match="ddx:function:computation"/>
  </admst:if>
  <admst:choose>
  </admst:choose>
  <admst:choose>
    <admst:when test="lhs/branch/nnode[grounded='no']">
      <admst:text format="//exec %e(rhs/tree))\n"/>
      <admst:text format=" $e\n"/>
    </admst:when>
    <admst:otherwise>
      <admst:value-of select="lhs/branch/pnode/name"/>
      <admst:text format="// exec residual1(%s,%e(rhs/tree))\n"/>
      <admst:text format="//founnde $e\n"/>
    </admst:otherwise>
  </admst:choose>
  <admst:for-each select="rhs/probe">
    <admst:text format="// Rhsp %(.)\n"/>

    <admst:variable name="probepnode" select="%(branch/pnode)"/>
    <admst:variable name="probennode" select="%(branch/nnode)"/>
    <admst:variable name="probepnodename" select="%($probepnode/name)"/>
    <admst:variable name="probennodename" select="%($probennode/name)"/>
    <admst:variable name="pprobe" select="%(.)"/>
    <admst:apply-templates select="../tree" match="%(adms/datatypename)"/>
    <admst:choose>
    </admst:choose>
  </admst:for-each>
</admst:template>

<!-- ----------------------------------------------------------------- -->
<!-- -------------------- analog//contribution ------------------------ -->
<admst:template match="contribution">
// MODULE contribution %(.)
  <admst:choose>
    <admst:when test="[whitenoise='yes' or flickernoise='yes']">
      <admst:variable name="SkipFVariable" select="y"/>
        <admst:variable name="dependency" select="%(math/dependency)"/>
        <admst:if test="[whitenoise='yes']">
          <admst:if test="[$dependency='constant']">
            <admst:text format="tnoise%(index($tnoise/item,.))="/>
          </admst:if>
          <admst:if test="[$dependency!='constant']">
            <admst:text format="wnoise%(index($wnoise/item,.))="/>
          </admst:if>
          <admst:text format=" // things 2...  \n"/>
          <admst:apply-templates select="rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
          <admst:text format="$e;\n"/>
        </admst:if>
        <admst:if test="[flickernoise='yes']">
          <admst:text format="fpnoise%(index($fnoise/item,.))="/>
          <admst:apply-templates select="rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
          <admst:text format="$e;\n"/>
          <admst:text format="fenoise%(index($fnoise/item,.))="/>
          <admst:apply-templates select="rhs/tree/arguments[2]" match="%(adms/datatypename)"/>
          <admst:text format="$e;\n"/>
        </admst:if>
      <admst:variable name="SkipFVariable" select="n"/>
    </admst:when>
    <admst:otherwise>
// nonoise contribution: %(.)

      <admst:apply-templates select="." match="contribution:nonoise"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- compute node arguments of noise routines -->
<admst:template match="noisebranch">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)noisebranch*\n"/>
  <admst:variable name="n1" select=""/>
  <admst:choose>
    <admst:when test="lhs[grounded='yes']">
      <admst:return name="noisebranch" value="%(lhs/branch/pnode/name),%(lhs/branch/pnode/name)"/>
    </admst:when>
    <admst:otherwise>
      <admst:return name="noisebranch" value="%(lhs/branch/pnode/name),%(lhs/branch/nnode/name)"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- compute range of variables -->
<admst:template match="variable:range:foreach">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:range:foreach*\n"/>
  <admst:choose>
    <admst:when test="infexpr[hasspecialnumber='NO']">
      <admst:apply-templates select="infexpr/tree" match="%(adms/datatypename)"/>
      <admst:variable name="lower" select="$e"/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="lower" select="-inf"/>
    </admst:otherwise>
  </admst:choose>
  <admst:choose>
    <admst:when test="supexpr[hasspecialnumber='NO']">
      <admst:apply-templates select="supexpr/tree" match="%(adms/datatypename)"/>
      <admst:variable name="upper" select="$e"/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="upper" select="inf"/>
    </admst:otherwise>
  </admst:choose>
  <admst:choose>
    <admst:when test="[type='include']">
      <admst:variable name="rangetype" select="from"/>
    </admst:when>
    <admst:when test="[type='exclude']">
      <admst:variable name="rangetype" select="exclude"/>
    </admst:when>
  </admst:choose>
  <admst:variable name="interval" select="$lower:$upper"/>
  <admst:choose>
    <admst:when test="[infboundtype='range_bound_include' and supboundtype='range_bound_include']">
      <admst:variable name="interval" select="[$interval]"/>
    </admst:when>
    <admst:when test="[infboundtype='range_bound_include' and supboundtype='range_bound_exclude']">
      <admst:variable name="interval" select="[$interval)"/>
    </admst:when>
    <admst:when test="[infboundtype='range_bound_exclude' and supboundtype='range_bound_include']">
      <admst:variable name="interval" select="($interval]"/>
    </admst:when>
    <admst:when test="[infboundtype='range_bound_exclude' and supboundtype='range_bound_exclude']">
      <admst:variable name="interval" select="($interval)"/>
    </admst:when>
  </admst:choose>
  <admst:text format="$rangetype $interval "/>
</admst:template>
<admst:template match="variable:range">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:range*\n"/>
  <admst:choose>
    <admst:when test="range">
      <admst:text format="        &quot;"/>
      <admst:apply-templates select="range" match="variable:range:foreach"/>
      <admst:text format="&quot;"/>
    </admst:when>
    <admst:otherwise>
      <admst:text format="        NULL"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- expression//function: mapping verilog-name == C-name of function -->
<!-- ----------------------------------------------------------------- -->
<admst:template match="af:print:expression">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)af:print:expression*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='expression']">
      <admst:apply-templates select="tree" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="expression" select="%s"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:value-of select="name"/>
          <admst:value-of select="returned('dx.%s')/value"/>
          <admst:value-of select="name"/>
          <admst:variable name="dx_%s" select="%s"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:return name="x" value="$expression"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:value-of select="name"/>
        <admst:value-of select="name"/>
        <admst:return name="dx.%s" value="$(dx_%s)"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='probe']">
      <admst:fatal format="probe not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='variable']">
      <admst:value-of select="name"/>
      <admst:variable name="variable" select="%s"/>
      <admst:return name="x" value="$variable"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:value-of select="name"/>
        <admst:variable name="ddx" select="%s"/>
        <admst:choose>
          <admst:when test="[$variable='$ddx']">
            <admst:return name="dx.$ddx" value="1.0"/>
          </admst:when>
          <admst:when test="../..[input='yes']">
            <admst:return name="dx.$ddx" value="0.0"/>
          </admst:when>
          <admst:otherwise>
            <admst:return name="dx.$ddx" value="$(variable)_$ddx"/>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='mapply_unary']">
      <admst:if test="[name='plus']">
        <admst:variable name="op" select="+"/>
      </admst:if>
      <admst:if test="[name='minus']">
        <admst:variable name="op" select="-"/>
      </admst:if>
      <admst:if test="[name='not']">
        <admst:variable name="op" select="!"/>
      </admst:if>
      <admst:if test="[name='bw_not']">
        <admst:variable name="op" select="~"/>
      </admst:if>
      <admst:apply-templates select="arg1" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/><admst:variable name="arg1" select="%s"/>
      </admst:apply-templates>
      <admst:return name="x" value="($op$arg1)"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:value-of select="name"/>
        <admst:return name="dx.%s" value="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='mapply_binary']">
      <admst:apply-templates select="arg1" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="x" select="%s"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:value-of select="name"/>
          <admst:value-of select="returned('dx.%s')/value"/>
          <admst:value-of select="name"/>
          <admst:variable name="dx_%s" select="%s"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:apply-templates select="arg2" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="y" select="%s"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:value-of select="name"/>
          <admst:value-of select="returned('dx.%s')/value"/>
          <admst:value-of select="name"/>
          <admst:variable name="dy_%s" select="%s"/>
        </admst:for-each>
      </admst:apply-templates>
      <admst:choose>
        <admst:when test="[name='addp']">
          <admst:choose>
            <admst:when test="[$x='0.0' and $y='0.0']">
              <admst:return name="x" value="0.0"/>
            </admst:when>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" value="(+$y)"/>
            </admst:when>
            <admst:when test="[$y='0.0']">
              <admst:return name="x" value="$x"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" value="($x+$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/><admst:variable name="df" select="%s"/>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:when>
              <admst:when test="[$y='0.0']">
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" value="(+$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" value="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" value="($dx+$dy)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='addm']">
          <admst:choose>
            <admst:when test="[$x='0.0' and $y='0.0']">
              <admst:return name="x" value="0.0"/>
            </admst:when>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" value="(-$y)"/>
            </admst:when>
            <admst:when test="[$y='0.0']">
              <admst:return name="x" value="$x"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" value="($x-$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/><admst:variable name="df" select="%s"/>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:when>
              <admst:when test="[$y='0.0']">
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" value="(-$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" value="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" value="($dx-$dy)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='multtime']">
          <admst:choose>
            <admst:when test="[$x='0.0' or $y='0.0']">
              <admst:return name="x" value="0.0"/>
            </admst:when>
            <admst:when test="[$x='1.0' and $y='1.0']">
              <admst:return name="x" value="1.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" value="($x*$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/><admst:variable name="df" select="%s"/>
            <admst:choose>
              <admst:when test="[$x='0.0' or $y='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0' and $y='1.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$df)"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$x='0.0' and $y='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='1.0']">
                <admst:return name="dx.$df" value="($x)"/>
              </admst:when>
              <admst:when test="[$dx='1.0' and $dy='0.0']">
                <admst:return name="dx.$df" value="($y)"/>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:return name="dx.$df" value="($x*$dy)"/>
              </admst:when>
              <admst:when test="[$dy='0.0']">
                <admst:return name="dx.$df" value="$dx*$y"/>
              </admst:when>
              <admst:when test="[$dx='1.0' and $dy='1.0']">
                <admst:return name="dx.$df" value="($x+$y)"/>
              </admst:when>
              <admst:when test="[$dx='1.0']">
                <admst:return name="dx.$df" value="($y+($dy*$x))"/>
              </admst:when>
              <admst:when test="[$dy='1.0']">
                <admst:return name="dx.$df" value="($dx*$y)+$x"/>
              </admst:when>
              <admst:when test="[$x='1.0']">
                <admst:return name="dx.$df" value="$dy"/>
              </admst:when>
              <admst:when test="[$y='1.0']">
                <admst:return name="dx.$df" value="$dx"/>
              </admst:when>
              <admst:otherwise>
                <admst:return name="dx.$df" value="(($dx*$y)+($x*$dy))"/>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:when test="[name='multdiv']">
          <admst:choose>
            <admst:when test="[$x='0.0']">
              <admst:return name="x" value="0.0"/>
            </admst:when>
            <admst:when test="[$x='1.0' and $y='1.0']">
              <admst:return name="x" value="1.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:return name="x" value="($x/$y)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/><admst:variable name="df" select="%s"/>
            <admst:choose>
              <admst:when test="[$x='0.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0' and $y='1.0']">
                <admst:variable name="dx" select="0.0"/>
                <admst:variable name="dy" select="0.0"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="dx" select="$(dx_$df)"/>
                <admst:variable name="dy" select="$(dy_$(df))"/>
              </admst:otherwise>
            </admst:choose>
            <admst:choose>
              <admst:when test="[$x='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$dx='0.0' and $dy='0.0']">
                <admst:return name="dx.$df" value="0.0"/>
              </admst:when>
              <admst:when test="[$x='1.0']">
                <admst:choose>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" value="(-1/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" value="(-$dy/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:when test="[$dx='0.0']">
                <admst:choose>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" value="(-$x/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" value="(-($x*$dy)/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:when test="[$dx='1.0']">
                <admst:choose>
                  <admst:when test="[$dy='0.0']">
                    <admst:return name="dx.$df" value="(1/$y)"/>
                  </admst:when>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" value="(($y-$x)/($y*$y))"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" value="(($y-($x*$dy))/($y*$y))"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:when>
              <admst:otherwise>
                <admst:choose>
                  <admst:when test="[$y='1.0']">
                    <admst:return name="dx.$df" value="$dx"/>
                  </admst:when>
                  <admst:when test="[$dy='0.0']">
                    <admst:return name="dx.$df" value="$dx/$y"/>
                  </admst:when>
                  <admst:when test="[$dy='1.0']">
                    <admst:return name="dx.$df" value="(($dx*$y)-$x)/($y*$y)"/>
                  </admst:when>
                  <admst:otherwise>
                    <admst:return name="dx.$df" value="($dx*$y-$x*$dy)/($y*$y)"/>
                  </admst:otherwise>
                </admst:choose>
              </admst:otherwise>
            </admst:choose>
          </admst:for-each>
        </admst:when>
        <admst:otherwise>
          <admst:choose>
            <admst:when test="[name='bw_equr']">
              <admst:return name="x" value="($x^~$y)"/>
            </admst:when>
            <admst:when test="[name='bw_equl']">
              <admst:return name="x" value="($x~^$y)"/>
            </admst:when>
            <admst:when test="[name='bw_xor']">
              <admst:return name="x" value="($x^$y)"/>
            </admst:when>
            <admst:when test="[name='bw_or']">
              <admst:return name="x" value="($x|$y)"/>
            </admst:when>
            <admst:when test="[name='bw_and']">
              <admst:return name="x" value="($x&amp;$y)"/>
            </admst:when>
            <admst:when test="[name='or']">
              <admst:return name="x" value="($x||$y)"/>
            </admst:when>
            <admst:when test="[name='and']">
              <admst:return name="x" value="($x&amp;&amp;$y)"/>
            </admst:when>
            <admst:when test="[name='equ']">
              <admst:return name="x" value="($x==$y)"/>
            </admst:when>
            <admst:when test="[name='multmod']">
              <admst:return name="x" value="((int)$x%%(int)$y)"/>
            </admst:when>
            <admst:when test="[name='notequ']">
              <admst:return name="x" value="($x!=$y)"/>
            </admst:when>
            <admst:when test="[name='lt']">
              <admst:return name="x" value="($x&lt;$y)"/>
            </admst:when>
            <admst:when test="[name='lt_equ']">
              <admst:return name="x" value="($x&lt;=$y)"/>
            </admst:when>
            <admst:when test="[name='gt']">
              <admst:return name="x" value="($x&gt;$y)"/>
            </admst:when>
            <admst:when test="[name='gt_equ']">
              <admst:return name="x" value="($x&gt;=$y)"/>
            </admst:when>
            <admst:when test="[name='shiftr']">
              <admst:return name="x" value="($x&gt;&gt;$y)"/>
            </admst:when>
            <admst:when test="[name='shiftl']">
              <admst:return name="x" value="($x&lt;&lt;$y)"/>
            </admst:when>
            <admst:otherwise>
              <admst:value-of select="name"/>
              <admst:error format="%s: function not handled\n"/>
            </admst:otherwise>
          </admst:choose>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/>
            <admst:return name="dx.%s" value="0.0"/>
          </admst:for-each>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="adms[datatypename='mapply_ternary']">
      <admst:apply-templates select="arg1" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="x" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="arg2" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="y" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="arg3" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="z" select="%s"/>
      </admst:apply-templates>
      <admst:if test="[name='conditional']">
        <admst:return name="x" value="($x?$y:$z)"/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:value-of select="name"/>
          <admst:value-of select="name"/>
          <admst:value-of select="name"/>
          <admst:return name="dx.%s" value="($x?$dy_%s:$dz_%s)"/>
        </admst:for-each>
      </admst:if>
    </admst:when>
    <admst:when test="adms[datatypename='function']">
      <admst:apply-templates select="." match="function:getname">
        <admst:value-of select="returned('function:getname')/value"/>
        <admst:variable name="function" select="%s"/>
      </admst:apply-templates>
      <admst:variable name="args" select=""/>
      <admst:for-each select="arguments">
        <admst:if test="[not($args='')]">
          <admst:variable name="args" select="$args,"/>
        </admst:if>
        <admst:apply-templates select="." match="af:print:expression">
          <admst:value-of select="index(../arguments,.)"/>
          <admst:variable name="index" select="%s"/>
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="args" select="$args%s"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/>
            <admst:value-of select="returned('dx.%s')/value"/>
            <admst:value-of select="name"/>
            <admst:variable name="arg$(index)_%s" select="%s"/>
          </admst:for-each>
        </admst:apply-templates>
      </admst:for-each>
      <admst:choose>
        <admst:when test="[ name='cos' or name='sin' or name='tan' or name='cosh'
                            or name='sinh' or name='tanh' or name='acos' or name='asin'
                            or name='atan' or name='ln' or name='log' or name='exp'
                            or name='sqrt' or name='abs' or name='limexp'
                            or name='div' or name='pow' or name='hypot' or name='min' or name='max' ]">
          <admst:return name="x" value="_$function($args)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/> <admst:variable name="name" select="%s"/>
            <admst:variable name="ret" select=""/>
            <admst:for-each select="../../arguments">
              <admst:if test="[not($ret='')]">
                <admst:variable name="ret" select="$ret+"/>
              </admst:if>
              <admst:value-of select="index(../arguments,.)"/>
              <admst:variable name="index" select="%s"/>
              <admst:variable name="ret" select="$(ret)_d$(index)_$function($args)*($(arg$(index)_$name))"/>
            </admst:for-each>
            <admst:return name="dx.$name" value="$ret"/>
          </admst:for-each>
        </admst:when>
        <admst:otherwise>
          <admst:return name="x" value="$function($args)"/>
          <admst:for-each select="$globalanalogfunction/variable[input='yes']">
            <admst:value-of select="name"/> <admst:variable name="name" select="%s"/>
            <admst:variable name="darg" select=""/>
            <admst:for-each select="../../arguments">
              <admst:value-of select="index(../arguments,.)"/>
              <admst:variable name="index" select="%s"/>
              <admst:variable name="darg" select="$darg,($(arg$(index)_$name))"/>
            </admst:for-each>
            <admst:return name="dx.$name" value="d_$function($args$darg)"/>
          </admst:for-each>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="adms[datatypename='string']">
      <admst:value-of select="value"/>
      <admst:return name="x" value="&quot;%s&quot;"/>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:value-of select="name"/>
        <admst:return name="dx.%s" value="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:when test="adms[datatypename='number']">
      <admst:choose>
        <admst:when test="[scalingunit='1']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="%s"/>
        </admst:when>
        <admst:when test="[scalingunit='E']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+18)"/>
        </admst:when>
        <admst:when test="[scalingunit='P']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+15)"/>
        </admst:when>
        <admst:when test="[scalingunit='T']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+12)"/>
        </admst:when>
        <admst:when test="[scalingunit='G']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+9)"/>
        </admst:when>
        <admst:when test="[scalingunit='M']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+6)"/>
        </admst:when>
        <admst:when test="[scalingunit='k']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+3)"/>
        </admst:when>
        <admst:when test="[scalingunit='h']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+2)"/>
        </admst:when>
        <admst:when test="[scalingunit='D']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e+1)"/>
        </admst:when>
        <admst:when test="[scalingunit='d']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-1)"/>
        </admst:when>
        <admst:when test="[scalingunit='c']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-2)"/>
        </admst:when>
        <admst:when test="[scalingunit='m']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-3)"/>
        </admst:when>
        <admst:when test="[scalingunit='u']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-6)"/>
        </admst:when>
        <admst:when test="[scalingunit='n']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-9)"/>
        </admst:when>
        <admst:when test="[scalingunit='A']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-10)"/>
        </admst:when>
        <admst:when test="[scalingunit='p']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-12)"/>
        </admst:when>
        <admst:when test="[scalingunit='f']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-15)"/>
        </admst:when>
        <admst:when test="[scalingunit='a']">
          <admst:value-of select="value"/>
          <admst:return name="x" value="(%s*1.0e-18)"/>
        </admst:when>
        <admst:otherwise>
          <admst:value-of select="scalingunit"/>
          <admst:fatal format="scaling unit not supported: %s\n"/>
        </admst:otherwise>
      </admst:choose>
      <admst:for-each select="$globalanalogfunction/variable[input='yes']">
        <admst:value-of select="name"/>
        <admst:return name="dx.%s" value="0.0"/>
      </admst:for-each>
    </admst:when>
    <admst:otherwise>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- ------------------------------------------------------ -->
<admst:template match="af:print:derivate">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)af:print:derivate*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='callfunction']">
      <admst:choose>
        <admst:when test="function[name='\$strobe']">
          <admst:variable name="outputfile" select="stdout"/>
        </admst:when>
      </admst:choose>
      <admst:variable name="args" select=""/>
      <admst:for-each select="function/arguments">
        <admst:apply-templates select="." match="af:print:expression">
          <admst:value-of select="index(../arguments,.)"/>
          <admst:variable name="index" select="%s"/>
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="args" select="$args,%s"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="fprintf($outputfile$args); fprintf($outputfile,&quot;\n&quot;);\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="whileblock" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="while" match="af:print:derivate">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="while" select="%s"/>
      </admst:apply-templates>
      <admst:return name="x" value="while($whileblock)\n$while"/>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="if" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="if" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="then" match="af:print:derivate">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="then" select="%s"/>
        </admst:apply-templates>
      <admst:if test="else">
        <admst:apply-templates select="else" match="af:print:derivate">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="then" select="$(then)else\n%s"/>
        </admst:apply-templates>
      </admst:if>
      <admst:return name="x" value="if($if)\n$then"/>
    </admst:when>
    <admst:when test="adms[datatypename='case']">
      <admst:apply-templates select="case" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="case" select="switch ((int)%s) {\n"/>
      </admst:apply-templates>
      <admst:for-each select="caseitem">
        <admst:variable name="condition" select=""/>
        <admst:for-each select="condition">
          <admst:value-of select="."/>
          <admst:variable name="condition" select="$condition case %s:"/>
        </admst:for-each>
        <admst:variable name="case" select="$case $condition"/>
        <admst:if test="[defaultcase='yes']">
          <admst:variable name="case" select="$case default:"/>
        </admst:if>
        <admst:variable name="case" select="$case \n"/>
        <admst:apply-templates select="code" match="af:print:derivate">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="case" select="$case%s break;\n"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="$case }"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']">
      <admst:fatal format="contribution not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='assignment']">
      <admst:value-of select="lhs/name"/>
      <admst:variable name="lhs" select="%s"/>
      <admst:apply-templates select="rhs" match="af:print:expression">
        <admst:variable name="rhs" select=""/>
        <admst:for-each select="$globalanalogfunction/variable[input='yes']">
          <admst:value-of select="name"/>
          <admst:value-of select="returned('dx.%s')/value"/>
          <admst:value-of select="name"/>
          <admst:variable name="rhs" select="$rhs$(lhs)_%s=%s;\n"/>
        </admst:for-each>
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="rhs" select="$rhs$lhs=%s;\n"/>
      </admst:apply-templates>
      <admst:return name="x" value="{$rhs}\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='nilled']">
      <admst:return name="x" value=";\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:variable name="block" select=""/>
      <admst:for-each select="item">
        <admst:apply-templates select="." match="af:print:derivate">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="block" select="$block%s"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="{$block}"/>
    </admst:when>
    <admst:otherwise>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- ------------------------------------------------------------ -->
<admst:template match="af:print">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)af:print*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='callfunction']">
      <admst:choose>
        <admst:when test="function[name='\$strobe']">
          <admst:variable name="outputfile" select="stdout"/>
        </admst:when>
      </admst:choose>
      <admst:variable name="args" select=""/>
      <admst:for-each select="function/arguments">
        <admst:apply-templates select="." match="af:print:expression">
          <admst:value-of select="index(../arguments,.)"/>
          <admst:variable name="index" select="%s"/>
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="args" select="$args,%s"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="fprintf($outputfile$args); fprintf($outputfile,&quot;\n&quot;);\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="whileblock" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="while" match="af:print">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="while" select="%s"/>
      </admst:apply-templates>
      <admst:return name="x" value="while($whileblock)\n$while"/>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="if" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="if" select="%s"/>
      </admst:apply-templates>
      <admst:apply-templates select="then" match="af:print">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="then" select="%s"/>
        </admst:apply-templates>
      <admst:if test="else">
        <admst:apply-templates select="else" match="af:print">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="then" select="$(then)else\n%s"/>
        </admst:apply-templates>
      </admst:if>
      <admst:return name="x" value="if($if)\n$then"/>
    </admst:when>
    <admst:when test="adms[datatypename='case']">
      <admst:apply-templates select="case" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:variable name="case" select="switch ((int)%s) {\n"/>
      </admst:apply-templates>
      <admst:for-each select="caseitem">
        <admst:variable name="condition" select=""/>
        <admst:for-each select="condition">
          <admst:value-of select="."/>
          <admst:variable name="condition" select="$condition case %s:"/>
        </admst:for-each>
        <admst:variable name="case" select="$case $condition"/>
        <admst:if test="[defaultcase='yes']">
          <admst:variable name="case" select="$case default:"/>
        </admst:if>
        <admst:variable name="case" select="$case \n"/>
        <admst:apply-templates select="code" match="af:print">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="case" select="$case%s break;\n"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="$case }"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']">
      <admst:fatal format="contribution not allowed inside analog functions\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='assignment']">
      <admst:apply-templates select="rhs" match="af:print:expression">
        <admst:value-of select="returned('x')/value"/>
        <admst:value-of select="../lhs/name"/>
        <admst:return name="x" value="%s=%s;\n"/>
      </admst:apply-templates>
    </admst:when>
    <admst:when test="adms[datatypename='nilled']">
      <admst:return name="x" value=";\n"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:variable name="block" select=""/>
      <admst:for-each select="item">
        <admst:apply-templates select="." match="af:print">
          <admst:value-of select="returned('x')/value"/>
          <admst:variable name="block" select="$block%s"/>
        </admst:apply-templates>
      </admst:for-each>
      <admst:return name="x" value="{$block}"/>
    </admst:when>
    <admst:otherwise>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!--
* Returns the type of a variable. The returned type
* is either int, double, or char *. This template is
* used to create the mint:instance and model data structures.
-->
<admst:template match="vtype">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)vtype*\n"/>
  <admst:choose>
    <admst:when test="[type='integer']">int</admst:when>
    <admst:when test="[type='real']">double</admst:when>
    <admst:when test="[type='string']">char*</admst:when>
    <admst:otherwise>
      <admst:fatal format="variable type unknown\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- ----------------------------------------- -->
<admst:template match="bname">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)bname*\n"/>
  <admst:choose>
    <admst:when test="[name='bw_equr']">
      <admst:return name="bname" value="^~"/>
    </admst:when>
    <admst:when test="[name='bw_equl']">
      <admst:return name="bname" value="~^"/>
    </admst:when>
    <admst:when test="[name='bw_xor']">
      <admst:return name="bname" value="^"/>
    </admst:when>
    <admst:when test="[name='bw_or']">
      <admst:return name="bname" value="|"/>
    </admst:when>
    <admst:when test="[name='bw_and']">
      <admst:return name="bname" value="&amp;"/>
    </admst:when>
    <admst:when test="[name='or']">
      <admst:return name="bname" value="||"/>
    </admst:when>
    <admst:when test="[name='and']">
      <admst:return name="bname" value="&amp;&amp;"/>
    </admst:when>
    <admst:when test="[name='equ']">
      <admst:return name="bname" value="=="/>
    </admst:when>
    <admst:when test="[name='notequ']">
      <admst:return name="bname" value="!="/>
    </admst:when>
    <admst:when test="[name='lt']">
      <admst:return name="bname" value="&lt;"/>
    </admst:when>
    <admst:when test="[name='lt_equ']">
      <admst:return name="bname" value="&lt;="/>
    </admst:when>
    <admst:when test="[name='gt']">
      <admst:return name="bname" value="&gt;"/>
    </admst:when>
    <admst:when test="[name='gt_equ']">
      <admst:return name="bname" value="&gt;="/>
    </admst:when>
    <admst:when test="[name='shiftr']">
      <admst:return name="bname" value="&gt;&gt;"/>
    </admst:when>
    <admst:when test="[name='shiftl']">
      <admst:return name="bname" value="&lt;&lt;"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="variable type unknown\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- expression//function: get function name -->
<admst:template match="funcname">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)funcname*\n"/>
  <admst:choose>
    <admst:when test="[name='abs']"><admst:return name="fname" value="fabs"/></admst:when>
    <admst:when test="[name='\$shrinkl']"><admst:return name="fname" value="shrinkl"/></admst:when>
    <admst:when test="[name='\$shrinka']"><admst:return name="fname" value="shrinka"/></admst:when>
    <admst:when test="[name='log']"><admst:return name="fname" value="log10"/></admst:when>
    <admst:when test="[name='ln']"><admst:return name="fname" value="logE"/></admst:when>
    <admst:when test="[name='limexp']"><admst:return name="fname" value="limexp"/></admst:when>
    <admst:when test="[name='\$limexp']"><admst:return name="fname" value="limexp"/></admst:when>
    <admst:otherwise><admst:return name="fname" value="%(name)"/></admst:otherwise>
  </admst:choose>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="fname">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)fname*\n"/>
  <admst:choose>
    <admst:when test="[name='abs']">fabs</admst:when>
    <admst:when test="[name='\$shrinkl']">shrinkl</admst:when>
    <admst:when test="[name='\$shrinka']">shrinka</admst:when>
    <admst:when test="[name='log']">log10</admst:when>
    <admst:when test="[name='ln']">logE</admst:when>
    <admst:when test="[name='limexp']">limexp</admst:when>
    <admst:when test="[name='\$limexp']">limexp</admst:when>
    <admst:otherwise>%(name)</admst:otherwise>
  </admst:choose>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="e">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)e %(adms/datatypename)*\n"/>
  <admst:apply-templates select="." match="%(adms/datatypename)">$e
  </admst:apply-templates>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="ddx">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)ddx*\n"/>
  <admst:apply-templates select="." match="%(adms/datatypename)"/>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="ddxname">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)ddxname*\n"/>
  <admst:return name="ddxname"
    value="%(name)_%($pprobe/nature/access)%($pprobe/branch/pnode/name)_%($pprobe/branch/nnode/name)_%($qprobe/nature/access)%($qprobe/branch/pnode/name)_%($qprobe/branch/nnode/name)"/>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="dxname">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)dxname*\n"/>
  <admst:return name="dxname"
    value="%(name)_%($pprobe/nature/access)%($pprobe/branch/pnode/name)_%($pprobe/branch/nnode/name)"/>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="variable">
  <admst:choose>
    <admst:when test="[input='yes' and parametertype='model']">
      <admst:variable name="e" select="m-&gt;$(var_prefix)%(name)"/>
    </admst:when>
    <admst:when test="[input='yes' and parametertype='instance']">
      <admst:variable name="e" select="$(instancevar)$(var_prefix)%(name)"/>
    </admst:when>
    <admst:when test="[input='no' and scope='global_model']">
      <admst:variable name="e" select="m-&gt;$(var_prefix)%(name)"/>
    </admst:when>
    <admst:when test="[input='no' and scope='global_instance']">
      <admst:variable name="e" select="$(var_prefix)%(name)"/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="e" select="%(name)"/>
    </admst:otherwise>
  </admst:choose>

  <admst:variable name="ep" select="0.0"/>
  <admst:if test="[insource='yes']">
    <admst:if-inside select="$pprobe" list="%(probe)">
      <admst:variable name="ep" select="%(name)_%($pprobe/nature/access)%($pprobe/branch/pnode/name)_%($pprobe/branch/nnode/name)"/>
    </admst:if-inside>
  </admst:if>
  <admst:if test="$qprobe">

  <admst:variable name="eq" select="0.0"/>
  <admst:if test="[insource='yes']">
    <admst:if-inside select="$qprobe" list="%(probe)">
      <admst:variable name="eq" select="%(name)_%($qprobe/nature/access)%($qprobe/branch/pnode/name)_%($qprobe/branch/nnode/name)"/>
    </admst:if-inside>
  </admst:if>

  <admst:variable name="epq" select="0.0"/>
  <admst:if test="[insource='yes']">
    <admst:if-inside select="$pprobe" list="%(probe)">
      <admst:if-inside select="$qprobe" list="%(probe)">
        <admst:variable name="epq" select="%(ddxname(.)/[name='ddxname']/value)"/>
      </admst:if-inside>
    </admst:if-inside>
  </admst:if>
  </admst:if>
</admst:template>

<!-- -------------------------------------------------- -->
<admst:template match="probe">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)probe*\n"/>
  <admst:choose>
    <admst:when test="branch/nnode[grounded='no']">
      <admst:variable name="e" select="(BP(%(branch/pnode/name), %(branch/nnode/name)))"/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="e" select="NP(%(branch/pnode/name)) /*%(nature)*/"/>
    </admst:otherwise>
  </admst:choose>

  <admst:choose>
    <admst:when test="[.=$pprobe]">
      <admst:variable name="ep" select=" /*pprobe*/ $(ddtp) "/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="ep" select="0. /*pprobe*/"/>
    </admst:otherwise>
  </admst:choose>

  <admst:if test="$qprobe">
  <admst:choose>
    <admst:when test="[.=$qprobe]">
      <admst:variable name="eq" select=" /*qprobe*/$(ddtp) "/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="eq" select="0. /*qprobe*/ "/>
    </admst:otherwise>
  </admst:choose>

  <admst:variable name="epq" select="0.0"/>
  </admst:if>

</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="node">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)node*\n"/>
  <admst:fatal format="module node not expected here ... %(name)\n"/>

  <admst:fatal format="module node not expected here ... %(name)\n"/>
  <admst:if test="$qprobe">

  <admst:fatal format="module node not expected here ... %(name)\n"/>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="string">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)string*\n"/>
  <admst:variable name="e" select="&quot;%(value)&quot;"/>

  <admst:variable name="ep" select="0.0"/>
  <admst:if test="$qprobe">

  <admst:variable name="eq" select="0.0"/>

  <admst:variable name="epq" select="0.0"/>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="number">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)number*\n"/>
  <admst:choose>
    <admst:when test="[scalingunit='1']">
      <admst:variable name="e" select="%(value)"/>
    </admst:when>
    <admst:when test="[scalingunit='E']">
      <admst:variable name="e" select="(%(value)*1.0e+18)"/>
    </admst:when>
    <admst:when test="[scalingunit='P']">
      <admst:variable name="e" select="(%(value)*1.0e+15)"/>
    </admst:when>
    <admst:when test="[scalingunit='T']">
      <admst:variable name="e" select="(%(value)*1.0e+12)"/>
    </admst:when>
    <admst:when test="[scalingunit='G']">
      <admst:variable name="e" select="(%(value)*1.0e+9)"/>
    </admst:when>
    <admst:when test="[scalingunit='M']">
      <admst:variable name="e" select="(%(value)*1.0e+6)"/>
    </admst:when>
    <admst:when test="[scalingunit='k']">
      <admst:variable name="e" select="(%(value)*1.0e+3)"/>
    </admst:when>
    <admst:when test="[scalingunit='h']">
      <admst:variable name="e" select="(%(value)*1.0e+2)"/>
    </admst:when>
    <admst:when test="[scalingunit='D']">
      <admst:variable name="e" select="(%(value)*1.0e+1)"/>
    </admst:when>
    <admst:when test="[scalingunit='d']">
      <admst:variable name="e" select="(%(value)*1.0e-1)"/>
    </admst:when>
    <admst:when test="[scalingunit='c']">
      <admst:variable name="e" select="(%(value)*1.0e-2)"/>
    </admst:when>
    <admst:when test="[scalingunit='m']">
      <admst:variable name="e" select="(%(value)*1.0e-3)"/>
    </admst:when>
    <admst:when test="[scalingunit='u']">
      <admst:variable name="e" select="(%(value)*1.0e-6)"/>
    </admst:when>
    <admst:when test="[scalingunit='n']">
      <admst:variable name="e" select="(%(value)*1.0e-9)"/>
    </admst:when>
    <admst:when test="[scalingunit='A']">
      <admst:variable name="e" select="(%(value)*1.0e-10)"/>
    </admst:when>
    <admst:when test="[scalingunit='p']">
      <admst:variable name="e" select="(%(value)*1.0e-12)"/>
    </admst:when>
    <admst:when test="[scalingunit='f']">
      <admst:variable name="e" select="(%(value)*1.0e-15)"/>
    </admst:when>
    <admst:when test="[scalingunit='a']">
      <admst:variable name="e" select="(%(value)*1.0e-18)"/>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="scaling unit not supported: %(scalingunit)\n"/>
    </admst:otherwise>
  </admst:choose>

  <admst:variable name="ep" select="0.0"/>
  <admst:if test="$qprobe">

  <admst:variable name="eq" select="0.0"/>

  <admst:variable name="epq" select="0.0"/>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="mapply_unary">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)mapply_unary*\n"/>
  <admst:choose>
    <admst:when test="[name='plus']">
      <admst:choose>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="ddx"/>
          <admst:variable name="e" select="(+$e)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='minus']">
      <admst:choose>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="ddx"/>
          <admst:variable name="e" select="(-$e)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='not']">
      <admst:choose>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="ddx"/>
          <admst:variable name="e" select="(!$e)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='bw_not']">
      <admst:choose>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:apply-templates select="arg1" match="ddx"/>
          <admst:variable name="e" select="(~$e)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:fatal format="%(name): function not handled\n"/>
    </admst:otherwise>
  </admst:choose>

  <admst:choose>
    <admst:when test="[$e='0.0']">
      <admst:variable name="ep" select="0.0"/>
    </admst:when>
    <admst:when test="[$ep='0.0']">
      <admst:variable name="ep" select="0.0"/>
    </admst:when>
    <admst:otherwise>
      <admst:choose>
        <admst:when test="[name='plus']">
          <admst:variable name="ep" select="(+$ep)"/>
        </admst:when>
        <admst:when test="[name='minus']">
          <admst:variable name="ep" select="(-$ep)"/>
        </admst:when>
        <admst:when test="[name='not']">
          <admst:variable name="ep" select="(!$ep)"/>
        </admst:when>
        <admst:when test="[name='bw_not']">
          <admst:variable name="ep" select="(~$ep)"/>
        </admst:when>
      </admst:choose>
    </admst:otherwise>
  </admst:choose>
  <admst:if test="$qprobe">

  <admst:choose>
    <admst:when test="[$e='0.0']">
      <admst:variable name="eq" select="0.0"/>
    </admst:when>
    <admst:when test="[$eq='0.0']">
      <admst:variable name="eq" select="0.0"/>
    </admst:when>
    <admst:otherwise>
      <admst:choose>
        <admst:when test="[name='plus']">
          <admst:variable name="eq" select="(+$eq)"/>
        </admst:when>
        <admst:when test="[name='minus']">
          <admst:variable name="eq" select="(-$eq)"/>
        </admst:when>
        <admst:when test="[name='not']">
          <admst:variable name="eq" select="(!$eq)"/>
        </admst:when>
        <admst:when test="[name='bw_not']">
          <admst:variable name="eq" select="(~$eq)"/>
        </admst:when>
      </admst:choose>
    </admst:otherwise>
  </admst:choose>

  <admst:variable name="epq" select="0.0"/>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="mapply_binary">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)mapply_binary*\n"/>
  <admst:apply-templates select="arg1" match="ddx"/>
  <admst:variable name="x" select="$e"/>
  <admst:variable name="xp" select="$ep"/>
  <admst:variable name="xq" select="$eq"/>
  <admst:variable name="xpq" select="$epq"/>
  <admst:apply-templates select="arg2" match="ddx"/>
  <admst:variable name="y" select="$e"/>
  <admst:variable name="yp" select="$ep"/>
  <admst:variable name="yq" select="$eq"/>
  <admst:variable name="ypq" select="$epq"/>
  <admst:choose>
    <admst:when test="[name='addp']">
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)and(arg2/math/value=0.0)]">
          <admst:variable name="e" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="(+$y)"/>
          <admst:variable name="xp" select="0.0"/>
        </admst:when>
        <admst:when test="[arg2/math/value=0.0]">
          <admst:variable name="e" select="%($x)"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:otherwise>
/* ... */
          <admst:variable name="e" select="($x+$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='addm']">
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)and(arg2/math/value=0.0)]">
          <admst:variable name="e" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="(-$y)"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="$ep"/>
        </admst:when>
        <admst:when test="arg2/math[value=0.0]">
          <admst:variable name="e" select="$x"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="e" select="($x-$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='multtime']">
      <admst:choose>
        <admst:when test="[(arg1/math/value=0.0)or(arg2/math/value=0.0)]">
          <admst:variable name="e" select="0.0"/>
          <admst:variable name="x" select="0.0"/>
          <admst:variable name="y" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:when test="[(arg1/math/value=1.0)and(arg2/math/value=1.0)]">
          <admst:variable name="e" select="1.0"/>
          <admst:variable name="x" select="0.0"/>
          <admst:variable name="y" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="e" select="($x*$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='multdiv']">
      <admst:choose>
        <admst:when test="[arg1/math/value=0.0]">
          <admst:variable name="e" select="0.0"/>
          <admst:variable name="x" select="0.0"/>
          <admst:variable name="y" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:when test="[(arg1/math/value=1.0)and(arg2/math/value=1.0)]">
          <admst:variable name="e" select="1.0"/>
          <admst:variable name="x" select="0.0"/>
          <admst:variable name="y" select="0.0"/>
          <admst:variable name="xp" select="0.0"/>
          <admst:variable name="yp" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="e" select="($x/$y)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="e" select="($x%(bname(.)/[name='bname']/value)$y)"/>
    </admst:otherwise>
  </admst:choose>

    <admst:choose>
      <admst:when test="[name='addp']">
        <admst:choose>
          <admst:when test="[$xp='0.0' and $yp='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$xp='0.0']">
            <admst:variable name="ep" select="$yp"/>
          </admst:when>
          <admst:when test="[$yp='0.0']">
            <admst:variable name="ep" select="$xp"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="ep" select="($xp+$yp)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='addm']">
        <admst:choose>
          <admst:when test="[$xp='0.0' and $yp='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$xp='0.0']">
            <admst:variable name="ep" select="(-$yp)"/>
          </admst:when>
          <admst:when test="[$yp='0.0']">
            <admst:variable name="ep" select="$xp"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="ep" select="($xp-$yp)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multtime']">
        <admst:choose>
          <admst:when test="[$x='0.0' and $y='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$xp='0.0' and $yp='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$xp='0.0' and $yp='1.0']">
            <admst:variable name="ep" select="($x)"/>
          </admst:when>
          <admst:when test="[$xp='1.0' and $yp='0.0']">
            <admst:variable name="ep" select="($y)"/>
          </admst:when>
          <admst:when test="[$xp='0.0']">
            <admst:variable name="ep" select="($x*$yp)"/>
          </admst:when>
          <admst:when test="[$yp='0.0']">
            <admst:variable name="ep" select="$xp*$y"/>
          </admst:when>
          <admst:when test="[$xp='1.0' and $yp='1.0']">
            <admst:variable name="ep" select="($x+$y)"/>
          </admst:when>
          <admst:when test="[$xp='1.0']">
            <admst:variable name="ep" select="($y+($yp*$x))"/>
          </admst:when>
          <admst:when test="[$yp='1.0']">
            <admst:variable name="ep" select="(($xp*$y)+$x)"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:variable name="ep" select="($yp)"/>
          </admst:when>
          <admst:when test="[$y='1.0']">
            <admst:variable name="ep" select="( $xp)"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="ep" select="(($xp*$y)+($x*$yp))"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multdiv']">
        <admst:choose>
          <admst:when test="[$x='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$xp='0.0' and $yp='0.0']">
            <admst:variable name="ep" select="0.0"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:choose>
              <admst:when test="[$yp='1.0']">
                <admst:variable name="ep" select="(-1/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="ep" select="(-$yp/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$xp='0.0']">
            <admst:choose>
              <admst:when test="[$yp='1.0']">
                <admst:variable name="ep" select="(-$x/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="ep" select="(-$x*$yp/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$xp='1.0']">
            <admst:choose>
              <admst:when test="[$yp='0.0']">
                <admst:variable name="ep" select="(1/$y)"/>
              </admst:when>
              <admst:when test="[$yp='1.0']">
                <admst:variable name="ep" select="(($y-$x)/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="ep" select="(($y-($x*$yp))/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:otherwise>
            <admst:choose>
              <admst:when test="[$y='1.0']">
                <admst:variable name="ep" select="( $xp )"/>
              </admst:when>
              <admst:when test="[$yp='0.0']">
                <admst:variable name="ep" select="($xp/$y)"/>
              </admst:when>
              <admst:when test="[$yp='1.0']">
                <admst:variable name="ep" select="(($xp*$y-$x)/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="ep" select="(($xp*$y-$x*$yp)/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:otherwise>
        <admst:variable name="ep" select="0.0"/>
      </admst:otherwise>
    </admst:choose>
  <admst:if test="$qprobe">

    <admst:choose>
      <admst:when test="[name='addp']">
        <admst:choose>
          <admst:when test="[$xq='0.0' and $yq='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$xq='0.0']">
            <admst:variable name="eq" select="$yq"/>
          </admst:when>
          <admst:when test="[$yq='0.0']">
            <admst:variable name="eq" select="$xq"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="eq" select="($xq+$yq)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='addm']">
        <admst:choose>
          <admst:when test="[$xq='0.0' and $yq='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$xq='0.0']">
            <admst:variable name="eq" select="(-$yq)"/>
          </admst:when>
          <admst:when test="[$yq='0.0']">
            <admst:variable name="eq" select="$xq"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="eq" select="($xq-$yq)"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multtime']">
        <admst:choose>
          <admst:when test="[$x='0.0' and $y='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$xq='0.0' and $yq='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$xq='0.0' and $yq='1.0']">
            <admst:variable name="eq" select="($x)"/>
          </admst:when>
          <admst:when test="[$xq='1.0' and $yq='0.0']">
            <admst:variable name="eq" select="($y)"/>
          </admst:when>
          <admst:when test="[$xq='0.0']">
            <admst:variable name="eq" select="($x*$yq)"/>
          </admst:when>
          <admst:when test="[$yq='0.0']">
            <admst:variable name="eq" select="$xq*$y"/>
          </admst:when>
          <admst:when test="[$xq='1.0' and $yq='1.0']">
            <admst:variable name="eq" select="($x+$y)"/>
          </admst:when>
          <admst:when test="[$xq='1.0']">
            <admst:variable name="eq" select="($y+($yq*$x))"/>
          </admst:when>
          <admst:when test="[$yq='1.0']">
            <admst:variable name="eq" select="(($xq*$y)+$x)"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:variable name="eq" select="$yq"/>
          </admst:when>
          <admst:when test="[$y='1.0']">
            <admst:variable name="eq" select="$xq"/>
          </admst:when>
          <admst:otherwise>
            <admst:variable name="eq" select="(($xq*$y)+($x*$yq))"/>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:when test="[name='multdiv']">
        <admst:choose>
          <admst:when test="[$x='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$xq='0.0' and $yq='0.0']">
            <admst:variable name="eq" select="0.0"/>
          </admst:when>
          <admst:when test="[$x='1.0']">
            <admst:choose>
              <admst:when test="[$yq='1.0']">
                <admst:variable name="eq" select="(-1/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="eq" select="(-$yq/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$xq='0.0']">
            <admst:choose>
              <admst:when test="[$yq='1.0']">
                <admst:variable name="eq" select="(-$x/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="eq" select="(-$x*$yq/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:when test="[$xq='1.0']">
            <admst:choose>
              <admst:when test="[$yq='0.0']">
                <admst:variable name="eq" select="(1/$y)"/>
              </admst:when>
              <admst:when test="[$yq='1.0']">
                <admst:variable name="eq" select="(($y-$x)/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="eq" select="(($y-($x*$yq))/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:when>
          <admst:otherwise>
            <admst:choose>
              <admst:when test="[$y='1.0']">
                <admst:variable name="eq" select="$xq"/>
              </admst:when>
              <admst:when test="[$yq='0.0']">
                <admst:variable name="eq" select="($xq/$y)"/>
              </admst:when>
              <admst:when test="[$yq='1.0']">
                <admst:variable name="eq" select="(($xq*$y-$x)/$y/$y)"/>
              </admst:when>
              <admst:otherwise>
                <admst:variable name="eq" select="(($xq*$y-$x*$yq)/$y/$y)"/>
              </admst:otherwise>
            </admst:choose>
          </admst:otherwise>
        </admst:choose>
      </admst:when>
      <admst:otherwise>
        <admst:variable name="eq" select="0.0"/>
      </admst:otherwise>
    </admst:choose>

    <admst:choose>
      <admst:when test="[name='addp']">
          <admst:variable name="t1" select="+$xpq"/>
          <admst:variable name="t2" select="+$ypq"/>
          <admst:variable name="epq" select="$t1$t2"/>
          <admst:choose>
            <admst:when test="[$epq='']">
              <admst:variable name="epq" select="0.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="epq" select="($epq)"/>
            </admst:otherwise>
          </admst:choose>
      </admst:when>
      <admst:when test="[name='addm']">
          <admst:variable name="t1" select="+$xpq"/>
          <admst:variable name="t2" select="-$ypq"/>
          <admst:variable name="epq" select="$t1$t2"/>
          <admst:choose>
            <admst:when test="[$epq='']">
              <admst:variable name="epq" select="0.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="epq" select="($epq)"/>
            </admst:otherwise>
          </admst:choose>
      </admst:when>
      <admst:when test="[name='multtime']">
          <admst:variable name="t1" select="+$xpq*$y"/>
          <admst:variable name="t2" select="+$xp*$yq"/>
          <admst:variable name="t3" select="+$xq*$yp"/>
          <admst:variable name="t4" select="+$x*$ypq"/>
          <admst:variable name="epq" select="$t1$t2$t3$t4"/>
          <admst:choose>
            <admst:when test="[$eq='']">
              <admst:variable name="eq" select="0.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="eq" select="($eq)"/>
            </admst:otherwise>
          </admst:choose>
          <admst:choose>
            <admst:when test="[$epq='']">
              <admst:variable name="epq" select="0.0"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="epq" select="($epq)"/>
            </admst:otherwise>
          </admst:choose>
      </admst:when>
      <admst:when test="[name='multdiv']">
        <admst:variable name="epq" select="($xpq/$y-($xp*$yq+$xq*$yp+$x*$ypq)/$y/$y+2*$x*$yp*$yq/$y/$y/$y)"/>
      </admst:when>
    </admst:choose>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="mapply_ternary">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)mapply_ternary*\n"/>
  <admst:apply-templates select="arg1" match="ddx"/>
  <admst:variable name="x" select="$e"/>
  <admst:apply-templates select="arg2" match="ddx"/>
  <admst:variable name="y" select="$e"/>
  <admst:variable name="yp" select="$ep"/>
  <admst:variable name="yq" select="$eq"/>
  <admst:apply-templates select="arg3" match="ddx"/>
  <admst:variable name="z" select="$e"/>
  <admst:variable name="zp" select="$ep"/>
  <admst:variable name="zq" select="$eq"/>
  <admst:variable name="e" select="($x?$y:$z)"/>

  <admst:variable name="ep" select="($x?$yp:$zp)"/>
  <admst:if test="$qprobe">

  <admst:variable name="ep" select="($x?$yp:$zp)"/>

  <admst:variable name="epq" select="fixme"/>
  </admst:if>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="DDT">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)DDT*\n"/>
  <admst:variable name="ddt_cur" select="$(ddt_cur)_"/>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="function">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)function*\n"/>
  <admst:choose>
    <admst:when test="[name='absdelay']">
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="e" select="$e"/>
    </admst:when>
    <admst:when test="[name='ddt']">


      <!-- --------------i
      <admst:text  format="//before  DDT: ep $ep\n"/>
      <admst:variable name="xp"/>
      <admst:apply-templates select="." match="DDT"/> ------ -->

			/* function ddt %(.)*/

      <admst:variable name="ddtp" select="keep_%(adms/id(.))"/>
      <admst:variable name="pddt" select="%(adms/id(.))"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="ddtp" select="1."/>
      <admst:variable name="x" select="$e"/>
      <admst:variable name="xp" select="$ep"/>
      <admst:variable name="xq" select="$eq"/>
      <admst:variable name="xpq" select="$epq"/>
#if 0
			* x $x
			* xp $xp 
			* xq $xq
			* xpq $xpq
#endif

      <admst:variable name="e" select="_DDT($x, _ddt_%(adms/id(.)) ) "/>

      <!-- <admst:apply-templates select="." match="DDT"/>  -->

      <admst:variable name="ddt_cur" select="$(ddt_cur)_"/>


      <!-- ------------------------ -->

    </admst:when>
    <admst:when test="[name='\$given']">
      <admst:variable name="arg1" select="%(arguments[1])"/>
      <admst:assert test="$arg1/adms[datatypename='variable']" format="\$given: argument is not a variable\n"/>
      <admst:assert test="$arg1/[input='yes']" format="\$given(%(name)): argument is not a parameter\n"/>
      <admst:choose>
        <admst:when test="$arg1/[parametertype='model']">
          <admst:variable name="e" select="_mpg(%(name),%(name)Given)"/>
        </admst:when>
        <admst:when test="$arg1/[parametertype='instance']">
          <admst:variable name="e" select="_ipg(%(name),%(name)Given)"/>
        </admst:when>
        <admst:otherwise>
          <admst:fatal format="$given(): should not be reached\n"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='\$temperature']">
      <admst:assert test="[nilled(arguments)]" format="%(name): should not have arguments\n"/>
      <admst:variable name="e" select="(_ambient_temp)"/>
    </admst:when>
    <admst:when test="[name='\$mfactor']">
      <admst:assert test="[nilled(arguments)]" format="%(name): should not have arguments\n"/>
      <admst:variable name="e" select="MFACTOR"/>
    </admst:when>
    <admst:when test="[name='\$vt']">
      <admst:choose>
        <admst:when test="[nilled(arguments)]">
          <admst:variable name="e" select="(BOLTZMANN*(_circuit_tnom)/ELECTRON_CHARGE)"/>
        </admst:when>
        <admst:when test="arguments[count(.)=1]">
          <admst:apply-templates select="arguments[1]" match="ddx"/>
          <admst:variable name="x" select="$e"/>
          <admst:variable name="xp" select="$ep"/>
          <admst:variable name="xq" select="$eq"/>
          <admst:variable name="xpq" select="$epq"/>
          <admst:variable name="e" select="(BOLTZMANN*$x/ELECTRON_CHARGE)"/>
        </admst:when>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='\$scale']">
      <admst:assert test="[nilled(arguments)]" format="%(name): should not have arguments\n"/>
      <admst:variable name="e" select="_scale"/>
    </admst:when>
    <admst:when test="[name='\$abstime']">
      <admst:assert test="[nilled(arguments)]" format="%(name): should not have arguments\n"/>
      <admst:variable name="e" select="_abstime"/>
    </admst:when>
    <admst:when test="[name='ddx']">
      <admst:assert test="arguments[count(.)=2]" format="%(name): should have two arguments exactly\n"/>
      <admst:assert test="arguments[2]/adms[datatypename='probe']" format="%(name): second argument is not a probe\n"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="e" select="0.0 /*ddx should be top node of expression!*/"/>
    </admst:when>
    <admst:when test="[name='floor']">
      <admst:assert test="arguments[count(.)=1]" format="%(name): should have one argument exactly\n"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="e" select="floor($e)"/>
    </admst:when>
    <admst:when test="[name='ceil']">
      <admst:assert test="arguments[count(.)=1]" format="%(name): should have one argument exactly\n"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="e" select="ceil($e)"/>
    </admst:when>
    <admst:when test="[$SkipFVariable='y']">
      <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
      <admst:variable name="args" select=""/>
      <admst:for-each select="arguments">
        <admst:if test="[$args!='']">
          <admst:variable name="args" select="$args,"/>
        </admst:if>
        <admst:apply-templates select="." match="ddx"/>
        <admst:variable name="args" select="$args$e"/>
      </admst:for-each>
      <admst:variable name="e" select="$(fname)($args)"/>
    </admst:when>
    <admst:when test="[name='abs' or name='acos' or name='asin' or name='atan'
                    or name='cos' or name='cosh' or name='exp' or name='hypot'
                    or name='limexp' or name='ln' or name='log' or name='sin'
                    or name='sinh' or name='sqrt' or name='tan' or name='tanh']">
      <admst:assert test="arguments[count(.)=1]" format="%(name): should have one argument exactly\n"/>
      <admst:variable name="index" select="%(index(subexpression/expression/function,.))"/>
      <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="x" select="$e"/>
      <admst:variable name="xp" select="$ep"/>
      <admst:variable name="xq" select="$eq"/>
      <admst:variable name="xpq" select="$epq"/>
      <admst:variable name="e" select="d00_$(fname)$index"/>
    </admst:when>
    <admst:when test="[name='div' or name='pow' or name='hypot' or name='min' or name='max']">
      <admst:assert test="arguments[count(.)=2]" format="%(name): should have two argument exactly\n"/>
      <admst:variable name="index" select="%(index(./subexpression/expression/function,.))"/>
      <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
      <admst:apply-templates select="arguments[1]" match="ddx"/>
      <admst:variable name="x" select="$e"/>
      <admst:variable name="xp" select="$ep"/>
      <admst:variable name="xq" select="$eq"/>
      <admst:variable name="xpq" select="$epq"/>
      <admst:apply-templates select="arguments[2]" match="ddx"/>
      <admst:variable name="y" select="$e"/>
      <admst:variable name="yp" select="$ep"/>
      <admst:variable name="yq" select="$eq"/>
      <admst:variable name="ypq" select="$epq"/>
      <admst:variable name="e" select="d00_$(fname)$index"/>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
      <admst:variable name="args" select=""/>
      <admst:for-each select="arguments">
        <admst:if test="[$args!='']">
          <admst:variable name="args" select="$args,"/>
        </admst:if>
        <admst:apply-templates select="." match="ddx"/>
        <admst:variable name="args" select="$args$e"/>
      </admst:for-each>
      <admst:variable name="dargs" select=""/>
      <admst:for-each select="arguments">
        <admst:if test="[$dargs!='']">
          <admst:variable name="dargs" select="$dargs,"/>
        </admst:if>
        <admst:apply-templates select="." match="ddx"/>
        <admst:variable name="dargs" select="$dargs$ep"/>
      </admst:for-each>
      <admst:variable name="e" select="$(fname)($args)"/>
    </admst:otherwise>
  </admst:choose>

  <!--  ???? ----- -->
  <admst:variable name="ep" select="0.0"/>
  <admst:choose>
    <admst:when test="[name='absdelay']">
    </admst:when>
    <admst:when test="[name='\$given']">
    </admst:when>
    <admst:when test="[name='\$temperature']">
    </admst:when>
    <admst:when test="[name='\$mfactor']">
    </admst:when>
    <admst:when test="[name='\$vt']">
      <admst:choose>
        <admst:when test="[nilled(arguments)]">
          <admst:variable name="ep" select="0.0"/>
        </admst:when>
        <admst:when test="arguments[count(.)=1]">
          <admst:variable name="ep" select="(BOLTZMANN*$xp/ELECTRON_CHARGE)"/>
        </admst:when>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='\$scale']">
    </admst:when>
    <admst:when test="[name='\$abstime']">
    </admst:when>
    <admst:when test="[name='ddx']">
    </admst:when>
    <admst:when test="[name='floor']">
    </admst:when>
    <admst:when test="[name='ceil']">
    </admst:when>
    <admst:when test="[$SkipFVariable='y']">
    </admst:when>
    <admst:when test="[name='ddt']">
      <admst:variable name="ep" select="$xp"/> <!-- ??? -->
    </admst:when>
    <admst:when test="[name='abs' or name='acos' or name='asin' or name='atan'
                    or name='cos' or name='cosh' or name='exp' or name='hypot'
                    or name='limexp' or name='ln' or name='log' or name='sin'
                    or name='sinh' or name='sqrt' or name='tan' or name='tanh']">
      <admst:variable name="index" select="%(index(subexpression/expression/function,.))"/>
      <admst:choose>
        <admst:when test="[$xp='0.0']">
          <admst:variable name="ep" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="ep" select="$xp*d10_$(fname)$index"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='div' or name='pow' or name='hypot' or name='min' or name='max']">
      <admst:variable name="index" select="%(index(./subexpression/expression/function,.))"/>
      <admst:choose>
        <admst:when test="[$xp='0.0' and $yp='0.0']">
          <admst:variable name="ep" select="0.0"/>
        </admst:when>
        <admst:when test="[$xp='0.0']">
          <admst:variable name="ep" select="(d11_$(fname)$index*$yp)"/>
        </admst:when>
        <admst:when test="[$yp='0.0']">
          <admst:variable name="ep" select="(d10_$(fname)$index*$xp)"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="ep" select="(d10_$(fname)$index*$xp+d11_$(fname)$index*$yp)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="ep" select="d_$(fname)($args,$dargs)"/>
    </admst:otherwise>
  </admst:choose>
  <admst:if test="$qprobe">

  <admst:variable name="eq" select="0.0"/>
  <admst:choose>
    <admst:when test="[name='absdelay']">
    </admst:when>
    <admst:when test="[name='\$given']">
    </admst:when>
    <admst:when test="[name='\$temperature']">
    </admst:when>
    <admst:when test="[name='\$mfactor']">
    </admst:when>
    <admst:when test="[name='\$vt']">
      <admst:choose>
        <admst:when test="[nilled(arguments)]">
          <admst:variable name="eq" select="0.0"/>
        </admst:when>
        <admst:when test="arguments[count(.)=1]">
          <admst:variable name="eq" select="(BOLTZMANN*$xq/ELECTRON_CHARGE)"/>
        </admst:when>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='\$scale']">
    </admst:when>
    <admst:when test="[name='\$abstime']">
    </admst:when>
    <admst:when test="[name='ddx']">
    </admst:when>
    <admst:when test="[name='floor']">
    </admst:when>
    <admst:when test="[name='ceil']">
    </admst:when>
    <admst:when test="[$SkipFVariable='y']">
    </admst:when>
    <admst:when test="[name='ddt']">
      <admst:variable name="eq" select="$xq"/>
    </admst:when>
    <admst:when test="[name='abs' or name='acos' or name='asin' or name='atan'
                       or name='cos' or name='cosh' or name='exp' or name='hypot'
                       or name='limexp' or name='ln' or name='log' or name='sin' or name='sinh'
                       or name='sqrt' or name='tan' or name='tanh']">
      <admst:variable name="index" select="%(index(subexpression/expression/function,.))"/>
      <admst:choose>
        <admst:when test="[$xq='0.0']">
          <admst:variable name="eq" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="eq" select="$xq*d10_$(fname)$index"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='div' or name='pow' or name='hypot' or name='min' or name='max']">
      <admst:variable name="index" select="%(index(./subexpression/expression/function,.))"/>
      <admst:choose>
        <admst:when test="[$xq='0.0' and $yq='0.0']">
          <admst:variable name="eq" select="0.0"/>
        </admst:when>
        <admst:when test="[$xq='0.0']">
          <admst:variable name="eq" select="(d11_$(fname)$index*$yq)"/>
        </admst:when>
        <admst:when test="[$yq='0.0']">
          <admst:variable name="eq" select="(d10_$(fname)$index*$xq)"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="eq" select="(d10_$(fname)$index*$xq+d11_$(fname)$index*$yq)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="eq" select="d_$(fname)($args,$dargs)"/>
    </admst:otherwise>
  </admst:choose>

  <admst:choose>
    <admst:when test="[name='abs' or name='acos' or name='asin' or name='atan'
                       or name='cos' or name='cosh' or name='exp' or name='hypot'
                       or name='limexp' or name='ln' or name='log' or name='sin'
                       or name='sinh' or name='sqrt' or name='tan' or name='tanh']">
      <admst:variable name="index" select="%(index(./subexpression/expression/function,.))"/>
      <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
      <admst:choose>
        <admst:when test="[$x='0.0']">
          <admst:variable name="epq" select="0.0"/>
        </admst:when>
        <admst:when test="[$xp='0.0']">
          <admst:variable name="epq" select="0.0"/>
        </admst:when>
        <admst:when test="[$xq='0.0']">
          <admst:variable name="epq" select="0.0"/>
        </admst:when>
        <admst:otherwise>
          <admst:variable name="epq" select="(m20_$(fname)($x)*$xq*$xp+d10_$(fname)$index*$xpq)"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:when test="[name='div' or name='pow' or name='hypot' or name='min' or name='max']">
      <admst:variable name="epq" select="fixme"/>
      <admst:if test="[$requiredderivateforddx='yes']">
        <admst:warning format="%(name): ddx dependency not implemented\n"/>
      </admst:if>
    </admst:when>
  </admst:choose>
  </admst:if>
</admst:template>

<!-- ------------------ analog//block -------------------- -->
<admst:template match="block">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)block*\n"/>
	{ //block %(datatypename) %(name)

  <admst:apply-templates select="item" match="%(adms/datatypename)"/>
	} //block

</admst:template>
<!-- --------------------- analog//blockvariable ------------ -->
<admst:template match="blockvariable">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)blockvariable*\n"/>
  <admst:text format="// { //blockvariable \n"/>
  <admst:text select="variable" format="%(vtype(.)) %(name);\n"/>
  <admst:if test="variable[insource='yes']/probe">

    <admst:for-each select="variable">
      <admst:variable name="myvariable" select="%(.)"/>
      <admst:for-each select="probe">
          <admst:variable name="pprobe" select="%(.)"/>
          <admst:variable name="ddxinsidethisprobe" select="no"/>
          <admst:if test="$myvariable/ddxprobe/branch/pnode[.=$pprobe/branch/pnode or .=$pprobe/branch/nnode]">
            <admst:variable name="ddxinsidethisprobe" select="yes"/>
          </admst:if>

          <admst:if test="[($ddxinsidederivate='yes' and $_DERIVATEFORDDX='yes'  )
             or ($ddxinsidederivate='no'  and $_DERIVATE='yes' )]">
            <admst:text format="double %($myvariable/name)_%(nature/access)%(branch/pnode/name)_%(branch/nnode/name);\n"/>
          </admst:if>

      </admst:for-each>
    </admst:for-each>

    <admst:for-each select="variable">
      <admst:variable name="myvariable" select="%(.)"/>
      <admst:new datatype="list" arguments="list of ddx probes">
        <admst:variable name="ddxprobes" select="%(.)"/>
        <admst:for-each select="$myvariable/probe">
          <admst:variable name="pprobe" select="%(.)"/>
          <admst:push into="$ddxprobes/item"
            select="$myvariable/ddxprobe/branch/pnode[.=$pprobe/branch/pnode
                                     or .=$pprobe/branch/nnode]/$pprobe" onduplicate="ignore"/>
        </admst:for-each>
      </admst:new>
      <admst:text test="$ddxprobes/item" format="#if defined(_DERIVATE) // FIXME C\n"/>
      <admst:for-each select="$ddxprobes/item">
        <admst:variable name="pprobe" select="%(.)"/>
        <admst:for-each select="$myvariable/probe">
          <admst:variable name="qprobe" select="%(.)"/>
          <admst:text format="  double %(ddxname($myvariable)/[name='ddxname']/value);\n"/>
        </admst:for-each>
      </admst:for-each>
      <admst:text test="$ddxprobes/item" format="#endif\n"/>
    </admst:for-each>
 </admst:if>

 // } // lockvariable  FIXME \\n

</admst:template>

<!-- ------ analog//function: ddx handling ---- -->
<admst:template match="function:precomputationyes">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)function:precomputationyes*\n"/>
  <admst:text format="//precq\n"/>
  <admst:variable name="index" select="%(index(../function,.))"/>
  <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
  //precompyes
  %(e(.))
</admst:template>
<!-- ------ analog//function precomp ---- -->
<admst:template match="function:precomputation">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)function:precomputation*\n"/>
  <admst:text format="//prec \n"/>
  <admst:variable name="index" select="%(index(../function,.))"/>
  <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
  <admst:choose>
    <admst:when test="[name='add']">
      <admst:text format="m00_add(d00_add$index"/>
    </admst:when>
    <admst:when test="[name='div']">
      <admst:text format="m00_div(d00_div$index,d10_$(fname)$index"/>
    </admst:when>
    <admst:when test="[name='mult']">
      <admst:text format="m00_mult(d00_mult$index,d10_mult$index,d11_mult$index"/>
    </admst:when>
    <admst:otherwise>
      <admst:text format="double m00_$(fname)(d00_$(fname)$index"/>
    </admst:otherwise>
  </admst:choose>
  <admst:text select="arguments" format=",%(e(.))"/>
  <admst:text format=");\n"/>
</admst:template>
<!-- ------ analog//function: ddx handling ---- -->
<admst:template match="function:derivate:precomputationyes">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)function:derivate:precomputationyes*\n"/>
  <admst:if test="[hasVoltageDependentFunction='yes']">
    <admst:if test="[($ddxinsidederivate='yes' and $_DERIVATEFORDDX='yes'  )
                  or ($ddxinsidederivate='no'  and $_DERIVATE='yes' )]">
      <admst:for-each select="function">
        <admst:variable name="index" select="%(index(../function,.))"/>
        <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
        <admst:choose>
          <admst:when test="[name='exp']">
          </admst:when>
          <admst:when test="[name='add']">
          </admst:when>
          <admst:when test="[name='mult']"/>
          <admst:when test="[name='add']"/>
          <admst:when test="[name='div']">
            <admst:for-each select="arguments">
              <admst:variable name="position" select="%(position(.)-1)"/>
              <admst:if test="math[dependency!='constant']">
                <admst:text select="../arguments" format="%(e(.))"/>
              </admst:if>
            </admst:for-each>
          </admst:when>
          <admst:otherwise>
            <admst:for-each select="arguments">
              <admst:variable name="position" select="%(position(.)-1)"/>
              <admst:if test="math[dependency!='constant']">
                <admst:text select="../arguments" format="%(e(.))"/>
              </admst:if>
            </admst:for-each>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
    </admst:if>
  </admst:if>
</admst:template>
<!-- ------------ -->
<admst:template match="function:derivate:precomputation">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)function:derivate:precomputation*\n"/>
  <admst:if test="[hasVoltageDependentFunction='yes']">
    <!--  <admst:choose>
      <admst:when test="[$ddxinsidederivate='yes']">
        <admst:text format="#if defined(_DERIVATEFORDDX)\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:text format="#if defined(_DERIVATE) // BLA\n"/>
      </admst:otherwise>
    </admst:choose> FIXME -->
    <admst:if test="[($ddxinsidederivate='yes' and $_DERIVATEFORDDX='yes'  )
                  or ($ddxinsidederivate='no'  and $_DERIVATE='yes' )]">
      <admst:for-each select="function">
        <admst:variable name="index" select="%(index(../function,.))"/>
        <admst:variable name="fname" select="%(funcname(.)/[name='fname']/value)"/>
        <admst:choose>
          <admst:when test="[name='exp']">
            <admst:if test="arguments/math[dependency!='constant']">
              <admst:text format="#define d10_exp$index d00_exp$index\n"/>
            </admst:if>
          </admst:when>
          <admst:when test="[name='add']">
            <admst:if test="arguments/math[dependency!='constant']">
              <admst:text format="#define d10_add$index 1\n"/>
              <admst:text format="#define d11_add$index 1\n"/>
            </admst:if>
          </admst:when>
          <admst:when test="[name='mult']"/>
          <admst:when test="[name='add']"/>
          <admst:when test="[name='div']">
            <admst:for-each select="arguments">
              <admst:variable name="position" select="%(position(.)-1)"/>
              <admst:if test="math[dependency!='constant']">
                <admst:text
                  format="m1$(position)_$(fname)(d1$(position)_$(fname)$index,d00_$(fname)$index,d10_$(fname)$index"/>
                <admst:text select="../arguments" format=",%(e(.))"/>
                <admst:text format=")\n"/>
              </admst:if>
            </admst:for-each>
          </admst:when>
          <admst:otherwise>
            <admst:for-each select="arguments">
              <admst:variable name="position" select="%(position(.)-1)"/>
              <admst:if test="math[dependency!='constant']">
                <admst:text format="double m1%(position(.)-1)_$(fname)(d1%(position(.)-1)_$(fname)$index,d00_$(fname)$index"/>
                <admst:text select="../arguments" format=",%(e(.))"/>
                <admst:text format="); // 7539\n"/>
              </admst:if>
            </admst:for-each>
          </admst:otherwise>
        </admst:choose>
      </admst:for-each>
    </admst:if>
  </admst:if>
</admst:template>
<!-- ------------------ analog//assignment:ddt ------------------------ -->
<admst:template match="assignment:ddt">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)assignment:ddt*\n"/>
</admst:template>
<!-- ---------------------------------------------------------- -->
<admst:template match="assignment">
  <admst:message test="[/dbg_xml='yes']" format=" assignment from ngspiceMODULE.hxx.xml?"/>
  /* MODULE assignment */
  <admst:variable name="assignment" select="%(.)"/>
  <admst:variable name="rhs" select="%(rhs)"/>
  <admst:variable name="lhs" select="%(lhs)"/>
  <!--  <admst:text test="[dynamic='yes']" format="#if defined(_DYNAMIC)\n"/> -->
  <admst:if test="[dynamic='no' or $_DYNAMIC='yes']" >
    <admst:text test="rhs/function" format="{ // rhs function \n"/>
    <admst:for-each select="lhs/probe">
      <admst:variable name="pprobe" select="%(.)"/>
      <admst:if test="$lhs/ddxprobe/branch/pnode[.=$pprobe/branch/pnode or .=$pprobe/branch/nnode]">
        <admst:variable name="ddxinsidederivate" select="yes"/>
      </admst:if>
    </admst:for-each>
    <admst:text format="// rhsf %(rhs/function) \n"/>
    <admst:apply-templates select="rhs/function" match="function:precomputation"/>
    <admst:text format="/// rhsf %(rhs/function)\n"/>
    <admst:if test="lhs[insource='yes']">
      <admst:text format="// fixme: 7581 dprec %(rhs[not(nilled(function))])\n"/>
      <admst:apply-templates select="rhs[not(nilled(function))]" match="function:derivate:precomputation"/>
    </admst:if>
    <admst:choose>
      <admst:when test="rhs/tree/adms[datatypename='function']/..[name='ddx']">
        <admst:text format="// #if defined(_DDX) line? \n"/>
        <admst:if test="[$_DDX='yes']">
          <admst:variable name="ddxprobe" select="%(rhs/tree/arguments[2])"/>
          <admst:text test="lhs[insource='yes']/probe" format="// #if defined(_DERIVATE)\n"/>
          <admst:if test="[$_DERIVATE='yes']">
            <admst:for-each select="lhs[insource='yes']/probe">
              <admst:variable name="qprobe" select="%(.)"/>
              <admst:variable name="allepq"/>
              <admst:for-each select="$lhs/probe">
                <admst:variable name="pprobe" select="%(.)"/>
                <admst:choose>
                  <admst:when test="$pprobe/branch/pnode[.=$ddxprobe/branch/pnode]">
                    <admst:apply-templates select="$rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
                    <admst:variable name="allepq" select="$allepq+($epq)"/>
                  </admst:when>
                  <admst:when test="$pprobe/branch/nnode/[.=$ddxprobe/branch/pnode]">
                    <admst:apply-templates select="$rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
                    <admst:variable name="allepq" select="$allepq-($epq)"/>
                  </admst:when>
                </admst:choose>
              </admst:for-each>
              <admst:variable name="pprobe" select="%($qprobe)"/>
              <admst:if test="[$requiredderivateforddx='yes']">
                <admst:text format="// rdddx\n%(dxname($lhs)/[name='dxname']/value)=$allepq;\n"/>
              </admst:if>
            </admst:for-each>
          </admst:if>
        </admst:if>
		  <admst:text test="lhs[insource='yes']/probe" format="// #endif DERIVATE \n"/>
        <admst:variable name="allep"/>
        <admst:variable name="qprobe"/>
        <admst:for-each select="$lhs/probe">
          <admst:variable name="pprobe" select="%(.)"/>
          <admst:choose>
            <admst:when test="$pprobe/branch/pnode[.=$ddxprobe/branch/pnode]">
              <admst:text format="// rns? %(adms/datatypename)\n"/>
              <admst:apply-templates select="$rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
              <admst:variable name="allep" select="$allep+($ep)"/>
            </admst:when>
            <admst:when test="$pprobe/branch/nnode/[.=$ddxprobe/branch/pnode]">
              <admst:apply-templates select="$rhs/tree/arguments[1]" match="%(adms/datatypename)"/>
              <admst:variable name="allep" select="$allep-($ep)"/>
            </admst:when>
          </admst:choose>
        </admst:for-each>
        <admst:text format="%(lhs/name) = $allep; //allep\n"/>
        <admst:text format="EXIT_IF_ISNAN(%(lhs/name)); //7622\n"/>
		  <admst:text format="// #endif _DDX\n"/>
      </admst:when>
      <admst:otherwise>
        <admst:if test="lhs[insource='yes']">
          <admst:variable name="definedrequired" select="yes"/>
          <admst:variable name="ddxreq" select="no"/>
          <admst:variable name="derreq" select="no"/>
          <admst:choose>
            <admst:when test="[$ddxinsidederivate='yes']">
              <!--   <admst:text format="#if defined(_DERIVATEFORDDX) \n
                // probe=%($lhs/probe) ddxprobe=%($lhs/ddxprobe) \n"/> -->
              <admst:variable name="ddxreq" select="yes"/>
            </admst:when>
            <admst:when test="lhs/probe">
                <!--     <admst:text format="#if defined(_DERIVATE)
                // FIXME E // probe=%($lhs/probe)  ddxprobe=%($lhs/ddxprobe) \n"/> -->
              <admst:variable name="dderreq" select="yes"/>
            </admst:when>
            <admst:otherwise>
              <admst:variable name="definedrequired" select="no"/>
            </admst:otherwise>
          </admst:choose>

          <admst:if test="[$ddxreq='no' or $_DERIVATEFORDDX='yes' or (( $derreq='no' or $_DERIVATE='yes') and $ddxreq='no' )] " >
            <admst:for-each select="lhs/probe">
              <admst:variable name="pprobe" select="%(.)"/>
              <admst:variable name="ddxinsidethisprobe" select="no"/>
              <admst:if test="$lhs/ddxprobe/branch/pnode[.=$pprobe/branch/pnode or .=$pprobe/branch/nnode]">
                <admst:variable name="ddxinsidethisprobe" select="yes"/>
              </admst:if>
              <admst:variable name="isinside" select="0"/>
              <admst:if test="$rhs/probe[.=$pprobe]">
                <admst:variable name="isinside" select="1"/>
              </admst:if>
              <admst:variable name="qprobe"/>
              <admst:variable name="ep" select="0.0"/>
              <admst:apply-templates select="[$isinside='1']/$rhs/tree" match="%(adms/datatypename)"/>
              <!-- <admst:text test="[$ddxinsidederivate='yes' and $ddxinsidethisprobe='no']"
                   format="#if defined(_DERIVATE) // FIXME 8\n"/> -->
              <admst:if  test="[$ddxinsidederivate='no' and $ddxinsidethisprobe='yes' or $_DERIVATE='yes']" >
            <admst:text format="%(dxname($lhs)/[name='dxname']/value) = /*ep*/ $ep; // 7672"/>
		assert(is_number(%(dxname($lhs)/[name='dxname']/value)));

              </admst:if>
              <!-- <admst:text test="[$ddxinsidederivate='yes' and $ddxinsidethisprobe='no']"
                   format="#endif //3\n"/>
              <admst:text test="[$ddxinsidethisprobe='yes']" format="#if defined(_DERIVATE2)\n"/> -->
              <admst:if test="[$ddxinsidethisprobe='no' or $_DERIVATE2='yes' ]">
                <admst:for-each select="$lhs[$ddxinsidethisprobe='yes']/probe">
                  <admst:variable name="epq" select="0."/>
                  <admst:variable name="qprobe" select="%(.)"/>
                  <admst:apply-templates select="[$isinside='1']/$rhs/tree" match="%(adms/datatypename)"/>
                  <admst:text format="  %(ddxname($lhs)/[name='ddxname']/value) = $epq; // epq\n"/>
                </admst:for-each>
                <!-- <admst:text test="$lhs[$ddxinsidethisprobe='yes']" format="#endif // 5\n"/> -->
              </admst:if>
            </admst:for-each>
            <!-- <admst:text test="[$definedrequired='yes']" format="#endif // _DERIVATE \n"/> -->
          </admst:if>
          <admst:variable name="ddxinsidederivate" select="no"/>
        </admst:if>
        <admst:variable name="qprobe"/>
        <admst:apply-templates select="lhs" match="variable"/>


          <admst:if test="[$eval_kept='yes']">
            <admst:text format="$e = "/> %(e(rhs/tree)) <admst:text format=" /*tree, incomplete.*/ \n"/>
          </admst:if>
          <admst:if test="[$eval_kept='no']">
          <admst:text format="double x%(adms/id(.)) = $e = "/> %(e(rhs/tree)) ;
		assert(is_number(x%(adms/id(.)))); /*tree7704*/
          </admst:if>
        <!--<admst:text format="EXIT_IF_ISNAN($e)\n"/>-->
      </admst:otherwise>
    </admst:choose>
    <admst:text test="rhs/function" format="} //rhs \n"/>
    <!--   <admst:text test="[dynamic='yes']" format="#endif // _DYNAMIC \n"/> -->
  </admst:if>
  <admst:text format="// FIXME done here\n"/>
</admst:template>
<!-- ---------------------------------------------------------- -->
<admst:template match="contribution:nonoise:mint">
  <admst:message format="*(tr_eval)contribution:nonoise:mint in use?*\n"/>
  <admst:text test="[dynamic='yes']" format="#if defined(_DYNAMIC) // fixm4\n"/>
  <admst:text test="rhs/function" format="{\n"/>
  <admst:apply-templates select="rhs/function" match="function:precomputation"/>
  <admst:apply-templates select="rhs[not(nilled(function))]" match="function:derivate:precomputation"/>
  <admst:variable name="sourcepnode" select="%(lhs/branch/pnode)"/>
  <admst:variable name="sourcennode" select="%(lhs/branch/nnode)"/>
  <admst:variable name="sourcepnodename" select="%($sourcepnode/name)"/>
  <admst:variable name="sourcennodename" select="%($sourcennode/name)"/>
  <admst:choose>
    <admst:when test="[dynamic='yes']">
      <admst:variable name="jname" select="dQ_dV"/>
      <admst:choose>
        <admst:when test="$sourcennode[grounded='no']">
          <admst:text format="charges_$(sourcepnodename)_$(sourcennodename) += %(e(rhs/tree));\n"/>
        </admst:when>
        <admst:otherwise>
          <admst:text format="charges_$(sourcepnodename)_$(sourcepnodename) += %(e(rhs/tree));\n"/>
        </admst:otherwise>
      </admst:choose>
    </admst:when>
    <admst:otherwise>
      <admst:variable name="jname" select="dI_dV"/>
      <admst:choose>
        <admst:when test="$sourcennode[grounded='no']">
          <admst:text format="currents[$sourcepnodename*TotalNodes+$sourcennodename] += %(e(rhs/tree));\n"/>
        </admst:when>
        <admst:otherwise>
          <admst:text format="currents[$sourcepnodename*TotalNodes+$sourcepnodename] += %(e(rhs/tree));\n"/>
        </admst:otherwise>
      </admst:choose>
    </admst:otherwise>
  </admst:choose>
  <admst:text format="#if defined(_DERIVATE) // FOO\n"/>
  <admst:for-each select="rhs/probe">
    <admst:variable name="probepnode" select="%(branch/pnode)"/>
    <admst:variable name="probennode" select="%(branch/nnode)"/>
    <admst:variable name="probepnodename" select="%($probepnode/name)"/>
    <admst:variable name="probennodename" select="%($probennode/name)"/>
    <admst:variable name="pprobe" select="%(.)"/>
    <admst:apply-templates select="../tree" match="%(adms/datatypename)"/>
    <admst:choose>
      <admst:when test="$probennode[grounded='no']">
        <admst:if test="$sourcennode[grounded='no']">
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probepnodename]+=$ep;\n"/>
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probennodename]-=$ep;\n"/>
          <admst:text format="  $jname[$sourcennodename*TotalNodes+$probepnodename]-=$ep;\n"/>
          <admst:text format="  $jname[$sourcennodename*TotalNodes+$probennodename]+=$ep;\n"/>
        </admst:if>
        <admst:if test="$sourcennode[grounded='yes']">
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probepnodename]+=$ep;\n"/>
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probennodename]-=$ep;\n"/>
        </admst:if>
      </admst:when>
      <admst:otherwise>
        <admst:if test="$sourcennode[grounded='no']">
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probepnodename]+=$ep;\n"/>
          <admst:text format="  $jname[$sourcennodename*TotalNodes+$probepnodename]-=$ep;\n"/>
        </admst:if>
        <admst:if test="$sourcennode[grounded='yes']">
          <admst:text format="  $jname[$sourcepnodename*TotalNodes+$probepnodename]+=$ep;\n"/>
        </admst:if>
      </admst:otherwise>
    </admst:choose>
  </admst:for-each>
  <admst:text format="#endif //88\n"/>
  <admst:text test="rhs/function" format="}\n"/>
  <admst:text select="[dynamic='yes']" format="#endif //98\n"/>
</admst:template>
<!-- ---------------------------------------------------------- -->
<admst:template match="contribution:nonoise">
  <admst:text format="// contribution:nonoise %(.) \n"/>

  <admst:if test="rhs[not(nilled(function[class='builtin']))]">
	{ // %(.) 7804

    <admst:for-each select="rhs/function">
      <admst:value-of select="position(.)-1"/>
		<admst:apply-templates select="." match="function:getname_push"/>
		<admst:text format="double __%s_%s = 0.0; // forward declaration...?\n"/>
    </admst:for-each>
// calling ddx:function:computation %(.)

    <admst:apply-templates select="." match="ddx:function:computation"/>
// done ddxfunctioncomp
  </admst:if>

  <admst:apply-templates select="." match="residual_"/>

  <admst:for-each select="rhs/probe">
    <admst:text format="\t// Rhs/probe: %(.)\n"/>
    <admst:text format="\t// dynamic: %(../dynamic)\n"/>
    <admst:variable name="probepnode" select="%(branch/pnode)"/>
    <admst:variable name="probennode" select="%(branch/nnode)"/>
    <admst:variable name="probepnodename" select="%($probepnode/name)"/>
    <admst:variable name="probennodename" select="%($probennode/name)"/>
    <admst:variable name="pprobe" select="%(.)"/>

	 <admst:apply-templates select="../tree" match="subexpression:differentiate"/>

	 // e: $e
	 // ep: $ep

    <admst:if test="branch/pnode[grounded='no']"> <!-- Pp != GND, but why?-->
      <admst:if test="../../lhs/branch/pnode[grounded='no']"> <!-- Sp != GND, but why? -->
        <admst:apply-templates select="." match="jacobian_"/>
      </admst:if>
    </admst:if>

	  <admst:if test="branch/pnode[grounded!='no']">
      <admst:warning format="untested: excluding Pp == GND"/>
    </admst:if>
	  <admst:if test="[../../lhs/branch/pnode[grounded!='no']]">
      <admst:warning format="untested: excluding Sp == GND"/>
    </admst:if>
	</admst:for-each>

  <admst:if test="rhs[not(nilled(function[class='builtin']))]">
	}
  </admst:if>

  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)contribution:nonoise<-*\n"/>
</admst:template>

<!-- ----------------------  ic vs  ------------------ -->
<admst:template match="ic_vs">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)ic_vs*\n"/>
  <admst:variable name="Pp" select="%(lhs/branch/pnode/name)"/>
  <admst:variable name="Pn" select="%(lhs/branch/nnode/name)"/>
	{ // ic_vs
	double vic = 0; // old: ic_$(Pp)_$(Pn);
		assert(is_number(vic));
		assert( is_number(it0[$(node_prefix)$(Pp)] ));
		assert( is_number(it0[$(node_prefix)$(Pn)] ));
		it0[$(node_prefix)$(Pp)] -= vic/OPT::shortckt;
		it0[$(node_prefix)$(Pn)] += vic/OPT::shortckt;
		assert( is_number(it0[$(node_prefix)$(Pp)] ));
		assert( is_number(it0[$(node_prefix)$(Pn)] ));

		<admst:text test="[$eval_mode='tr']" format="\t\t_write_ptr($Pp, $Pp, 1/OPT::shortckt);\n"/>
		<admst:text test="[$eval_mode='tr']" format="\t\t_write_ptr($Pn, $Pn, 1/OPT::shortckt);\n"/>
		<admst:text test="[$eval_mode='tr']" format="\t\t_write_ptr($Pp, $Pn, -1/OPT::shortckt);\n"/>
		<admst:text test="[$eval_mode='tr']" format="\t\t_write_ptr($Pn, $Pp, -1/OPT::shortckt);\n"/>
	}

</admst:template>
<!-- ----------------------  residuals ------------------ -->
<admst:template match="residual_">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)residual2*\n"/>
  <admst:variable name="Pp" select="%(lhs/branch/pnode/name)"/>
  <admst:variable name="Pn" select="%(lhs/branch/nnode/name)"/>
	// residual2($(Pp),$(Pn),%(rhs/tree))
	{  // recheck: static is same as dynamic?
		double vx = %(e(rhs/tree)); // %(rhs/tree)
		assert(is_number(vx));
		it0[$(node_prefix)$(Pp)] -= vx;
	<admst:if test="lhs/branch/nnode[grounded='no']">
		it0[$(node_prefix)$(Pn)] += vx;
	</admst:if>
	}

</admst:template>

<admst:template match="residual1">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)residual1*\n"/>
  <admst:variable name="Pp" select="%(lhs/branch/pnode/name)"/>
	// residual1($(Pp),%(rhs/tree))
	{ // recheck: static is same as dynamic?
		double vx = %(e(rhs/tree));
		assert(is_number(vx));
		it0[$node_prefix$Pp] -= vx;
	}

</admst:template>

<!-- ----------------------  jacobian stamp ------------------ -->
<admst:template match="jacobian_">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)jacobian_*\n"/>
  <admst:variable name="Sp" select="%(../../lhs/branch/pnode/name)"/>
  <admst:variable name="Sn" select="%(../../lhs/branch/nnode/name)"/>
  <admst:variable name="Pp" select="%(branch/pnode/name)"/>
  <admst:variable name="Pn" select="%(branch/nnode/name)"/>

  <admst:choose>

		<admst:when test="..[dynamic='yes']">
#ifndef NOJAC
	{ // jacobian4($(Sp), $(Sn), $(Pp), $(Pn)) - dynamic
		double eac = $ep;
		double etr = eac * CKTag0;
		double ceq;

#ifdef HAVE_KEEPCOEFF
		if(uic_now( _ddt_$pddt ) /* HACK. argh.*/){
			ceq = eac * BR($Pp, $Pn);

			keep_$pddt=1;	
			eac = $ep;
			etr = - eac * OPT::keepcoeff;
		} else
#endif
		{
			ceq = etr * BR($Pp, $Pn);
		}

<admst:if test="../../lhs/branch/pnode[grounded='no']">
		it0[$(node_prefix)$Sp] += ceq;
</admst:if>
<admst:if test="../../lhs/branch/nnode[grounded='no']">
		it0[$(node_prefix)$Sn] -= ceq;
</admst:if>

		_write_d$(eval_mode)($Sp, $Pp, e$(eval_mode)); // 1
<admst:if test="../../lhs/branch/nnode[grounded='no']">
	<admst:if test="branch/nnode[grounded='no']">
		_write_d$(eval_mode)($Sn, $Pn, e$(eval_mode)); // 2
	</admst:if>
		_write_d$(eval_mode)($Sn, $Pp, -e$(eval_mode)); // 3
</admst:if>
<admst:if test="branch/nnode[grounded='no']">
		_write_d$(eval_mode)($Sp, $Pn, -e$(eval_mode)); // 4
</admst:if>
	}
#endif

    </admst:when>

    <admst:otherwise>
	{ // jacobian4($(Sp), $(Sn), $(Pp), $(Pn), $(ep)) - static
		double ep = $ep;
		it0[$(cont_prefix)$Sp] += ep*BR($Pp,$Pn);
<admst:if test="../../lhs/branch/nnode[grounded='no']">
		it0[$(cont_prefix)$Sn] -= ep*BR($Pp,$Pn);
</admst:if>

		_write_p$(eval_mode)($Sp, $Pp, ep); // 1
<admst:if test="../../lhs/branch/nnode[grounded='no']">
	<admst:if test="branch/nnode[grounded='no']">
		_write_p$(eval_mode)($Sn, $Pn, ep); // 2
	</admst:if>
		_write_p$(eval_mode)($Sn, $Pp, -ep); // 3
</admst:if>
<admst:if test="branch/nnode[grounded='no']">
		_write_p$(eval_mode)($Sp, $Pn, -ep); // 4
</admst:if>
	}
    </admst:otherwise>

  </admst:choose>
</admst:template>

<!-- ---------------------------------------------------------- -->
<!-- analog//conditional -->
<admst:template match="conditional">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)conditional*\n"/>
  <admst:text select="if[dynamic='yes']" format="#ifdef _DYNAMIC // BUG182\n"/>
  <admst:text test="if/function" format="{ // conditional\n"/>
  <admst:apply-templates select="if/function" match="function:precomputation"/>
  <admst:text format="if\n(%(e(if/tree))) //conditional \n"/>
  <admst:text select="then/adms[datatypename!='block']" format="{\n"/>
  <admst:apply-templates select="then" match="%(adms/datatypename)"/>
  <admst:text select="then/adms[datatypename!='block']" format="}//conditional\n"/>
  <admst:if test="else">
    <admst:text format="else // conditional\n"/>
    <admst:text test="else/adms[datatypename!='block']" format="{ // conditional.\n"/>
    <admst:apply-templates select="else" match="%(adms/datatypename)"/>
    <admst:text test="else/adms[datatypename!='block']" format="} // conditional.\n"/>
  </admst:if>
  <admst:text test="if/function" format="} // conditional\n"/>
  <admst:text select="if[dynamic='yes']" format="#endif /* if(...) */\n"/>
</admst:template>

<!-- ---------------------------------------------------------- -->
<!-- analog//nilled -->
<admst:template match="nilled">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)nilled*\n"/>
</admst:template>

<!-- ---------------------------------------------------------- -->
<!-- analog//whileloop -->
<admst:template match="whileloop">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)whileloop*\n"/>
  <admst:text select="while[dynamic='yes']" format="#ifdef _DYNAMIC\n"/>
  <admst:variable name="SkipFVariable" select="y"/>
  <admst:text format="while\n(%(e(while/tree)))\n"/>
  <admst:variable name="SkipFVariable" select="n"/>
  <admst:text select="whileblock/adms[datatypename!='block']" format="{\n"/>
  <admst:apply-templates select="whileblock" match="%(adms/datatypename)"/>
  <admst:text select="whileblock/adms[datatypename!='block']" format="}\n"/>
  <admst:text select="while[dynamic='yes']" format="#endif /*&lt;/dynamic_while&gt;*/\n"/>
</admst:template>

<!-- ---------------------------------------------------------- -->
<!-- analog//callfunctions -->
<admst:template match="callfunction">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)callfunction*\n"/>
  <admst:choose>
    <admst:when test="function[name='\$strobe']">_strobe(</admst:when>
    <admst:when test="function[name='\$warning']">_warning(</admst:when>
    <admst:when test="function[name='\$error']">_error(</admst:when>
    <admst:when test="function[name='\$finish']">_finish(</admst:when>
    <admst:when test="function[name='\$stop']">_stop(</admst:when>
    <admst:otherwise>
      <admst:fatal format="function not supported: %(function/name)\n"/>
    </admst:otherwise>
  </admst:choose>
  <admst:join select="function/arguments" separator=",">
    %(e(tree))
    <admst:if test="[position(.)=1]">&quot;\\n&quot;</admst:if>
  </admst:join>
  <admst:text format=");\n"/>
</admst:template>

<!-- ---------------------------------------------------------- -->
<!-- analog/code -->
<!-- save all variables used for local declaration -->
<admst:variable name="ddxinsidederivate" select="no"/>

<!-- ---------------------------------------------------------- -->
<admst:template match="variable:declaration">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:declaration*\n"/>
  <admst:variable name="myvariable" select="%(.)"/>
  <admst:new datatype="list" arguments="list of ddx probes">
    <admst:variable name="ddxprobes" select="%(.)"/>
    <admst:for-each select="$myvariable/probe">
      <admst:variable name="pprobe" select="%(.)"/>
      <admst:push into="$ddxprobes/item"
        select="$myvariable/ddxprobe/branch/pnode[.=$pprobe/branch/pnode or .=$pprobe/branch/nnode]/$pprobe"
        onduplicate="ignore"/>
    </admst:for-each>
  </admst:new>
  <admst:if test="block/adms[datatypename='module']">
    <admst:text test="[static='no' and dynamic='yes']" format="// #if defined(_DYNAMIC) // fixme\n"/>
    <admst:if   test="[static='yes' or dynamic='no' or $_DYNAMIC='yes']">

    <admst:text test="[scope='local']" format="%(vtype(.)) %(name);\n"/>
    <admst:if test="[insource='yes']/probe">
      <admst:for-each select="probe">
        <admst:variable name="pprobe" select="%(.)"/>
        <admst:if test="../ddxprobe/branch/pnode[.=$pprobe/branch/pnode or .=$pprobe/branch/nnode]">
          <admst:variable name="ddxinsidederivate" select="yes"/>
        </admst:if>
      </admst:for-each>
      <!--       <admst:choose>
        <admst:when test="[$ddxinsidederivate='yes']">
          <admst:text format="#if defined(_DERIVATEFORDDX)\n"/>
        </admst:when>
        <admst:otherwise>
          <admst:text format="#if defined(_DERIVATE) //FIXME A /n"/>
        </admst:otherwise>
      </admst:choose> -->
      <admst:if test="[ $ddxinsidederivate='no' or $_DERIVATEFORDDX='yes' or ( $_DERIVATEFORDDX='no' and  $_DERIVATE='yes'   )    ]" >

        <admst:text select="probe"
          format="double %(../name)_%(nature/access)%(branch/pnode/name)_%(branch/nnode/name);\n"/>
        <!-- <admst:text test="[$ddxinsidederivate='yes']" format="#if defined(_DERIVATE2)\n"/> -->
        <admst:if test="[ $ddxinsidederivate='no'  or $_DERIVATE2='yes']" >
          <admst:for-each select="$ddxprobes/item">
            <admst:variable name="pprobe" select="%(.)"/>
            <admst:for-each select="$myvariable/probe">
              <admst:variable name="qprobe" select="%(.)"/>
              <admst:text format="  double %(ddxname($myvariable)/[name='ddxname']/value);\n"/>
            </admst:for-each>
          </admst:for-each>
        </admst:if>
      <!-- <admst:text test="[$ddxinsidederivate='yes']" format="#endif\n"/> -->
      <admst:variable name="ddxinsidederivate" select="no"/>
      <!-- <admst:text format="#endif // FIX 5\n"/> -->
      </admst:if>
    </admst:if>
    </admst:if>
    <admst:text test="[static='no' and dynamic='yes']" format="// #endif _DYNAMIC\n"/>
  </admst:if>
</admst:template>
<admst:template match="block:local:declaration">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)block:local:declaration*\n"/>
  <admst:choose>
    <admst:when test="adms[datatypename='assignment']">
      <admst:push into="module/evaluation/variable" select="lhs" onduplicate="ignore"/>
    </admst:when>
    <admst:when test="adms[datatypename='block']">
      <admst:apply-templates select="item" match="block:local:declaration"/>
    </admst:when>
    <admst:when test="adms[datatypename='conditional']">
      <admst:apply-templates select="then" match="block:local:declaration"/>
      <admst:apply-templates select="else" match="block:local:declaration"/>
    </admst:when>
    <admst:when test="adms[datatypename='whileloop']">
      <admst:apply-templates select="whileblock" match="block:local:declaration"/>
    </admst:when>
    <admst:when test="adms[datatypename='contribution']"/>
    <admst:when test="adms[datatypename='callfunction']"/>
    <admst:when test="adms[datatypename='nilled']"/>
    <admst:when test="adms[datatypename='blockvariable']"/>
    <admst:otherwise>
      <admst:fatal format="'datatypename=%(adms/datatypename)': should not be reached %s\n"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>
<!-- -------------------------------------------------- -->
<admst:template match="analog:evaluate">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)analog:evaluate*\n"/>
  <admst:text select="$fnoise/item" format="fpnoise%(index($fnoise/item,.))=0.0;
    fenoise%(index($fnoise/item,.))=0.0;\n"/>
  <admst:text select="$tnoise/item" format="tnoise%(index($tnoise/item,.))=0.0;\n"/>
  <admst:text select="$wnoise/item" format="wnoise%(index($wnoise/item,.))=0.0;\n"/>

  <!-- local variable declarations. move to template? -->
  <admst:for-each select="item">
    <admst:if test="adms[datatypename='block']/..[name!='initial_model' and name!='initial_instance']">
      <admst:text format="// declaration1\n"/>
      <admst:apply-templates select="." match="block:local:declaration"/>
    </admst:if>
    <admst:if test="adms[datatypename!='block']">
      <admst:text format="// declaration2\n"/>
      <admst:apply-templates select="." match="block:local:declaration"/>
    </admst:if>
  </admst:for-each>
  <admst:text format="// declaration3\n"/>
  <admst:apply-templates select="module/evaluation/variable" match="variable:declaration"/>
  <!-- local variable declarations -->

  <admst:reset select="module/evaluation/variable"/>
  <admst:for-each select="item">
    <admst:choose>
      <admst:when test="adms[datatypename!='block']">
        <admst:apply-templates select="." match="%(adms/datatypename)"/>
      </admst:when>
      <admst:otherwise>
        <admst:if test="[name!='initial_model' and name!='initial_instance']">
//bug -- parts of this (%(name), %(datatypename)) has been done in tr_begin.

          <admst:apply-templates select="." match="block"/>
          <admst:text format="// end bug?\n"/>
        </admst:if>
      </admst:otherwise>
    </admst:choose>
  </admst:for-each>
</admst:template>
<!-- -------------------------------------------------------------------- -->
<admst:template match="analog:initializeModel">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)analog:initializeModel*\n"/>
  <admst:for-each select="item">
    <admst:if test="adms[datatypename='block']/..[name='initial_model']">
      <admst:apply-templates select="." match="block:local:declaration"/>
    </admst:if>
  </admst:for-each>
  <admst:apply-templates select="module/evaluation/variable" match="variable:declaration"/>
  <admst:reset select="module/evaluation/variable"/>
  <admst:for-each select="item">
    <admst:if test="adms[datatypename='block']/..[name='initial_model']">
      <admst:apply-templates select="." match="block"/>
    </admst:if>
  </admst:for-each>
</admst:template>
<!-- -------------------------------------------------------------------- -->
<admst:template match="analog:initializeInstance">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)analog:initializeInstance*\n"/>
  <admst:for-each select="item">
    <admst:if test="adms[datatypename='block']/..[name='initial_instance']">
      <admst:apply-templates select="." match="block:local:declaration"/>
    </admst:if>
  </admst:for-each>
  <admst:apply-templates select="module/evaluation/variable" match="variable:declaration"/>
  <admst:reset select="module/evaluation/variable"/>
  <admst:for-each select="item">
    <admst:if test="adms[datatypename='block']/..[name='initial_instance']">
      <admst:apply-templates select="." match="block"/>
    </admst:if>
  </admst:for-each>
</admst:template>

<!--
* This template returns the description of an instance or
* a model parameter. It works for both formats of :
*   `ATTR(info="description"  ...)
* or
*   `ATTR(desc="description"  ...)
* This template is used in the creation of the mint:defineParameters
* routine. If there is no description given in the VerilogA file,
* it returns a NULL. The returned value is used in the 'description'
* field of the mint_param_xxx file, where xxx is [integer|real|string]
-->
<admst:template match="variable:desc">
  <admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:desc*\n"/>
  <admst:choose>
    <admst:when test="attribute[name='desc' or name='info']">
      <admst:return name="variable:desc" value="&quot;%(attribute[name='desc' or name='info']/value)&quot;"/>
    </admst:when>
    <admst:otherwise>
      <admst:return name="variable:desc" value="NULL"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!--
* This template returns the unit of an instance or
* a model parameter as given in the VA file.
* This template is used in the creation of the mint:defineParameters
* routine. If there is no 'unit' string given in the VerilogA file,
* it returns a NULL. The returned value is used in the 'unit'
* field of the mint_param_xxx file, where xxx is [integer|real|string]
-->
<admst:template match="variable:unit">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:unit*\n"/>
  <admst:choose>
    <admst:when test="attribute[name='unit']">
      <admst:return name="variable:unit" value="&quot;%(attribute[name='unit']/value)&quot;"/>
    </admst:when>
    <admst:otherwise>
      <admst:return name="variable:unit" value="NULL"/>
    </admst:otherwise>
  </admst:choose>
</admst:template>

<!--
* This template returns the default value of an instance or
* a model parameter as given in the VA file.
* This template is used in the creation of the mint:defineParameters
* routine. The returned value is used in the 'default'
* field of the mint_param_xxx file, where xxx is [integer|real|string]
-->

<admst:template match="variable:paramdef">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)variable:paramdef*\n"/>
  <admst:return name="name" value="%(name)"/>
  <admst:apply-templates select="." match="variable:desc">
    <admst:return name="desc" value="%(returned('variable:desc')/value)"/>
  </admst:apply-templates>
  <admst:apply-templates select="." match="variable:unit">
    <admst:return name="unit" value="%(returned('variable:unit')/value)"/>
  </admst:apply-templates>
  <admst:return name="minttype" value="%(type)"/>
</admst:template>

<!-- ------------------------------------ -->
<admst:template match="c:math_h">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)c:math_h*\n"/>
// *** arithmetics ***
<admst:text format="\n\n"/>
</admst:template>
<!-- ------------------------------------------- -->
<admst:template match="wrapper">
<admst:message test="[/dbg_xml='yes']" format="*(tr_eval)wrapper*\n"/>

<admst:if test="/module/variable[derivate='yes' and insource='yes']">
<admst:variable name="requiredderivateforddx" select="yes"/>
</admst:if>
<admst:text format="\n"/>
</admst:template>


<!-- 2 -->
</admst>