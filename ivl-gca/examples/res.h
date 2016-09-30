
#define abstol_potential_electrical 1e-6
#define abstol_flow_electrical 1e-12
#define abstol_potential_voltage 1e-6
#define abstol_flow_voltage 
#define abstol_potential_current 1e-12
#define abstol_flow_current 
#define abstol_potential_magnetic 1e-12
#define abstol_flow_magnetic 1e-9
#define abstol_potential_thermal 1e-4
#define abstol_flow_thermal 1e-9
#define abstol_potential_kinematic 1e-6
#define abstol_flow_kinematic 1e-6
#define abstol_potential_kinematic_v 1e-6
#define abstol_flow_kinematic_v 1e-6
#define abstol_potential_rotational 1e-6
#define abstol_flow_rotational 1e-6
#define abstol_potential_rotational_omega 1e-6
#define abstol_flow_rotational_omega 1e-6
#define abstol_potential_degradational 1e-9
#define abstol_flow_degradational 1e-9

#define VALUE_NAME "#"
#define TAIL_SIZE 1
#define _DDX
#define _DERIVATEFORDDX

#define SPICE_LETTER "\0"
#define STATIC_RHS 1

#define CONSTCtoK (273.15) // also in gnucap/constans.h P_CELSIUS0

// counted ddts 0

// there are nodes that work different... hmmm
// BUG: number of stamp nodes is topology dependent!

class DEV_RES;
class MODEL_RES : public MODEL_CARD {
	private:
		static std::map<string, PARA_BASE MODEL_RES::*> _param_dict;
		static std::map<string, PARA_BASE MODEL_RES::*> _param_dict_low;
	protected:
		explicit MODEL_RES(const MODEL_RES& p);	// for clone
	public:
		explicit MODEL_RES(const DEV_RES* p);	// for dispatcher
		~MODEL_RES();
	public: // override virtual
		MODEL_CARD* clone()const {return new MODEL_RES(*this);}
	public: // analog functions from va

	public:
		void precalc_first();
		void precalc_last();

	public: // type
		void set_dev_type(const std::string& nt);
		std::string dev_type()const	{ return _key;}

	public: // parameters
		int param_count() const;
		bool param_is_printable(int i) const;
		std::string param_name(int i) const;
		std::string param_name(int i, int j)const;
		std::string param_value(int i) const;
		void set_param_by_index(int, std::string&, int);
		void set_param_by_name(std::string Name, std::string Value);
		bool IsParamGiven(const char *name) const;
		PARAMETER<double> v_R; // yes global_model real 


	public: // global_model



	public:
		std::string _key;
		std::string _level;
};

/* ========================================================= */
/* ========================================================= */

class DEV_RES : public ADMS_BASE
{
	private:
		std::string _modelname;
		node_t _nodes[ 2 ];
		static std::map<string, PARA_BASE DEV_RES::*> _param_dict;
		static std::map<string, PARA_BASE DEV_RES::*> _param_dict_low;
	private:
		explicit DEV_RES(const DEV_RES& p);
	public:
		explicit DEV_RES();
		~DEV_RES();
	private:
	public: // global_instance input

	protected: // override virtual
		char id_letter()const	{untested();return SPICE_LETTER[0];}
		bool print_type_in_spice()const {return true;}
		std::string value_name()const {return VALUE_NAME;}
		uint_t max_nodes()const {return 2;}
		uint_t min_nodes()const {return 0;}
		uint_t matrix_nodes()const {return 2;}
		// uint_t matrix_nodes()const {return net_nodes() + int_nodes();}
		uint_t net_nodes()const {return _net_nodes;}
		uint_t int_nodes()const {return 0;}
		CARD* clone()const {return new DEV_RES(*this);}

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
		double tr_probe_num(const std::string&)const;

	public: // ac
		void ac_iwant_matrix();
		void ac_begin();
		void do_ac();
		void ac_eval(); // -> COMMON?
		void ac_load();
		COMPLEX ac_involts()const {unreachable();return NOT_VALID;}
		COMPLEX ac_amps()const {unreachable();return NOT_VALID;}
		XPROBE  ac_probe_ext(const std::string&)const {itested(); return XPROBE(NOT_VALID, mtNONE);}
		uint_t tail_size()const {return TAIL_SIZE;}

#ifdef HAVE_TT
	public: // tt
		void tt_advance();
#if HAVE_TT + 0 == 2
		void do_tt();
#else
		void stress_apply();
#endif
		void tt_accept();
		void tt_begin();
		void tr_stress_last(); // tt_stress??
#endif

	public: // misc
		void keep_ic();
		bool has_memory(){return false;}

	public: // ports
		std::string port_name(uint_t i)const {
			assert(i < max_nodes());
			return node_name(i);
		}
	protected: //nodes
		std::string node_name(uint_t i)const {
			static std::string names[] = {
				"p",
				"n"
			};
			assert(i < max_nodes()+int_nodes());
			return names[i];
		}

	private: // parameters
		void set_param_by_name(std::string Name, std::string Value);

	private:
		void init_nodemap();
		PARTITION<2>* _nodemap;
		void guesstopology();
		void collapse_nodes(unsigned a, unsigned b)
		{ untested();
			assert(_nodemap);
			error(bTRACE,"%s: collapse node %u and %u\n", long_label().c_str(), a, b);
			_nodemap->identify(a,b);
		}
	protected: // implementation
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

	private: // global_instance (ask?)

	public: // global_model... hmmm
		double v_R;


	private: // node numbers
		// stick to gnucap node ordering (external first)
		enum __Node {
			__p=0,
			__n=1,
			_num_nodes
		};
		// map nodenames to itl slots
		// some disciplines dont have itls
		enum __CNode {
			c_p=0,
			c_n=1
		};
	protected: // nodes sorted according to discipline

		enum Node_degradational {
		};
		enum Node_rotational_omega {
		};
		enum Node_rotational {
		};
		enum Node_kinematic_v {
		};
		enum Node_kinematic {
		};
		enum Node_thermal {
		};
		enum Node_magnetic {
		};
		enum Node_current {
		};
		enum Node_voltage {
		};
		 // REQUIRED, write out in any case (incomplete)...  electrical
		enum Node_electrical {
			n_GND=-1,
			n_p=__p /* electrical -- electrical */ ,
			n_n=__n /* electrical -- electrical */ 
		};

	private: // nodes / matrix indexes

// raw value
#define NR(p) _NR(n_ ## p)
inline double _NR(int p){ return _n[p].v0(); } 

#define BR(p,n) _BR(n_ ## p, n_ ## n)
inline double _BR(int p, int n) {
	if (n==-1)
		return _n[p].v0();
	return volts_limited(_n[p],_n[n]);
}

// discipline wrapped value
#define NP(p) _NP(n_ ## p)
inline voltage_t _NP(Node_electrical p){
	return std::min(1000., std::max(-1000., _n[p].v0()) );
}
inline double _NP(Node_magnetic p){
	assert(is_number(_n[p].v0()));
	return _n[p].v0();
}
inline double _NP(Node_rotational_omega p){
	assert(is_number(_n[p].v0()));
	return _n[p].v0();
}
inline double _NP(Node_thermal p){
	assert(is_number(_n[p].v0()));
	return _n[p].v0() + P_CELSIUS0;
}
#define BP(p,n) _BP(n_ ## p, n_ ## n)
inline voltage_t _BP(Node_electrical p, Node_electrical n) { return volts_limited(_n[p],_n[n]);}

#define _N(n) _n[n_ ## n]

	private:
		enum EMatrixEntries {
			m_p_p=0 /* electrical to electrical */,
			m_p_n=1 /* electrical to electrical */,
			m_n_p=2 /* electrical to electrical */,
			m_n_n=3 /* electrical to electrical */
		};
		enum ElossEntries { //?
		};

	enum EStates {
		_num_states //  = 0 // 0+0
	};

	double* _states[OPT::_keep_time_steps]; // for each step: q0,i0, q1,i1, q2,i2, etc.
	double* _states_q1;	// for step t1:  q0, q1, q2, q3, etc.
	double* _states_kept;

//even states: _y
// odd states: _i

	double _DDT(double qq, unsigned stateno){
		unsigned qcap = stateno*2;
		trace4("DEV_RES::_DDT", stateno, qq, _sim->_mode, _sim->_phase);
		_states[0][stateno*2] = qq;

		if (_sim->analysis_is_ac()) {
			// BUG? fix later.
			// shouldnt call DDT at all!
			return 0;
		}
		assert(_dt != 0. || _sim->analysis_is_restore());
		METHOD method;
		if( (_time[1]!=0) && (_method_a==mTRAP)) {
			method = mTRAP;
		} else {
			method = mEULER;
		}

		std::valarray<FPOLY1> q(OPT::_keep_time_steps);
		std::valarray<FPOLY1> i(OPT::_keep_time_steps);

		unsigned k = OPT::_keep_time_steps;
		// stupid copy. fixme.
		// gnucap stores things more usable.

		q[0] = FPOLY1(NOT_VALID, _states[0][qcap], 0.);

		for (unsigned ii = 0; ii < k; ii++)
		{
			trace5("DEV_RES::_DDT", ii, _states[ii][qcap], _states[ii][qcap+1], _time[0], _time[1] );
			assert( _states[ii][qcap] == _states[ii][qcap] );
			assert( _states[ii][qcap+1] == _states[ii][qcap+1] );
			q[ii] = FPOLY1(NOT_VALID, _states[ii][qcap], 0.);
			i[ii] = FPOLY1(NOT_VALID, _states[ii][qcap+1], 0.);
		}

		trace0("calling differentiate");
		i[0] = differentiate(&q[0], &i[0], _time, _method_a);
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
	double _IDT(double qq, unsigned stateno) {
		unsigned qcap = stateno * 2;
		trace4("DEV_RES::_DDT", stateno, qq, _sim->_mode, _sim->_phase);
		_states[0][qcap+1] = qq;

		if (_sim->analysis_is_ac()) {
		// BUG? fix later.
		// shouldnt call IDT at all!
			return 0;
		}
		assert(_dt != 0.0);
		METHOD method; USE(method);
		if( (_time[1]!=0) && (_method_a==mTRAP)) { // incomplete(); ... later
			method = mTRAP;
		} else { itested();
			method = mEULER;
		}

		double integral;

		if( _sim->analysis_is_static() ) { untested();
			integral = 0;
			return 0;
		} else { untested();
			integral = _states[1][qcap];
			integral += _dt * ( _states[1][qcap+1] + qq ) / 2;
		}

		_states[0][qcap] = integral;
		assert(integral == integral);
		return integral;
	}

	// noise: something is calculated, but never used

//		bool _need_accept;

	public: // initial conditions
	private: // discipline stuff
		// absolute checks FIXME: this is not the way gnucap does checks
		bool conchk( Node_electrical a, double reltol=OPT::reltol ) const {
			/// hmm really OPT::abstol?
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()], _n[a].v0(), OPT::abstol, reltol);
		}
		bool conchk( Node_magnetic a, double reltol=OPT::reltol ) const {
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()], _n[a].v0(), abstol_flow_magnetic, reltol);
		}
		bool conchk( Node_rotational_omega a, double reltol=OPT::reltol ) const {
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()], _n[a].v0(), abstol_flow_rotational_omega, reltol);
		}
		bool conchk( Node_thermal a, double reltol=OPT::reltol ) const {
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()], _n[a].v0(), abstol_flow_thermal, reltol);
		}
		bool conchk( Node_degradational a, double reltol=OPT::reltol ) const {
			return true;
		}
		bool conchk( Node_degradational, Node_degradational, double=0) const {
			return true;
		}
		// relative checks, FIXME: only where needed and keep history...
		bool conchk( Node_electrical a, Node_electrical b, double reltol=OPT::reltol ) const {
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()] - CKT_BASE::_sim->_vt1[_n[b].m_()],
			    _n[a].v0() - _n[b].v0(), OPT::abstol, reltol);
		}
		bool conchk( Node_thermal a, Node_thermal b, double reltol=OPT::reltol ) const {
			return ::conchk(CKT_BASE::_sim->_vt1[_n[a].m_()] - CKT_BASE::_sim->_vt1[_n[b].m_()],
			    _n[a].v0() - _n[b].v0(), abstol_flow_thermal, reltol);
		}
		// these are never inputs.
		bool conchk( Node_electrical, Node_degradational, double=0 ) const {
			return true;
		}
		bool conchk( Node_degradational, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_rotational_omega, double=0 ) const {
			return true;
		}
		bool conchk( Node_rotational_omega, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_rotational, double=0 ) const {
			return true;
		}
		bool conchk( Node_rotational, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_kinematic_v, double=0 ) const {
			return true;
		}
		bool conchk( Node_kinematic_v, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_kinematic, double=0 ) const {
			return true;
		}
		bool conchk( Node_kinematic, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_thermal, double=0 ) const {
			return true;
		}
		bool conchk( Node_thermal, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_magnetic, double=0 ) const {
			return true;
		}
		bool conchk( Node_magnetic, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_current, double=0 ) const {
			return true;
		}
		bool conchk( Node_current, Node_electrical, double=0 ) const {
			return true;
		}
		bool conchk( Node_electrical, Node_voltage, double=0 ) const {
			return true;
		}
		bool conchk( Node_voltage, Node_electrical, double=0 ) const {
			return true;
		}
// ----- netlist ------------ //
	private: // subdevices
	private: // branches
      // p,n
	private: // sources
		// I(p,n)
	private: // module_source_branches
		// p,n p n
	private: // module_probe_branches
		// p,n p n
	private: // probes
      // V(p,n)
	private: // all sources
// I(p,n)
	private: // unique (directed) branches
// p,n
	private: // sources by branch
		COMPONENT* _br_p_n;//I(p,n)


		double _br_p_n_xxx[
		    + 2 // constant and self admittance
		    + 0 // V transadmittance
		    + 2*0 // I transadmittance
      ];
		double _br_p_n_xxx_old[
		    + 2 // constant and self admittance
		    + 0 // V transadmittance
		    + 2*0 // I transadmittance
      ];

		enum{
			P__br_p_n__br_n_p_V=-1,
			P__br_p_n_const=0,
			P__br_p_n__br_p_n_V=1
		};

		enum{ // different polarity "aliases"
			P__br_n_p__br_p_n_V=-1,
			P__br_n_p__br_n_p_V=1
		/*.I.*/
		/*...*/
		};



#ifndef NDEBUG
		unsigned _loaditer;
#endif
	private: // ops
#ifndef NDEBUG
	private: // debugging
		int _ddtmindt_index;
#endif
}; // DEV_RES
// ------------------------------------------------------------------------ //
// ------------------------------------------------------------------------ //
class COMMON_RES : public COMMON_ADMS
{
	friend class DEV_RES;
  private:
//    static map<string, PARA_BASE COMMON_RES::*> _param_dict;
	public:
		void expand(const COMPONENT*);
		void precalc_first(const CARD_LIST*);
		void precalc_last(const CARD_LIST*);
	explicit COMMON_RES(const COMMON_RES& p);
	explicit COMMON_RES(int c=0);
	~COMMON_RES();
	bool operator==(const COMMON_COMPONENT&)const;
	COMMON_COMPONENT* clone()const {return new COMMON_RES(*this);}
	std::string name()const {itested();return "res";}
	void set_param_by_name(std::string name, std::string value);
	public: // initial conditions

};

//inline double DEV_RES::BP2( DEV_RES::Node_electrical p,  DEV_RES::Node_electrical n  ) const{
//	assert(is_number( _n[p].v0() - _n[n].v0() ));
//	return  min( 100. , max (  _n[p].v0() - _n[n].v0(), -100. ));
//}
// -------------------------------------------------------------------------- //

//#endif
// vim should work with any ts=sw


// -------------------------------------------------------------------------- //
// fixme: generate from disciplines.h
