/* adms base 2011 felix salfelder
 * This file is part of "Gnucap", the Gnu Circuit Analysis Package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 *------------------------------------------------------------------
 * Base class for adms models
 */
#ifndef E_ADMS_H
#define E_ADMS_H
#include <u_lang.h>
#include <e_node.h>
#include <m_cpoly.h>
#include <l_denoise.h>
#include <e_subckt.h>
#include <u_status.h>
#include <algorithm>
#include "gcuf_compat.h"
/*--------------------------------------------------------------------------*/
using namespace std;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class COMMON_ADMS : public COMMON_COMPONENT {
	public:
		COMMON_ADMS(int i) : COMMON_COMPONENT(i){}
	protected:
		void attach_model(const COMPONENT* d); // not const!
};
/*--------------------------------------------------------------------------*/
class ADMS_BASE : public COMPONENT {
	protected:
		explicit ADMS_BASE();
		explicit ADMS_BASE(const ADMS_BASE& p);

		~ADMS_BASE() {}

		virtual void	   store_values()=0;
		//void   reject_values()		{ _y0 = _y1;}
	public: // type
		bool	   skip_dev_type(CS&);
		string dev_type()const {return common()->modelname();} // from BASE_SUBCKT
		virtual bool	   print_type_in_spice()const {return false;}
	public:
		// virtual void precalc_last();
		virtual void tr_begin();
		virtual void tr_restore();
		virtual void dc_advance();
		virtual void tr_advance();
		virtual void tr_regress();
		virtual bool tr_needs_eval()const;
		virtual bool do_tr() { unreachable(); return false; }
			
		// TIME_PAIR tr_review();
		virtual void  tr_stress() {
			trace0( ("ADMS_BASE device " + short_label() + ": no stress").c_str() );
		} // calcul

		//void   map_nodes();
		virtual double   tr_probe_num(const std::string&)const;
	protected: // inline, below
		template <class T>
			T   dampdiff(T*, const T&);

		template <class T>
			void	   tr_load_extended(const node_t& no1, const node_t& no2,
					const node_t& ni1, const node_t& ni2,
					T* value, T* old_value);

		void	   ac_load_extended(const node_t& no1, const node_t& no2,
				const node_t& ni1, const node_t& ni2,
				COMPLEX value);

		void	   tr_load_source_point(node_t& no1, conductance_t* value, conductance_t* old_value);
		void	   ac_load_source_point(node_t& no1, COMPLEX new_value);

		template <class T>
			void	   tr_load_diagonal_point(const node_t& no1, T* value, T* old_value);

		void	   ac_load_diagonal_point(const node_t& no1, COMPLEX value);

		template <class T>
		void	   tr_load_point(const node_t& no1, const node_t& no2, T* value, T* old_value);
		void	   ac_load_point(const node_t& no1, const node_t& no2,
				COMPLEX value);

		bool	   conv_check()const;
		bool	   has_tr_eval()const;
		bool	   has_ac_eval()const;
		bool	   using_tr_eval()const;
		bool	   using_ac_eval()const;
	
	public:
		double   tr_review_trunc_error(const FPOLY1* q);
		double   tr_review_check_and_convert(double timestep);

		double   tr_outvolts()const	{return dn_diff(_n[OUT1].v0(), _n[OUT2].v0());}
		double   tr_outvolts_limited()const{return volts_limited(_n[OUT1],_n[OUT2]);}
		COMPLEX  ac_outvolts()const	{return _n[OUT1].vac() - _n[OUT2].vac();}

		int param_count()const {return (0 + COMPONENT::param_count());}
	protected:
		int      _loaditer;	// load iteration number
	private:
		node_t   _nodes[NODES_PER_BRANCH]; // nodes (0,1:out, 2,3:in)
	public: // commons
		COMPLEX  _ev;		// ac effective value (usually real)
		double   _dt;
		double   _time[OPT::_keep_time_steps];
	private: // from storage
		int	   order()const {
			const int o[] = {1, 1, 2, 1, 1};
			int ord = o[_method_a];
			assert(ord < OPT::_keep_time_steps);
			return ord;
		}
		double   error_factor()const {
			const double f[]={1./2., 1./2., 1./12., 1./6., 1./2.};
			return f[_method_a];
		}
	public: // used by commons
		method_t _method_u;	/* method to use for this part per user */
		METHOD   _method_a;	/* actual integration method (auto)	*/
	protected: // from storag
		static METHOD method_select[meNUM_METHODS][meNUM_METHODS];
	public: // from storag
		double   tr_c_to_g(double c, double g)const;

	protected: // callfunctions
		void do_strobe(const char* fmt, ... )
		{
			if(::status.control == 1){  // FIXME: if print_now?!
				va_list arg_ptr;
				va_start(arg_ptr,fmt);
				vfprintf(stdout,fmt,arg_ptr);
				va_end(arg_ptr);
				printf("\n");
			}
		}
		void do_warning(const char* fmt, ... )
		{
			if(::status.control == 1){  // FIXME: if print_now?!
				va_list arg_ptr;
				va_start(arg_ptr,fmt);
				vfprintf(stderr,fmt,arg_ptr);
				va_end(arg_ptr);
				fprintf(stderr,"\n");
			}
		}
		void do_error(const char* fmt, ... )
		{
			va_list arg_ptr;
			va_start(arg_ptr,fmt);
			vfprintf(stderr,fmt,arg_ptr);
			va_end(arg_ptr);
			fprintf(stderr,"\n");
			throw(Exception("va error"));
		}
		void do_stop()
		{
			throw(Exception("va stop"));
		}
};
/*--------------------------------------------------------------------------*/
void COMPONENT::set_port_by_index(uint_t num, std::string& ext_name)
{
	if (num <= max_nodes()) {

		trace3("COMPONENT::set_port_by_index ", short_label(), ext_name, num);
		_n[num].new_node(ext_name, this);
		_net_nodes = max(_net_nodes, num+1);

	}else{untested();
		throw Exception_Too_Many(num, max_nodes(), 0/*offset*/);
	}
}
/*--------------------------------------------------------------------------*/
template <class T>
inline T ADMS_BASE::dampdiff(T* v0, const T& v1)
{
	//double diff = v0 - v1;
	assert(v0);
	if(*v0 != *v0){ 
		error(bDANGER,"ADMS_BASE::dampdiff: %s\n", short_label().c_str());
		//    std::cerr << "ungood: " << *v0 << "\n";
		assert(false);
	}
	assert(v1 == v1);
	T diff = dn_diff(*v0, v1);
	assert(diff == diff);
	if (!_sim->is_advance_or_first_iteration()) {
		diff *= _sim->_damp;
		*v0 = v1 + diff;
	}else{
	}
	return mfactor() * ((_sim->is_inc_mode()) ? diff : *v0);
}
/*--------------------------------------------------------------------------*/
template <class T>
inline void ADMS_BASE::tr_load_extended(const node_t& no1, const node_t& no2,
		const node_t& ni1, const node_t& ni2,
		T* new_value, T* old_value)
{
	// untested();
	assert( new_value == new_value);
	T d = dampdiff(new_value, *old_value);
	if (d != 0.) {
		_sim->_aa.load_asymmetric(no1.m_(), no2.m_(), ni1.m_(), ni2.m_(), d);
	}else{
	}
	*old_value = *new_value;
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::ac_load_extended(const node_t& no1, const node_t& no2,
		const node_t& ni1, const node_t& ni2,
		COMPLEX new_value)
{
	_sim->_acx.load_asymmetric(no1.m_(), no2.m_(), ni1.m_(), ni2.m_(), mfactor() * new_value);
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::tr_load_source_point(node_t& no1,
		conductance_t* new_value, conductance_t* old_value)
{
	assert( new_value == new_value);
	double d = dampdiff(new_value, *old_value);
	if (d != 0.) {
		if (no1.m_() != 0) {
			no1.i() += d;
		}else{
		}
	}else{
	}
	*old_value = *new_value;
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::ac_load_source_point(node_t& no1, COMPLEX new_value)
{
	if (no1.m_() != 0) {
		no1->iac() += mfactor() * new_value;
	}else{itested();
	}
}
/*--------------------------------------------------------------------------*/
template <class T>
inline void ADMS_BASE::tr_load_diagonal_point(const node_t& no1,
		T* new_value, T* old_value)
{
	T d = dampdiff(new_value, *old_value);
	if (d != 0.) {
		_sim->_aa.load_diagonal_point(no1.m_(), d);
	}else{
	}
	*old_value = *new_value;
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::ac_load_diagonal_point(const node_t& no1, COMPLEX new_value)
{
	_sim->_acx.load_diagonal_point(no1.m_(), mfactor() * new_value);
}
/*--------------------------------------------------------------------------*/
template <class T>
inline void ADMS_BASE::tr_load_point(const node_t& no1, const node_t& no2,
		T* new_value, T* old_value)
{
	T d = dampdiff(new_value, *old_value);
	if (d != 0.) {
		_sim->_aa.load_point(no1.m_(), no2.m_(), d);
	}else{
	}
	*old_value = *new_value;
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::ac_load_point(const node_t& no1, const node_t& no2,
		COMPLEX new_value)
{itested();
	_sim->_acx.load_point(no1.m_(), no2.m_(), mfactor() * new_value);
}
/*--------------------------------------------------------------------------*/
inline bool ADMS_BASE::conv_check()const
{
	incomplete();
	unreachable();
	return false;
}
/*--------------------------------------------------------------------------*/
inline bool ADMS_BASE::has_tr_eval()const
{
	//  untested0(long_label().c_str());
	return (has_common() && common()->has_tr_eval());
}
/*--------------------------------------------------------------------------*/
inline bool ADMS_BASE::has_ac_eval()const
{
	return (has_common() && common()->has_ac_eval());
}
/*--------------------------------------------------------------------------*/
inline bool ADMS_BASE::using_tr_eval()const
{
	return (has_probes() || has_tr_eval());
}
/*--------------------------------------------------------------------------*/
inline bool ADMS_BASE::using_ac_eval()const
{
	return (has_probes() || has_ac_eval());
}
/*--------------------------------------------------------------------------*/
inline void ADMS_BASE::tr_advance()
{
	trace2("ADMS_BASE::tr_advance ", _sim->_time0, net_nodes());
#ifdef USE_DTIME
	double delta = _time[0];
#else
	double delta = 0;
	assert(_time[0] < _sim->_time0); // moving forward
#endif

	for (int i=OPT::_keep_time_steps-1; i>0; --i) {
		assert(delta || _time[i] < _time[i-1] || _time[i] == 0.);
		_time[i] = _time[i-1] - delta;
	}

#ifdef USE_DTIME
	_time[0] = _sim->_dt0;
	_dt = _time[0];
#else
	_time[0] = _sim->_time0;
	_dt = _time[0] - _time[1];
#endif

#ifdef HAVE_DISCONT
	method_t m = _method_u;
	trace1("tr_advance", long_label());

	DISCONT d = disNONE;
	for (unsigned i=0;i < net_nodes(); ++i) {
		d |= _n[i]->discont();
	}

	if (_time[1] == 0) {
		m = meEULER;
	} else if (d) {
		m = meEULER;
	} else {
	}
	_method_a = method_select[OPT::method][m];
#endif

}
/*--------------------------------------------------------------------------*/
inline double ADMS_BASE::tr_probe_num(const std::string& x)const
{
	if (Umatch(x, "conv{erged} ")) {
		return converged();
	} else if (Umatch(x, "met{hod} ")) { untested();
		return _method_a;
	}

	return COMPONENT::tr_probe_num(x);
}
/*--------------------------------------------------------------------------*/
inline double ADMS_BASE::tr_c_to_g(double c, double g)const
{
	if (_sim->analysis_is_static()) {
		assert(_time[0] == 0.);
		return 0.;
	}else if (_sim->analysis_is_restore()) {itested();
		assert(_time[0] > 0);
		return g;
		// no change, fake
	}else{
		assert(_sim->analysis_is_tran_dynamic());
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
			default: incomplete();
		}
		return g;
	}
}
/*--------------------------------------------------------------------------*/
#define _write_ptr(S1, S2, value) \
	m_entries[m_## S1 ## _ ## S2 ]+=value;

#define _write_pac(S1, S2, value) \
	m_entries[m_ ## S1 ## _ ## S2]+=value;	

#define _write_dtr(S1, S2, value) \
	m_entries[m_ ## S1 ## _ ## S2]+=value;	

#define _write_JS(S1, S2, value) \
	m_entries[m_ ## S1 ## _ ## S2]+=value;	

#define _write_JD(S1, S2, value) \
	m_entries_old[m_ ## S1 ## _ ## S2]+=value;

#define _write_dac(S1, S2, value) \
	m_entries_old[m_ ## S1 ## _ ## S2]+=value;
/*--------------------------------------------------------------------------*/
#define _circuit_gdev                (OPT::gmin)
#define _circuit_gmin                (OPT::gmin)
#define _circuit_imax                1.0
#define _circuit_imelt               1.0
#define _circuit_iteration           1.0
#define _circuit_scale               1.0
#define _circuit_shrink              1.0
#define _circuit_simulatorSubversion 0
#define _circuit_simulatorVersion    3.5
#define _circuit_sourceScaleFactor   1.0
#define _circuit_tnom                (OPT::tnom_c + CONSTCtoK)
#define _ambient_temp                (OPT::temp_c + CONSTCtoK)
#define _circuit_temp                (CKT_BASE::_sim->_temp_c + CONSTCtoK)
#define _vt_nom                      (BOLTZMANN*_ambient_temp/ELECTRON_CHARGE)
#define _scale                       1.0
/*--------------------------------------------------------------------------*/
#define EXP90 1.220403294317841e+039
#define m00_hypot(v00,x,y)      v00 = sqrt((x)*(x)+(y)*(y));
#define m10_hypot(v10,v00,x,y)  v10 = (x)/(v00);
#define m11_hypot(v11,v00,x,y)  v11 = (y)/(v00);
#define m00_max(v00,x,y)        v00 = ((x)>(y))?(x):(y);
#define m10_max(v10,v00,x,y)    v10 = ((x)>(y))?1.0:0.0;
#define m11_max(v11,v00,x,y)    v11 = ((x)>(y))?0.0:1.0;
#define m00_min(v00,x,y)        v00 = ((x)<(y))?(x):(y);
#define m10_min(v10,v00,x,y)    v10 = ((x)<(y))?1.0:0.0;
#define m11_min(v11,v00,x,y)    v11 = ((x)<(y))?0.0:1.0;
#define m00_pow(v00,x,y)        v00 = pow(x,y);
#define m10_pow(v10,v00,x,y)    v10 = (x==0.0)?0.0:(v00)*(y)/(x);
#define m11_pow(v11,v00,x,y)    v11 = (x==0.0)?0.0:(log(x)*(v00));
#define m00_div(v00,v10,x,y)    double v10=1/(y); double v00=(x)*v10;
#define m10_div(v10,v00,vv,x,y)
#define m11_div(v11,v00,vv,x,y) double v11 = -v00*vv;
#define m00_mult(v00,v10,v11,x,y) double v10=(x); double v11=(y); double v00=v10*v11;
#define m00_add(v00,x,y)        double v00=(x)+(y);
#define m00_cos(v00,x)          v00 = cos(x);
#define m10_cos(v10,v00,x)      v10 = (-sin(x));
#define m00_sin(v00,x)          v00 = sin(x);
#define m10_sin(v10,v00,x)      v10 = (cos(x));
#define m00_tan(v00,x)          v00 = tan(x);
#define m10_tan(v10,v00,x)      v10 = (1.0/cos(x)/cos(x));
#define m00_cosh(v00,x)         v00 = cosh(x);
#define m10_cosh(v10,v00,x)     v10 = (sinh(x));
#define m00_sinh(v00,x)         v00 = sinh(x);
#define m10_sinh(v10,v00,x)     v10 = (cosh(x));
#define m00_tanh(v00,x)         v00 = tanh(x);
#define m10_tanh(v10,v00,x)     v10 = (1.0/cosh(x)/cosh(x));
#define m00_acos(v00,x)         v00 = acos(x);
#define m10_acos(v10,v00,x)     v10 = (-1.0/sqrt(1-x*x));
#define m00_asin(v00,x)         v00 = asin(x);
#define m10_asin(v10,v00,x)     v10 = (+1.0/sqrt(1-x*x));
#define m00_atan(v00,x)         v00 = atan(x);
#define m10_atan(v10,v00,x)     v10 = (+1.0/(1+x*x));
#define m00_logE(v00,x)         v00 = log(x);
#define m10_logE(v10,v00,x)     v10 = (1.0/x);
#define m00_log10(v00,x)        v00 = log10(x);
#define m10_log10(v10,v00,x)    v10 = (1.0/x/log(10));
#define m00_sqrt(v00,x)         v00 = sqrt(x);
#define m10_sqrt(v10,v00,x)     v10 = (0.5/v00);
#define m00_fabs(v00,x)         v00 = fabs(x);
#define m10_fabs(v10,v00,x)     v10 = (((x)>=0)?(+1.0):(-1.0));
#define m00_exp(v00,x)          v00 = exp(x); assert(is_number(x));
#define m10_exp(v10,v00,x)      v10 = v00;
#define m00_abs(v00)            ((v00)<(0)?(-(v00)):(v00))
#define m00_limexp(v00,x)       v00 = ((x)<90.0?exp(x):EXP90*(x-89.0));
#define m10_limexp(v10,v00,x)   v10 = ((x)<90.0?(v00):EXP90);
#define m20_logE(v00)         (-1.0/v00/v00)
#define m20_exp(v00)          exp(v00)
#define m20_limexp(v00)       ((v00)<90.0?exp(v00):0.0)
#define m20_sqrt(v00)         (-0.25/(v00)/sqrt(v00))
#define m20_fabs(v00)         0.0

/*--------------------------------------------------------------------------*/
#define _cos(val,arg)            val = cos(arg);
#define _d_cos(val,dval,arg)     val = cos(arg);     dval = (-sin(arg));
#define _sin(val,arg)            val = sin(arg);
#define _d_sin(val,dval,arg)     val = sin(arg);     dval = (cos(arg));
#define _tan(val,arg)            val = tan(arg);
#define _d_tan(val,dval,arg)     val = tan(arg);     dval = (1.0/cos(arg)/cos(arg));
#define _hypot(xy,x,y)           xy = sqrt((x)*(x)+(y)*(y));
#define _dx_hypot(dx,xy,x,y)     dx = (x)/(xy);
#define _dy_hypot(dy,xy,x,y)     dy = (y)/(xy);
#define _max(xy,x,y)             xy = ((x)>(y))?(x):(y);
#define _dx_max(dx,xy,x,y)       dx = ((x)>(y))?1.0:0.0;
#define _dy_max(dy,xy,x,y)       dy = ((x)>(y))?0.0:1.0;
#define _min(xy,x,y)             xy = ((x)<(y))?(x):(y);
#define _dx_min(dx,xy,x,y)       dx = ((x)<(y))?1.0:0.0;
#define _dy_min(dy,xy,x,y)       dy = ((x)<(y))?0.0:1.0;
#define _cosh(val,arg)           val = cosh(arg);
#define _d_cosh(val,dval,arg)    val = cosh(arg);    dval = (sinh(arg));
#define _sinh(val,arg)           val = sinh(arg);
#define _d_sinh(val,dval,arg)    val = sinh(arg);    dval = (cosh(arg));
#define _tanh(val,arg)           val = tanh(arg);
#define _d_tanh(val,dval,arg)    val = tanh(arg);    dval = (1.0/cosh(arg)/cosh(arg));
#define _acos(val,arg)           val = acos(arg);
#define _d_acos(val,dval,arg)    val = acos(arg);    dval = (-1.0/sqrt(1-arg*arg));
#define _asin(val,arg)           val = asin(arg);
#define _d_asin(val,dval,arg)    val = asin(arg);    dval = (+1.0/sqrt(1-arg*arg));
#define _atan(val,arg)           val = atan(arg);
#define _d_atan(val,dval,arg)    val = atan(arg);    dval = (+1.0/(1+arg*arg));
#define _logE(val,arg)           val = log(arg);
#define _d_logE(val,dval,arg)    val = log(arg);     dval = (1.0/arg);
#define _log10(val,arg)          val = log10(arg);
#define _d_log10(val,dval,arg)   val = log10(arg);   dval = (1.0/arg/log(10));
#define _exp(val,arg)            val = exp(arg);
#define _d_exp(val,dval,arg)     val = exp(arg);     dval = val;
#define _sqrt(val,arg)           val = sqrt(arg);
#define _d_sqrt(val,dval,arg)    val = sqrt(arg);    dval = (1.0/val/2.0);
#define _pow(xy,x,y)             xy = pow(x,y);
#define _dx_pow(dx,xy,x,y)       dx = (x==0.0)?0.0:((y/x)*xy);
#define _dy_pow(dy,xy,x,y)       dy = (x==0.0)?0.0:((log(x)/exp(0.0))*xy);

#define _div1(x,y)               ((x)/(y))
#define _div0(xy,x,y)            xy=(x)/(y);
#define _div(xy,dx,x,y)          dx=1/(y); xy=(x)*dx;
#define _dx_div(dx,xy,x,y)
#define _dy_div(dy,dx,xy,x,y)    dy = -xy*dx;

#define _limexp(val,arg)         val = ((arg)<(90)) ? (exp(arg)) : (exp(90)*(1.0+(arg-90)));
#define _d_limexp(val,dval,arg)  val = ((arg)<(90)) ? (exp(arg)) : (exp(90)*(1.0+(arg-90))); dval = val;
#define _fabs(val,arg)           val = fabs(arg);
#define _d_fabs(val,dval,arg)    val = fabs(arg);    dval = (((val)>=0)?(+1.0):(-1.0));

inline void _abs(double& val, const double& arg){val=abs(arg);}
inline double _abs(const double& arg){return abs(arg);}

/*
#define _(val,arg)         val = ((arg)<(90)) ? (exp(arg)) : (exp(90)*(1.0+(arg-90)));
#define _d_limexp(val,dval,arg)  val = ((arg)<90)) ? (exp(arg)) : (exp(90)*(1.0+(arg-90))); dval = val;
#define _fabs(val,arg)           val = fabs(arg);
#define _d_fabs(val,dval,arg)    val = fabs(arg);    dval = (((val)>=0)?(+1.0):(-1.0));
#define _abs(val)                ((val)<(0) ? (-(val)):(val))
*/



/*--------------------------------------------------------------------------*/
#define jacobian(a,b) m_required[m_##a##_##b]=true;	
#define static_jacobian4(p,q,r,s)  jacobian(p,r) jacobian(p,s) jacobian(q,r) jacobian(q,s)
#define static_jacobian2s(p,q,r)   jacobian(p,r) jacobian(q,r)
#define static_jacobian2p(p,r,s)   jacobian(p,r) jacobian(p,s)
#define static_jacobian1(p,r)      jacobian(p,r)
#define dynamic_jacobian4(p,q,r,s) jacobian(p,r) jacobian(p,s) jacobian(q,r) jacobian(q,s)
#define dynamic_jacobian2s(p,q,r)  jacobian(p,r) jacobian(q,r)
#define dynamic_jacobian2p(p,r,s)  jacobian(p,r) jacobian(p,s)
#define dynamic_jacobian1(p,r)     jacobian(p,r)
#define whitenoise_jacobian4(p,q,r,s)
#define whitenoise_jacobian2s(p,q,r)
#define whitenoise_jacobian2p(p,r,s)
#define whitenoise_jacobian1(p)
#define flickernoise_jacobian4(p,q,r,s)
#define flickernoise_jacobian2s(p,q,r)
#define flickernoise_jacobian2p(p,r,s)
#define flickernoise_jacobian1(p)
/*--------------------------------------------------------------------------*/


// temp hacks
#define _stop() assert(false)
// Boltzmann's constant in joules/kelvin
#define BOLTZMANN 1.3806503e-23
// coulombs of an electron
#define ELECTRON_CHARGE 1.602176462e-19

inline string toLower(string s){
	transform (s.begin(), s.end(), s.begin(), ::tolower);
	return s;
}
#define EXIT_IF_ISNAN(var) assert(is_number((double)var))

// -------------------------------------------------------------------------- //
#endif
