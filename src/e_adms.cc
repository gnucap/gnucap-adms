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
 */
#include <m_divdiff.h>
#include <u_xprobe.h>
#include <e_aux.h>
#include "e_adms.h"
/*--------------------------------------------------------------------------*/
ADMS_BASE::ADMS_BASE():
	COMPONENT(),
	_loaditer(0),
	_ev(0.),
	_dt(0.),
	_method_u(meUNKNOWN), _method_a(mTRAPGEAR) 
{
	_n = _nodes;

	std::fill_n(_time, int(OPT::_keep_time_steps), 0.);
}
/*--------------------------------------------------------------------------*/
ADMS_BASE::ADMS_BASE(const ADMS_BASE& p)
: COMPONENT(p),
	_loaditer(0),
	_ev(0.),
	_dt(0.),
	_method_u(p._method_u), _method_a(p._method_a)
{
	trace0(long_label().c_str());
	_n = _nodes;
	if (p._n == p._nodes) {
		for (int ii = 0;  ii < NODES_PER_BRANCH;  ++ii) {
			_n[ii] = p._n[ii];
		}
	}else{
		assert(p._nodes);
		// the constructor for a derived class will take care of it
	}

	notstd::copy_n(p._time, int(OPT::_keep_time_steps), _time);
}
/*--------------------------------------------------------------------------*/
bool ADMS_BASE::skip_dev_type(CS& cmd)
{
	return cmd.umatch(dev_type() + ' ');
}
/*--------------------------------------------------------------------------*/
void ADMS_BASE::tr_begin()
{
	trace1("ADMS_BASE::tr_begin for", short_label());
	for (int i=OPT::_keep_time_steps-1; i>=0; --i) {
		_time[i] = 0;
	}
	_dt = NOT_VALID;

	//from st.
	_method_a = method_select[OPT::method][_method_u];
}
/*--------------------------------------------------------------------------*/
void ADMS_BASE::tr_restore()
{
	untested();
	trace0(("ADMS_BASE::tr_restore for " + short_label()).c_str());
	if (_time[0] > _sim->_time0) {itested();
		trace0("shift back");
		for (int i=0  ; i<OPT::_keep_time_steps-1; ++i) {itested();
			_time[i] = _time[i+1];
			// _y[i] = _y[i+1]; FIXME!
		}
		_time[OPT::_keep_time_steps-1] = 0.;
		// _y[OPT::_keep_time_steps-1]    = FPOLY1(0., 0., 0.);
	}else if (_time[0] == _sim->_time0) {
		trace2( "*no shift ", _time[0] , _sim->_time0 );
	}else{untested();
	}

	//if( _time[0] != _sim->_time0 )
	//  assert(_time[0] == _sim->_time0);
	if (_time[0] != _sim->_time0) {itested();
		error(bDANGER, "//BUG// restore time mismatch.  t0=%.12f, s->t=%.12f\n",
				_time[0], _sim->_time0);
		//BUG// happens when continuing after a ^c,
		// when the last step was not printed
		// _time[0] is the non-printed time.  _sim->_time0 is the printed time.
	}else{
	}

	for (int i=OPT::_keep_time_steps-1; i>0; --i) {
		assert(_time[i] < _time[i-1] || _time[i] == 0.);
	}

	_method_a = method_select[OPT::method][_method_u];
}
/*--------------------------------------------------------------------------*/
void ADMS_BASE::dc_advance() // from elt.
{
	trace2( "ADMS_BASE::dc_advance " , long_label(), net_nodes());
	assert(_sim->_time0 == 0.); // DC

	for (int i=OPT::_keep_time_steps-1; i>=0; --i) {
		_time[i] = _sim->_time0;
	}

	_dt = NOT_VALID;

	//from storag
	// for (int i = 1;  i < OPT::_keep_time_steps;  ++i) {
	//  _i[i] = _i[0];
	// }
	//
}
/*--------------------------------------------------------------------------*/
void ADMS_BASE::tr_regress()
{
  	trace2("ADMS_BASE::tr_regress", long_label(), net_nodes());
	assert(_time[0] >= _sim->_time0); // moving backwards
	assert(_time[1] <= _sim->_time0); // but not too far backwards

	for (int i=OPT::_keep_time_steps-1; i>0; --i) {
		assert(_time[i] < _time[i-1] || _time[i] == 0.);
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

	DISCONT d = disNONE;
	for (unsigned i=0;i < net_nodes(); ++i) {
		d |= _n[i]->discont();
	}

#ifdef USE_DTIME
	if (_time[1] == _time[0])
#else
	if (_time[1] == 0)
#endif
	{
		m = meEULER;
	} else if (d) {
		m = meEULER;
	} else {
	}
	_method_a = method_select[OPT::method][m];
#endif
}
/*--------------------------------------------------------------------------*/
double ADMS_BASE::tr_review_trunc_error(const FPOLY1* q)
{
	int error_deriv = order()+1;
	double timestep;
	trace3("ADMS_BASE::tr_review_trunc_error", order(), error_deriv, _time[0]);
	if (_time[0] <= 0.) {
		// DC, I know nothing
		timestep = NEVER;
#ifdef USE_DTIME
	}else if (_time[error_deriv-1] == _time[error_deriv]) {
#else
	}else if (_time[error_deriv] <= 0.) {
#endif
		// first few steps, I still know nothing
		// repeat whatever step was used the first time
		timestep = _dt;
	}else{
		for (int i=error_deriv; i>0; --i) {
			assert(_time[i] < _time[i-1]); // || _time[i] == 0.);
		}

		double c[OPT::_keep_time_steps];
		for (int i=0; i<OPT::_keep_time_steps; ++i) {
			c[i] = q[i].f0;
		}
		assert(error_deriv < OPT::_keep_time_steps);
		derivatives(c, OPT::_keep_time_steps, _time);
		// now c[i] is i'th derivative

		assert(OPT::_keep_time_steps >= 5);
		trace0(("ts" + long_label()).c_str());
		trace5("time", _time[0], _time[1], _time[2], _time[3], _time[4]);
		trace5("charge", q[0].f0, q[1].f0, q[2].f0, q[3].f0, q[4].f0);
		trace5("deriv", c[0], c[1], c[2], c[3], c[4]);

		if (c[error_deriv] == 0) {
			timestep = NEVER;
		}else{
			double chargetol = std::max(OPT::chgtol,
					OPT::reltol * std::max((double)std::abs(q[0].f0), (double) std::abs(q[1].f0)));
			double tol = OPT::trtol * chargetol;
			double denom = error_factor() * std::abs(c[error_deriv]);
			assert(tol > 0.);
			assert(denom > 0.);
			switch (error_deriv) { // pow is slow.
				case 1:	timestep = tol / denom;		break;
				case 2:	timestep = sqrt(tol / denom);	break;
				case 3:	timestep = cbrt(tol / denom);	break;
				default:	timestep = pow((tol / denom), 1./(error_deriv)); break;
			}
			trace4("", chargetol, tol, denom, timestep);
		}
	}
	assert(timestep > 0.);
	return timestep;
}
/*--------------------------------------------------------------------------*/
double ADMS_BASE::tr_review_check_and_convert(double timestep)
{
	double time_future;
	if (timestep == NEVER) {
		time_future = NEVER;
	}else{
		if (timestep < _sim->_dtmin) {
			timestep = _sim->_dtmin;
		}else{
		}

		if (timestep < _dt * OPT::trreject) {
			if (_time[order()] == 0) {
				error(bWARNING, "initial step rejected: " + long_label() + '\n');
				error(bWARNING, "new=%g  old=%g  required=%g\n",
						timestep, _dt, _dt * OPT::trreject);
			}else{
				error(bTRACE, "step rejected: " + long_label() + '\n');
				error(bTRACE, "new=%g  old=%g  required=%g\n",
						timestep, _dt, _dt * OPT::trreject);
			}
			time_future = _time[1] + timestep;
			trace3("reject", timestep, _dt, time_future);
		}else{
			time_future = _time[0] + timestep;
			trace3("accept", timestep, _dt, time_future);
		}
	}
	assert(time_future > 0.);
	assert(time_future > _time[1]);
	return time_future;
}
/*--------------------------------------------------------------------------*/
 // * ripped from e_storag.h
/*--------------------------------------------------------------------------*/
/* table for selecting local integraton method
 * Determines which one wins in a conflict.
 * "only" wins over non-only.  local (_method_u) wins over opt.
 */
//                     OPT::method    _method_u
METHOD ADMS_BASE::method_select[meNUM_METHODS][meNUM_METHODS] = {
	/*vv OPT vv*/
	//local>>>EULER,EULERONLY,TRAP,TRAPONLY,GEAR2,GEAR2ONLY,TRAPGEAR,TRAPEULER
	/*meUNKNOWN*/
	{mTRAPGEAR,mEULER,mEULER,mTRAP, mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER},
	/*meEULER*/
	{mEULER,   mEULER,mEULER,mTRAP, mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER},
	/*meEULERONLY*/
	{mEULER,   mEULER,mEULER,mEULER,mTRAP,mEULER,mGEAR,mEULER,   mEULER},
	/*meTRAP*/
	{mTRAP,    mEULER,mEULER,mTRAP, mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER},
	/*meTRAPONLY*/
	{mTRAP,    mTRAP, mEULER,mTRAP, mTRAP,mTRAP, mGEAR,mTRAP,    mTRAP},
	/*meGEAR*/
	{mGEAR,    mEULER,mEULER,mTRAP, mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER},
	/*meGEAR2ONLY*/
	{mGEAR,    mGEAR, mEULER,mGEAR, mTRAP,mGEAR, mGEAR,mGEAR,    mGEAR},
	/*meTRAPGEAR*/
	{mTRAPGEAR,mEULER,mEULER,mTRAP, mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER},
	/*meTRAPEULER*/
	{mTRAPEULER,mEULER,mEULER,mTRAP,mTRAP,mGEAR, mGEAR,mTRAPGEAR,mTRAPEULER}
};
/*--------------------------------------------------------------------------*/
/* tr_needs_eval: check to see if this device needs to be evaluated
 * this works, and saves significant time
 * but possible errors.
 * Errors do not seem significant, but I can't tell without more data.
 * To disable:  option nolcbypass
 */
bool ADMS_BASE::tr_needs_eval()const
{
	incomplete();
	return 0;
	//assert(!is_q_for_eval());
#if 0
	return (!OPT::lcbypass
			|| !converged() 
			|| _sim->is_advance_or_first_iteration()
			|| !conchk(_y[0].x, tr_input(), OPT::abstol)
			|| _sim->uic_now());
#endif
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

