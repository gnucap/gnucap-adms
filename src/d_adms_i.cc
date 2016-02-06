/*$Id: d_poly_g.cc,v 26.137 2010/04/10 02:37:05 al Exp $ -*- C++ -*-
 * Copyright (C) 2001 Albert Davis
 * Author: Albert Davis <aldavis@gnu.org>
 *
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
 * ADMS_I
 * number of nodes = 2*n_ports
 * number of val, ov = n_ports+1
 * val[0] is the constant part, val[1] is self admittance,
 *   val[2+] are transadmittances, up to n_ports
 * node[0] and node[1] are the output.
 * node[2] up are inputs.
 * node[2*i] and node[2*i+1] correspond to val[i+1]
 *
 */
#include <globals.h>
#include <io_trace.h>
#include "e_admsrc.h"
#include "gcuf_compat.h"
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class ADMS_I : public ADMS_SOURCE {
protected:
  double   _time;
protected:
  explicit ADMS_I(const ADMS_I& p);
public:
  explicit ADMS_I();
//  ~ADMS_I();
protected: // override virtual
  char	   id_letter()const	{unreachable(); return '\0';}
  std::string value_name()const	{incomplete(); return "";}
  std::string dev_type()const	{unreachable(); return "adms_i";}
  CARD*	   clone()const		{return new ADMS_I(*this);}
  bool	   has_iv_probe()const  {return true;}

  void	   tr_iwant_matrix();
  bool	   do_tr();
  void	   tr_load();
// void	   tr_unload();
  double   tr_involts()const	{unreachable(); return NOT_VALID;}
  double   tr_involts_limited()const {unreachable(); return NOT_VALID;}
  double   tr_amps()const;
  void	   ac_iwant_matrix()	{incomplete(); ac_iwant_matrix_extended();}
  void	   ac_load();
  COMPLEX  ac_involts()const	{itested(); return NOT_VALID;}
  COMPLEX  ac_amps()const	{itested(); return NOT_VALID;}

  std::string port_name(uint_t)const {untested();
    incomplete();
    unreachable();
    return "";
  }
protected:
  bool do_tr_con_chk_and_q();
};
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
ADMS_I::ADMS_I(const ADMS_I& p)
  :ADMS_SOURCE(p),
   _time(NOT_VALID)
{
  // not really a copy .. only valid to copy a default
  // too lazy to do it right, and that's all that is being used
  // to do it correctly requires a deep copy
  // just filling in defaults is better than a shallow copy, hence this:
  assert(!p._values);
  assert(!p._old_values);
  assert(p._n_vports == 0);
  assert(p._n_iports == 0);
  assert(!p._inputs);
}
/*--------------------------------------------------------------------------*/
ADMS_I::ADMS_I()
  :ADMS_SOURCE(),
   _time(NOT_VALID)
{
  trace2("ADMS_I::ADMS_I", _m0, _m1);
}
/*--------------------------------------------------------------------------*/
bool ADMS_I::do_tr_con_chk_and_q()
{
  q_load();

  assert(_old_values);
  set_converged(conchk(_time, _sim->_time0));
  _time = _sim->_time0;
  for (uint_t i=0; converged() && i<2+_n_vports; ++i) {
    // FIXME: abs/reltol
    set_converged(conchk(_old_values[i], _values[i]));
    trace5("ADMS_I::conchk_and_q", i, converged(), _values[i], _old_values[i], _n_vports);
  }
  return converged();
}
/*--------------------------------------------------------------------------*/
void ADMS_I::tr_iwant_matrix()
{
  // requests too much!
  tr_iwant_matrix_extended();
  tr_iwant_matrix_active(); // in ADMS_SOURCE

}
/*--------------------------------------------------------------------------*/
bool ADMS_I::do_tr()
{
  trace5("ADMS_I::do_tr", long_label(),
      _values[0], _old_values[0],
      _values[1], _old_values[1]);
  trace2("ADMS_I::do_tr", _m0, _m1);
  assert(_values);
  _m0 = CPOLY1(0., _values[0], _values[1]);
  bool P = converged();

  do_tr_con_chk_and_q(); // _values[2+]
  if(_n_iports){
    _sim->_late_evalq.push_back(this);
    trace1("ADMS_I::do_tr push", converged());
  }else{
    trace1("ADMS_I::do_tr", converged());
    P = converged();
  }
  trace5("ADMS_I::do_tr done", long_label(),
                               _values[0], _old_values[0],
                               _values[1], _old_values[1]);
  return P;
}
/*--------------------------------------------------------------------------*/
void ADMS_I::tr_load()
{
  trace4("ADMS_I::tr_load", _values[0], _old_values[0],
                            _values[1], _old_values[1]);
  trace4("ADMS_I::tr_load", _m0, _m1, _n_vports, _n_iports);
  trace1("ADMS_I::tr_load rhs", _m0.c0);
  tr_load_passive();
  _old_values[0] = _values[0];
  _old_values[1] = _values[1];

  for (uint_t i=2; i<=1+_n_vports; ++i) {
  // like d_poly_g
    tr_load_extended(_n[OUT1], _n[OUT2], _n[2*i-2], _n[2*i-1], &(_values[i]), &(_old_values[i]));
  }

  node_t* ni = _n+2+2*_n_vports;
  double* ival = _values+2+_n_vports;
  double* ioval = _old_values+2+_n_vports;

  for (uint_t i=0; i<matrix_nodes(); ++i) {
    trace2("MATRIX",i, _n[i].m_());
  }

  for (uint_t i=0; i<_n_iports; ++i) { itested();
    // something like tr_load_active
    if (_inputs[i]->has_iv_probe()) { itested();
      ELEMENT const* ie = prechecked_cast<ELEMENT* const>(_inputs[i]);
      assert(ie);
      double d = dampdiff(&(ival[i]), ioval[i]); //  * ie->_m0.c1;
      if (d != 0.) {
	_sim->_aa.load_asymmetric(_n[OUT1].m_(), _n[OUT2].m_(),  ni[2*i].m_(), ni[2*i+1].m_(), d);
      }
      ioval[i]=ival[i];
    }else if(_inputs[i]->has_inode()){ untested();
      incomplete();
    }
  }

/// cccs
//  double d = dampdiff(&_m0.c1, _m1.c1);
//  if (d != 0.) {
//    _sim->_aa.load_asymmetric(_n[OUT1].m_(), _n[OUT2].m_(),
//		       _n[IN1].m_(), _n[IN2].m_(), d);
}
/*--------------------------------------------------------------------------*/
double ADMS_I::tr_amps()const
{
  double amps = _m0.c0;
  double mul = 1.;
  for (uint_t i=1; i<2+_n_vports; ++i) {
    // don't use _values. they have just been zeroed!
    amps += dn_diff(_n[2*i-2].v0(),_n[2*i-1].v0()) * _old_values[i];
  }
  for (uint_t i=0; i<_n_iports; ++i) {
    if (_inputs[i]==this){ incomplete();
      assert(!i);
      mul = 1./(1.-_old_values[2+_n_vports]);
    }else{
      ADMS_SOURCE const* S = prechecked_cast<ADMS_SOURCE const*>(_inputs[i]);
      amps += S->tr_amps() * _old_values[i+2+_n_vports];
    }
  }
  return amps*mul;
}
/*--------------------------------------------------------------------------*/
void ADMS_I::ac_load()
{
  _acg = COMPLEX(_values[1], _old_values[1] * _sim->_jomega.imag());
  trace2("ADMS_I::ac_load", long_label(), _acg);
  ac_load_passive();
  for (uint_t i=2; i<2+_n_vports; ++i) {
    ac_load_extended(_n[OUT1], _n[OUT2], _n[2*i-2], _n[2*i-1],
	COMPLEX(_values[i], _old_values[i] * _sim->_jomega.imag()));
  }
  if(_n_iports){ incomplete();
  }
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
ADMS_I p4;
DISPATCHER<CARD>::INSTALL d4(&device_dispatcher, "adms_i", &p4);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
