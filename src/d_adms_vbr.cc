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
 * DEV_CPOLY_V, derived from DEV_CPOLY_G
 * number of nodes = boff + 2*(n_ports=n_vports+n_iports)
 * number of val, ov = n_ports+1
 * val[0] is the constant part, val[1] is self admittance,
 *   val[2+] are transadmittances, up to n_ports
 * node[0] and node[1] are the output (voltage)
 * node[2] is the branch node (if boff). this is a hack.
 * node[boff+3] up are inputs.
 * node[boff+2*i+1] and node[2*i+2] correspond to val[i+1]
 * nodes beyond boff+_n_ports*2 correspond to controlling currents
 */
#include <globals.h>
#include "e_admsrc.h"
#include "gcuf_compat.h"

// the voltage source has a branch node.
#define boff 1

#ifdef boff
// must be IN1, it is used as inode
#define BR IN1
#else
incomplete
#define boff 0
#endif
/*--------------------------------------------------------------------------*/
namespace {
/*--------------------------------------------------------------------------*/
class DEV_CPOLY_V : public ADMS_SOURCE {
protected:
  double   _time;
protected:
  explicit DEV_CPOLY_V(const DEV_CPOLY_V& p);
public:
  explicit DEV_CPOLY_V();
  ~DEV_CPOLY_V();
protected: // override virtual
  char	   id_letter()const	{unreachable(); return '\0';}
  std::string value_name()const	{incomplete(); return "";}
  std::string dev_type()const	{unreachable(); return "cpoly_g";}

  bool	   has_iv_probe()const  {return 1-boff;}
  bool	   has_inode()const  {return boff;}

  CARD*	   clone()const		{return new DEV_CPOLY_V(*this);}
  void	   tr_iwant_matrix();
  void	   tr_begin();
  bool	   do_tr();
  void	   tr_load();
  void	   tr_unload();
  double   tr_involts()const	{unreachable(); return NOT_VALID;}
  double   tr_involts_limited()const {unreachable(); return NOT_VALID;}
  double   tr_amps()const;
  void	   ac_iwant_matrix();
  void	   ac_load();
  COMPLEX  ac_involts()const	{itested(); return NOT_VALID;}
  COMPLEX  ac_amps()const	{itested(); return NOT_VALID;}

  std::string port_name(uint_t)const {untested();
    incomplete();
    unreachable();
    return "";
  }
public:
  void expand();
protected:
  bool do_tr_con_chk_and_q();
private:
  double _one0, _one1;
};
/*--------------------------------------------------------------------------*/
DEV_CPOLY_V::DEV_CPOLY_V(const DEV_CPOLY_V& p)
  :ADMS_SOURCE(p),
   _time(NOT_VALID)
{
  // not really a copy .. only valid to copy a default
  // too lazy to do it right, and that's all that is being used
  // to do it correctly requires a deep copy
  // just filling in defaults is better than a shallow copy, hence this:
  assert(!p._values);
  assert(!p._old_values);
  assert(p._n_ports == 0);
  assert(p._n_iports == 0);
  assert(!p._inputs);
  _boff = boff;
}
/*--------------------------------------------------------------------------*/
DEV_CPOLY_V::DEV_CPOLY_V()
  :ADMS_SOURCE(),
   _time(NOT_VALID)
{
  _boff = boff;
}
/*--------------------------------------------------------------------------*/
DEV_CPOLY_V::~DEV_CPOLY_V()
{
}
/*--------------------------------------------------------------------------*/
bool DEV_CPOLY_V::do_tr_con_chk_and_q()
{
  q_load();

  assert(_old_values);
  set_converged(conchk(_time, _sim->_time0));
  _time = _sim->_time0;
  for (uint_t i=0; converged() && i<=_n_ports; ++i) {
    set_converged(conchk(_old_values[i], _values[i]));
  }
  for (uint_t i=0; converged() && i<_n_iports; ++i) {
    // incomplete
  }
  return converged();
}
/*--------------------------------------------------------------------------*/
#if 0
void ELEMENT::tr_iwant_matrix_extended()
{
  trace4("ELEMENT::tr_iwant_matrix_extended", long_label(), dev_type(),  ext_nodes(), int_nodes());
  assert(is_device());
  assert(ext_nodes() + int_nodes() == matrix_nodes());

  for (uint_t ii = 0;  ii < matrix_nodes();  ++ii) {
      trace2("ELEMENT::tr_iwant_matrix_extended", ii, _n[ii].m_() );
  }
  for (uint_t ii = 0;  ii < matrix_nodes();  ++ii) {
    if (_n[ii].m_()  != INVALID_NODE) {
      for (uint_t jj = 0;  jj < ii ;  ++jj) {
	_sim->_aa.iwant(_n[ii].m_(),_n[jj].m_());
	_sim->_lu.iwant(_n[ii].m_(),_n[jj].m_());
      }
    }else{itested();
      // node 1 is grounded or invalid
    }
  }
}
#endif

void DEV_CPOLY_V::tr_iwant_matrix()
{
  // FIXME: this is likely too much.
  trace2("tr_iwant_matrix", long_label(), matrix_nodes());
  tr_iwant_matrix_extended(); // uses matrix_nodes

  _sim->_aa.iwant(_n[BR].m_(),_n[OUT1].m_());
  _sim->_lu.iwant(_n[BR].m_(),_n[OUT1].m_());
  _sim->_aa.iwant(_n[BR].m_(),_n[OUT2].m_());
  _sim->_lu.iwant(_n[BR].m_(),_n[OUT2].m_());
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::tr_begin()
{
  ELEMENT::tr_begin();
  _one0 = 1.;
  _one1 = 0.;
}
/*--------------------------------------------------------------------------*/
bool DEV_CPOLY_V::do_tr()
{
  assert(_values);

  _m0 = CPOLY1(0., _values[0], _values[1]); // hmm really?!

  do_tr_con_chk_and_q();
  if(_n_iports){
    _sim->_late_evalq.push_back(this);
    trace1("ADMS_V::do_tr push late eval", converged());
    return true;
  }else{
    trace1("ADMS_V::do_tr", converged());
    return converged();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::tr_load()
{
  node_t gnd(&ground_node);
  trace4("DEV_CPOLY_V::tr_load", _n_vports, _n_iports, _values[0], _values[1]);
  assert(_n_vports+1 == _n_ports);

  double d = dampdiff(&_m0.c1, _m1.c1);
  if (d != 0.) { untested();
    _sim->_aa.load_symmetric(_n[BR].m_(), 0, -d);
  }else{
  }
  //
  // // load source...
  //
  _old_values[0] = _values[0];
  _old_values[1] = _values[1];
  // inject current
  if(!boff) incomplete();

  // voltage control, like poly_g, but with branch node
  for (uint_t i=2; i<=_n_ports; ++i) {
    trace3("DEV_CPOLY_V::tr_load vc", long_label(), i, _values[i]);
    tr_load_extended(gnd, _n[BR], _n[2*i-2+boff], _n[2*i-1+boff], &(_values[i]), &(_old_values[i]));
  }

  double* ival = _values+2+_n_vports;
  double* ioval = _old_values+2+_n_vports;
  node_t* ni = _n+2*_n_ports + boff; // start of current input ports.

  trace4("loading control", long_label(), _n_iports, _n_vports, _n_ports);
  for (uint_t i=0; i<_n_iports; ++i) { untested();
    // something like tr_load_active
    trace4("loading control", i, ival[i], ioval[i], _inputs[i]->has_iv_probe());
    if (_inputs[i]->has_iv_probe()) { untested();
      ELEMENT const* ie = prechecked_cast<ELEMENT* const>(_inputs[i]);
      assert(ie);

      double d = dampdiff(&(ival[i]), ioval[i]);
      if(d!=0.){
	_sim->_aa.load_asymmetric(_n[BR].m_(), 0,  ni[2*i].m_(), ni[2*i+1].m_(), d);
      }
      ioval[i]=ival[i];
    }else{ incomplete();
      assert(_inputs[i]->has_inode());
      double d = dampdiff(&(ival[i]), ioval[i]);
      if(d!=0.){ untested();
	_sim->_aa.load_asymmetric(_n[BR].m_(), 0,  ni[2*i].m_(), ni[2*i+1].m_(), d);
      }
      ioval[i]=ival[i];
    }
  }

  d = _one0 - _one1;
  if (!_sim->is_advance_or_first_iteration()) {
    _one0 = _one1 + d;
  }
  d = ((_sim->is_inc_mode()) ? d : _one0);

  // convert current to voltage between OUT1 OUT2
  if (d != 0.) {
    _sim->_aa.load_point(_n[BR].m_(), _n[OUT1].m_(), d);
    _sim->_aa.load_point(_n[BR].m_(), _n[OUT2].m_(), -d);
    _sim->_aa.load_point(_n[OUT1].m_(), _n[BR].m_(), d);
    _sim->_aa.load_point(_n[OUT2].m_(), _n[BR].m_(), -d);
  }

  _one1 = _one0;

/*--load_source------------------------------------------------------------------------*/
  trace3("DEV_BVS::tr_load", long_label(), _m0.c0, _m1.c0);
  assert(_m0.c0 == _m0.c0);
  assert (_n[BR].m_() != 0);
  d = dampdiff(&_m0.c0, _m1.c0); // hmmm
  trace1("DEV_BVS::tr_load source",d);
  if (d != 0.) {
    _n[BR].i() += d;
  }else{
  }
  _m1 = _m0;
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::tr_unload()
{
  ADMS_SOURCE::tr_unload();
  _one0 = 0;
}
/*--------------------------------------------------------------------------*/
double DEV_CPOLY_V::tr_amps()const
{
  if(!boff) incomplete();
  double amps = 0; // _m0.c0; // == _values[0]?
  return dn_diff(_n[BR].v0(), 0.); // _values[1];
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::ac_iwant_matrix()
{ incomplete();
  ac_iwant_matrix_extended();

  _sim->_acx.iwant(_n[BR].m_(),_n[OUT1].m_());
  _sim->_acx.iwant(_n[BR].m_(),_n[OUT2].m_());
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::ac_load()
{ incomplete();
  _acg = _values[1];
  ac_load_passive();
  node_t gnd(&ground_node);
  if(!boff) incomplete();
  for (uint_t i=2; i<=_n_ports; ++i) { untested();
    ac_load_extended(_n[BR], gnd, _n[2*i-2], _n[2*i-1], _values[i]);
  }
  for (uint_t i=0; i<_n_iports; ++i) {
    incomplete();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::expand()
{
  if (!subckt()) {
    new_subckt(); // hmm probably not a good idea
  }else{ untested();
  }

  assert(matrix_nodes() == ext_nodes()+boff);
  assert(int_nodes() == boff);

#if boff
  assert(BR);
  if (!(_n[BR].n_())) { untested();
    _n[BR].new_model_node( "br", this);
  }else{ untested();
  }
#endif
  COMPONENT::expand();
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_CPOLY_V p4;
DISPATCHER<CARD>::INSTALL d4(&device_dispatcher, "adms_v", &p4);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
