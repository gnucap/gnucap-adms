/*
 * Copyright (C) 2015-2017 Felix Salfelder
 *               2001 Albert Davis
 * Author: Felix Salfelder <felix@salfelder.org>
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
 * node[2] up are inputs.
 * node[2*i+1] and node[2*i+2] correspond to val[i+1]
 * nodes beyond _n_ports*2 correspond to controlling currents
 */
#include <globals.h>
#include "e_admsrc.h"
#include "gcuf_compat.h"

// the voltage source has a branch node.
#define boff 1

#ifdef boff
// must be IN1, it is used as inode
// #define BR IN1
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
  std::string dev_type()const	{unreachable(); return "cpoly_v";}

  bool	   has_iv_probe()const  {return 1-boff;}
  bool	   has_inode()const  {return boff;}
  virtual node_t get_inode() const{ untested();
    return _n[BR()];
  }

  CARD*	   clone()const		{return new DEV_CPOLY_V(*this);}
  void	   tr_iwant_matrix();
  void	   tr_iwant_matrix_extended_branch();
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
public: //hmmm
  void map_nodes()
  { untested();
    COMPONENT::map_nodes();
    //_nBR.map();
  }
protected:
  bool do_tr_con_chk_and_q();
private:
  double _one0, _one1;
private:
  //node_t _nBR;
  unsigned BR() const{
    assert(ext_nodes());
    return ext_nodes();
  }

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
  // notstd::copy_n(&p._nBR, 1, &_nBR);
  _boff = p._boff;
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
void DEV_CPOLY_V::tr_iwant_matrix_extended_branch()
{ untested();
  trace4("ELEMENT::tr_iwant_matrix_extended", long_label(), dev_type(),  ext_nodes(), int_nodes());
  assert(is_device());
  assert(ext_nodes() + int_nodes() == matrix_nodes());

  for (uint_t ii = 0;  ii < matrix_nodes();  ++ii) { untested();
      trace2("ELEMENT::tr_iwant_matrix_extended", ii, _n[ii].m_() );
  }
  for (uint_t ii = 2;  ii < matrix_nodes();  ++ii) { untested();
    // connect all to branch..
    _sim->_aa.iwant(_n[BR()].m_(),_n[ii].m_());
    _sim->_lu.iwant(_n[BR()].m_(),_n[ii].m_());

    // is this too much?
    if (_n[ii].m_()  != INVALID_NODE) { untested();
      for (uint_t jj = 2;  jj < ii ;  ++jj) { untested();
	_sim->_aa.iwant(_n[ii].m_(),_n[jj].m_());
	_sim->_lu.iwant(_n[ii].m_(),_n[jj].m_());
      }
    }else{ untested();
      trace3("eek", ii, _n[ii].m_(), long_label() );
      // node 1 is grounded or invalid
    }
  }
}
/*--------------------------------------------------------------------------*/

void DEV_CPOLY_V::tr_iwant_matrix()
{
  trace3("tr_iwant_matrix", long_label(), matrix_nodes(), BR());
  assert(!subckt());
  if(boff){
    tr_iwant_matrix_extended_branch();
  }else{ untested();
    // FIXME: this is likely too much.
    tr_iwant_matrix_extended(); // uses matrix_nodes
  }

  _sim->_aa.iwant(_n[BR()].m_(),_n[OUT1].m_());
  _sim->_lu.iwant(_n[BR()].m_(),_n[OUT1].m_());
  _sim->_aa.iwant(_n[BR()].m_(),_n[OUT2].m_());
  _sim->_lu.iwant(_n[BR()].m_(),_n[OUT2].m_());
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
  trace2("DEV_CPOLY_V::tr_load", matrix_nodes(), BR());
  trace4("DEV_CPOLY_V::tr_load", _n_vports, _n_iports, _values[0], _values[1]);
  assert(_n_vports+1 == _n_ports);

  double d = dampdiff(&_m0.c1, _m1.c1);
  if (d != 0.) {
    _sim->_aa.load_symmetric(_n[BR()].m_(), 0, -d);
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
    tr_load_extended(gnd, _n[BR()], _n[2*i-2], _n[2*i-1], &(_values[i]), &(_old_values[i]));
  }

  double* ival = _values+2+_n_vports;
  double* ioval = _old_values+2+_n_vports;
  node_t* ni = _n+2*_n_ports; // start of current input ports.

  trace4("loading control", long_label(), _n_iports, _n_vports, _n_ports);
  for (uint_t i=0; i<_n_iports; ++i) {
    // something like tr_load_active
    trace4("loading control", i, ival[i], ioval[i], _inputs[i]->has_iv_probe());
    if (_inputs[i]->has_iv_probe()) { untested();
      ELEMENT const* ie = prechecked_cast<ELEMENT* const>(_inputs[i]);
      assert(ie);

      double d = dampdiff(&(ival[i]), ioval[i]);
      if(d!=0.){ untested();
	_sim->_aa.load_asymmetric(_n[BR()].m_(), 0,  ni[2*i].m_(), ni[2*i+1].m_(), d);
      }else{untested();
      }
      ioval[i] = ival[i];
    }else{ incomplete();
      assert(_inputs[i]->has_inode());
      double d = dampdiff(&(ival[i]), ioval[i]);
      if(d!=0.){ untested();
	_sim->_aa.load_asymmetric(_n[BR()].m_(), 0,  ni[2*i].m_(), ni[2*i+1].m_(), d);
      }else{
	untested();
      }
      ioval[i] = ival[i];
    }
  }

  d = _one0 - _one1;
  if (!_sim->is_advance_or_first_iteration()) {
    _one0 = _one1 + d;
  }else{
  }
  d = ((_sim->is_inc_mode()) ? d : _one0);

  // convert current to voltage between OUT1 OUT2
  trace2("hmm", BR(), _n[BR()].m_());
  if (d != 0.) {
    assert(_n[BR()].m_() != INVALID_NODE);
    assert(_n[OUT1].m_() != INVALID_NODE);
    assert(_n[OUT2].m_() != INVALID_NODE);
    trace5("load", long_label(), _n[BR()].m_(), _n[OUT1].m_(), _n[OUT2].m_(), d);
    _sim->_aa.load_point(_n[BR()].m_(), _n[OUT1].m_(), d);
    _sim->_aa.load_point(_n[BR()].m_(), _n[OUT2].m_(), -d);
    _sim->_aa.load_point(_n[OUT1].m_(), _n[BR()].m_(), d);
    _sim->_aa.load_point(_n[OUT2].m_(), _n[BR()].m_(), -d);
  }else{
  }

  _one1 = _one0;

/*--load_source------------------------------------------------------------------------*/
  trace5("DEV_BVS::tr_load", long_label(), _m0.c0, _m1.c0, ext_nodes(), BR());
  assert(_m0.c0 == _m0.c0);
  assert(ext_nodes());
  assert (_n[BR()].m_() != 0);
  d = dampdiff(&_m0.c0, _m1.c0); // hmmm
  trace1("DEV_BVS::tr_load source",d);
  if (d != 0.) {
    _n[BR()].i() += d;
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
  return dn_diff(_n[BR()].v0(), 0.); // _values[1];
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::ac_iwant_matrix()
{ incomplete();
  ac_iwant_matrix_extended();

  _sim->_acx.iwant(_n[BR()].m_(),_n[OUT1].m_());
  _sim->_acx.iwant(_n[BR()].m_(),_n[OUT2].m_());
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::ac_load()
{ incomplete();
  _acg = _values[1];
  ac_load_passive();
  node_t gnd(&ground_node);
  if(!boff) incomplete();
  for (uint_t i=2; i<=_n_ports; ++i) { untested();
    ac_load_extended(_n[BR()], gnd, _n[2*i-2], _n[2*i-1], _values[i]);
  }
  for (uint_t i=0; i<_n_iports; ++i) {
    incomplete();
  }
}
/*--------------------------------------------------------------------------*/
void DEV_CPOLY_V::expand()
{ untested();
  if (!subckt()) {
    // new_subckt(); // hmm probably not a good idea
  }else{ untested();
  }

  assert(matrix_nodes() == ext_nodes()+boff);
  assert(int_nodes() == boff);

#if boff
  assert(BR());
//  if (!(_n[BR()].n_())) 
  if (_sim->is_first_expand()) {
    _n[BR()].new_model_node( long_label() + ".br", this);
    trace2("newmodelnode", long_label(), _n[BR()].t_());
  }else{ untested();
    trace3("no newmodelnode", BR(), long_label(), _n[BR()].t_());
  }
#endif
  COMPONENT::expand();

  assert(_n[BR()].t_() != INVALID_NODE);
  assert(_n[OUT1].t_() != INVALID_NODE);
  assert(_n[OUT2].t_() != INVALID_NODE);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
DEV_CPOLY_V p4;
DISPATCHER<CARD>::INSTALL d4(&device_dispatcher, "adms_v", &p4);
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet:
