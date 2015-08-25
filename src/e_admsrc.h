/* adms sources 2015 Felix Salfelder
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
 * Base class for adms sources, sources are sources in the verilog sense.
 */

#include <e_elemnt.h>
//#include <io_misc.h>

// -------------------------------------------------------------------------- //
class ADMS_SOURCE : public ELEMENT {
public:
  ADMS_SOURCE() :
   ELEMENT(),
   _values(NULL),
   _old_values(NULL),
   _n_ports(0), // required?
   _n_vports(0),
   _n_iports(0),
   _inputs(NULL),
   _boff(0)
  { untested(); }
  ADMS_SOURCE(const ADMS_SOURCE& p) :
   ELEMENT(p),
   _values(NULL),
   _old_values(NULL),
   _n_ports(p._n_ports),
   _n_vports(p._n_vports),
   _n_iports(p._n_iports),
   _inputs(NULL),
   _boff(0)
  { untested(); }

  ~ADMS_SOURCE(){
//    delete [] _old_values;
    if (matrix_nodes() > NODES_PER_BRANCH) {
      delete [] _n;
    }else{
      // it is part of a base class
    }
  }
public:
  void set_port_by_index_(uint_t index, COMPONENT*);
  void set_parameters_va(const std::string& Label, CARD* Parent,
		      COMMON_COMPONENT* Common, double Value,
		      uint_t state_count, double state[], double old_state[],
		      uint_t node_count, const node_t nodes[],
		      uint_t branchnode_count, COMPONENT* branchnodes[]);
  void set_parameters(const std::string& Label, CARD* Parent,
		      COMMON_COMPONENT* Common, double Value,
		      uint_t state_count, double state[],
		      uint_t node_count, const node_t nodes[]){
	  untested(); incomplete();
	  set_parameters_va(Label, Parent, Common, Value, state_count, state, state /*incomplete*/,
			  node_count, nodes, 0, NULL);
  }
  void expand_last();
  bool do_tr_last();
//  void tr_restore();
  double tr_probe_num(const std::string&) const;

  uint_t max_nodes()const	{return net_nodes();}
  uint_t min_nodes()const	{return net_nodes();}
  // iports not count in MUTUAL_L. (here, they do.)
  uint_t matrix_nodes()const	{return (_n_ports+_n_iports)*2 + _boff;}
  uint_t net_nodes()const	{return _n_ports*2;}
  uint_t ext_nodes()const	{return (_n_ports+_n_iports)*2;}
  uint_t int_nodes()const	{return _boff;}
  void tr_iwant_matrix_active();
  void tr_unload()
  { untested();
    std::fill_n(_values, _n_ports+_n_iports, 0.);
    _m0.c0 = _m0.c1 = 0.;
    _sim->mark_inc_mode_bad();
    tr_load();
  }
/*--------------------------------------------------------------------------*/
protected:
  double*  _values;
  double*  _old_values;
  uint_t _n_ports;
  uint_t _n_vports;
  uint_t _n_iports;
  COMPONENT** _inputs;
  std::vector<COMPONENT*> _iports;
protected:
  unsigned _boff; //a hack
};
/*--------------------------------------------------------------------------*/
/* set current port
 */
inline void ADMS_SOURCE::set_port_by_index_(uint_t n, COMPONENT* p)
{ untested();
  if(n+1>_n_iports){untested();
    assert(net_nodes() == 0); // first time.
    _n_iports = n+1;
    _iports.resize(_n_iports);
  }else{untested();
  }
  _iports[n] = p;
}
/*--------------------------------------------------------------------------*/
/* set: set parameters, used in model building
 */
inline void ADMS_SOURCE::set_parameters_va(const std::string& Label, CARD *Owner,
				 COMMON_COMPONENT *Common, double Value,
				 uint_t n_states, double states[], double old_states[],
				 uint_t n_nodes, const node_t nodes[],
				 uint_t n_currents, COMPONENT* currents[])
{ itested();
  trace4("set_parameters_va", Label, n_states, n_nodes, n_currents);
  bool first_time = (net_nodes() == 0);
  assert(currents || !n_currents);

  set_label(Label);
  set_owner(Owner);
  set_value(Value);
  attach_common(Common);

  if (first_time) { itested();
    _n_ports = n_nodes/2;  // also sets num_nodes() = _n_ports*2
    _n_vports = _n_ports - 1;
    _n_iports = n_currents;
    trace4("ADMS_SOURCE::sp_va ports", _n_ports, _n_vports, _n_iports, n_nodes);
    assert(_n_ports + _n_iports == n_states - 1);

//    assert(!_old_values);
//    _old_values = new double[n_states];
    _inputs = new COMPONENT*[_n_iports];

    if (matrix_nodes() > NODES_PER_BRANCH) {
      // allocate a bigger node list
      _n = new node_t[matrix_nodes()];
    }else{
      // use the default node list, already set
    }
  }else{ itested();
    assert(_n_ports + _n_iports == n_states - 1);// ?
    assert(net_nodes() == n_nodes);
  }
  assert(matrix_nodes() == 2*n_states + _boff - 2); // used by iwant_matrix.
  assert(ext_nodes() == 2*n_states - 2); // used by map_nodes
  trace2("ADMS_SOURCE::sp_va now", _m0, _m1);

  _values = states;
  _old_values = old_states;
  assert(states);
  std::fill_n(_values, n_states, 0.);
  std::fill_n(_old_values, n_states, 0.);
  notstd::copy_n(nodes, 2, _n); // output nodes
  // leave gap for current branch
  assert(net_nodes()>1);
  assert(net_nodes() == _n_ports * 2);

  trace4("ADMS_SOURCE::sp_va DEBUG", _n_ports, _n_vports, net_nodes(), _boff);
  // skip branch node (IN1) if boff
  notstd::copy_n(nodes+2, net_nodes()-2, _n+_boff+2);
  notstd::copy_n(currents, _n_iports, _inputs);

}

void ADMS_SOURCE::expand_last()
{
  trace2("ADMS_SOURCE::expand_last", long_label(), mfactor());
  for(unsigned i=0; i<_n_iports; ++i) { untested();
    assert(_inputs[i]);

    node_t* ni = _n+_boff+2+_n_vports*2;
    assert(net_nodes() == 2+_n_vports*2);

    if (_inputs[i]->has_inode()) {untested();
      trace3("ADMS_SOURCE::expand_last inode", _inputs[i]->long_label(), i, net_nodes());
      ni[2*i] = _inputs[i]->n_(IN1);
      ni[2*i+1].set_to_ground(this);
    }else if (_inputs[i]->has_iv_probe()) { untested();
      trace4("ADMS_SOURCE::expand_last colleting nodes", _inputs[i]->long_label(), i, net_nodes(), INT_MAX);
      ni[2*i]   = _inputs[i]->n_(OUT1);
      ni[2*i+1] = _inputs[i]->n_(OUT2);
      trace2("ADMS_SOURCE::expand_last colleting nodes", ni[2*i].m_(), ni[2*i+1].m_());
    }
  }
  trace3("ADMS_SOURCE::sp_va done", _n_ports, _n_vports, _n_iports);
  assert(_n_vports+1 == _n_ports);
}

/*--------------------------------------------------------------------------*/
// the current controlled current part...
inline void ADMS_SOURCE::tr_iwant_matrix_active()
{ untested();
  assert(matrix_nodes() == 2+_boff+2*(_n_vports + _n_iports));
  assert(is_device());
  assert(!subckt());
  assert(!_boff);

  assert(_n[OUT1].m_() != INVALID_NODE);
  assert(_n[OUT2].m_() != INVALID_NODE);
//  assert(_n[IN1].m_() != INVALID_NODE);
//  assert(_n[IN2].m_() != INVALID_NODE);

  node_t* ni = _n+_boff+2+2*_n_vports;

  //_sim->_aa.iwant(_n[OUT1].m_(),_n[OUT2].m_());
  for(unsigned i=0; i<_n_vports; ++i) {
    if (!_boff){
      _sim->_aa.iwant(_n[OUT1].m_(),_n[2*i+0].m_());
      _sim->_aa.iwant(_n[OUT1].m_(),_n[2*i+1].m_());
      _sim->_lu.iwant(_n[OUT1].m_(),_n[2*i+0].m_());
      _sim->_lu.iwant(_n[OUT1].m_(),_n[2*i+1].m_());

    }
    _sim->_aa.iwant(_n[OUT2+_boff].m_(),_n[2*i+0+_boff].m_());
    _sim->_aa.iwant(_n[OUT2+_boff].m_(),_n[2*i+1+_boff].m_());
    _sim->_lu.iwant(_n[OUT2+_boff].m_(),_n[2*i+0+_boff].m_());
    _sim->_lu.iwant(_n[OUT2+_boff].m_(),_n[2*i+1+_boff].m_());
  }
}
/*--------------------------------------------------------------------------*/

inline bool ADMS_SOURCE::do_tr_last()
{ itested();
  double* ival = _values+2+_n_vports;
  double* ctrl = _values+2+_n_vports+_n_iports;

  for (uint_t i=0; i<_n_iports; ++i) { untested();
    trace3("ADMS_SOURCE::do_tr_last", long_label(), _inputs[i]->long_label(), _inputs[i]->has_iv_probe());
    trace1("ADMS_SOURCE::do_tr_last", ival[i]);
    if (_inputs[i]->has_iv_probe()) { untested();
      ELEMENT const* ie = prechecked_cast<ELEMENT* const>(_inputs[i]);
      assert(ie);
      //cccs do_tr_last.
      // _m0.c1  = _y[0].f1 * (_input->_loss0 + _input->_m0.c1); ?
      ival[i] = ctrl[i] * ie->_m0.c1;
    }else{ itested();
      assert(_inputs[i]->has_inode());
      ival[i] = - ctrl[i];
    }
  }
  trace1("ADMS_SOURCE::do_tr_last", converged());
  return converged();
}
/*--------------------------------------------------------------------------*/

double ADMS_SOURCE::tr_probe_num(const std::string& x) const
{ untested();
  unsigned nval = _n_vports+_n_iports+2;
  std::string* vnames = new std::string[nval];
  for(unsigned i=0; i<nval; i++)
  {
    std::stringstream ss;
    ss << i;
    vnames[i] = "_val" + ss.str();
  }

  // node probes
  for(int i=0; i<2; i++) {
    if( Umatch(x, vnames[i] + ' ') ) {
      delete[] vnames;
      return _values[i];
    }
  }
  delete[] vnames;

  return ELEMENT::tr_probe_num(x);
}

// vim:ts=8:sw=2:noet
