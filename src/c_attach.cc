/*                               -*- C++ -*-
 * Copyright (C) 2007 Albert Davis
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
 */
#include "e_cardlist.h"
#include "c_comand.h"
#include "constant.h"
#include "globals.h"
#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#include <sys/stat.h>
#include <libgen.h>
#include <sys/types.h>
#include <sys/wait.h>

/*--------------------------------------------------------------------------*/
namespace {

/*--------------------------------------------------------------------------*/
using std::string;
/*--------------------------------------------------------------------------*/
std::map<const std::string, void*> attach_list;
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
class CMD_ATTACH : public CMD {
  static void compile(string& filename, string source);
  static void* do_attach(string filename, int flags, bool force=false);
public:
  void do_it(CS& cmd, CARD_LIST*)
  {
    string make = OS::getenv("GNUCAP_ADMS_MAKE", GNUCAP_ADMS_MAKE);
    cmd >> "load_va";
    string command = string("load")
             + " makefile=" + make
	     + " " + cmd.tail();
    trace1("attach", command);
    CMD::command(command, &CARD_LIST::card_list);
  }
} p1;
DISPATCHER<CMD>::INSTALL d1(&command_dispatcher, "load_va|ahdl_include", &p1);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:noet
