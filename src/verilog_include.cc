/*                         -*- C++ -*-
 * Copyright (C) 2016 Felix Salfelder
 * Author: Felix Salfelder
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
 */
#include <globals.h>
#include <c_comand.h>
#include <u_lang.h>

namespace{

class CMD_INCLUDE : public CMD {
public:
  void do_it(CS& cmd, CARD_LIST* Scope)
  {untested();
    trace1("`include", cmd.tail());
    //    getmerge(cmd, NO_HEADER, Scope);
    std::string file_name, section_name;
    cmd >> file_name;
    CS file(CS::_INC_FILE, file_name);
    LANGUAGE* lang=OPT::language;

    try { untested();
      for (;;) { untested();
	file.get_line("");
	  //skip_pre_stuff(cmd);
	lang->new__instance(file, NULL, Scope);
      }
    }catch (Exception_End_Of_Input& e) { untested();
    }
  }
} p3;
DISPATCHER<CMD>::INSTALL d3(&command_dispatcher, "`include", &p3);
}

// vim:ts=8:sw=2:noet
