/*                              -*- C++ -*-
 * Copyright (C) 2013 Felix Salfelder
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

#define ADD_VERSION

#include "config.h"
#include <c_comand.h>
#include <d_dot.h>
#include <d_coment.h>
#include <d_subckt.h>
#include <e_model.h>
#include <u_lang.h>
#include <l_lib.h>
#include <globals.h>
extern "C"{
#include "md5.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/wait.h>
}
#include "gcuf_compat.h"
/*--------------------------------------------------------------------------*/
#ifndef GNUCAP_ADMS_IMPLICIT
#define GNUCAP_ADMS_IMPLICIT ""
#endif
/*--------------------------------------------------------------------------*/
namespace { //
/*--------------------------------------------------------------------------*/
using std::string;
using std::vector;
using std::ostream;
/*--------------------------------------------------------------------------*/
class LANG_ADMS : public LANGUAGE { //
public:
  LANG_ADMS();
  std::string name()const {return "adms";}
  bool case_insensitive()const {return false;}
  UNITS units()const {return uSI;}

public: // override virtual, called by commands
  void		parse_top_item(CS&, CARD_LIST*);
//  DEV_COMMENT*	parse_comment(CS&, DEV_COMMENT*);
    DEV_DOT*	parse_command(CS&, DEV_DOT*);
//  MODEL_CARD*	parse_paramset(CS&, MODEL_CARD*);
//  MODEL_SUBCKT* parse_module(CS&, MODEL_SUBCKT*);
//  COMPONENT*	parse_instance(CS&, COMPONENT*);

private: // override virtual, called by print_item
  void print_paramset(OMSTREAM&, const MODEL_CARD*){}
  void print_module(OMSTREAM&, const MODEL_SUBCKT*){}
  void print_instance(OMSTREAM&, const COMPONENT*){}
  void print_comment(OMSTREAM&, const DEV_COMMENT*){}
  void print_command(OMSTREAM& o, const DEV_DOT* c){}
private: // local
  void print_args(OMSTREAM&, const MODEL_CARD*);
  void print_args(OMSTREAM&, const COMPONENT*);
 
private: // unnecessary, make compiler hapopy
  virtual std::string arg_front()const{return "";}
  virtual std::string arg_mid()const{return "";}
  virtual std::string arg_back()const{return "";}
  virtual DEV_COMMENT*	parse_comment(CS&, DEV_COMMENT*){return NULL;}
//  virtual DEV_DOT*	parse_command(CS&, DEV_DOT*){return NULL;}
  virtual MODEL_CARD*	parse_paramset(CS&, MODEL_CARD*){return NULL;}
  virtual MODEL_SUBCKT* parse_module(CS&, MODEL_SUBCKT*){return NULL;}
  virtual COMPONENT*	parse_instance(CS&, COMPONENT*){return NULL;}
  std::string	find_type_in_string(CS&) GCUF_CONST {return "";}

  md5_ctx _md5_ctx;
  ostream* vastream;
  vector<string> _lines;
  string _adms_tmpdir;
  void init();
  void attachadms(vector<string> lines, string p, string f);
  void admsXml(string p, string n);
} lang_adms;

DISPATCHER<LANGUAGE>::INSTALL
	d(&language_dispatcher, lang_adms.name(), &lang_adms);
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
LANG_ADMS::LANG_ADMS()
{
  init();
  string user = getenv("USER", "anonymous");
  _adms_tmpdir = getenv("GNUCAP_ADMS_TMPDIR", "/tmp/gnucap_adms_" + user);
  mkdir(_adms_tmpdir.c_str(), 0777);
}
/*--------------------------------------------------------------------------*/
void LANG_ADMS::init()
{
  _lines.clear();
  _lines.push_back("//" PATCHLEVEL "\n");
  __md5_init_ctx (&_md5_ctx);
}
/*--------------------------------------------------------------------------*/
DEV_DOT* LANG_ADMS::parse_command(CS& cmd, DEV_DOT* x)
{ untested();
  assert(x);
  x->set(cmd.fullstring());
  CARD_LIST* scope = (x->owner()) ? x->owner()->subckt() : &CARD_LIST::card_list;

  cmd.reset();
  CMD::cmdproc(cmd, scope);
  delete x;
  return NULL;
}
/*--------------------------------------------------------------------------*/
void LANG_ADMS::admsXml(string path, string file)
{
  error(bTRACE, "running admsXml" + path + " " + file);
  int childret;
  pid_t p = vfork();
  if (p) {
    waitpid(p, &childret, 0);
  } else {
    string adms_implicit = getenv("GNUCAP_ADMS_IMPLICIT", GNUCAP_ADMS_IMPLICIT);
    if (adms_implicit != "") {
      setenv("adms_implicit_transforms", adms_implicit.c_str(), 1);
    } else { untested();
    }

    string sp = getenv("GNUCAP_ADMS_DATADIR", ADMS_DATADIR);
    string id = getenv("GNUCAP_ADMS_INCLUDEDIR", ADMS_INCLUDEDIR);
    chdir(path.c_str());

    error(bDEBUG, string("calling cd ") + path + "; admsXml "
        + "-I" + id + " "
        + "-e" + sp + "/gnucap_0.xml "
        + "-e" + sp + "/gnucap_1.xml "
        + "-e" + sp + "/gnucap_2.xml "
        + file+".va "
        + "-o"+ file);

    int ret = execlp( "admsXml", "admsXml",
        ("-I" + id).c_str(),
        ("-e" + sp + "/gnucap_0.xml").c_str(),
        ("-e" + sp + "/gnucap_1.xml").c_str(),
        ("-e" + sp + "/gnucap_2.xml").c_str(),
        (file+".va").c_str(),
        "-o", file.c_str(),
        (char*) NULL);
    _exit(ret);
  }
}
/*--------------------------------------------------------------------------*/
void LANG_ADMS::attachadms(vector<string> lines, string path, string file)
{
  trace2("attachadms", path, file);
  struct stat ccattrib;
  int ccstat = stat((path+"/"+file+".va").c_str(), &ccattrib);
  if (ccstat) {
    FILE* f = fopen((path+"/"+file+".va").c_str(),"w");
    for (vector<string>::const_iterator i=lines.begin(); i!=lines.end(); ++i) {
      fprintf(f, "%s\n", i->c_str());
    }
    fclose(f);
    admsXml(path, file);
  } else {
  }

  char* cppf_save = getenv("CPPFLAGS");
  char* libs_save = getenv("LIBS");
  char* ldfl_save = getenv("LDFLAGS");

  string cppf = getenv("GNUCAP_ADMS_CPPFLAGS", string(ADMS_CPPFLAGS));
  string ldfl = getenv("GNUCAP_ADMS_LDFLAGS", string(ADMS_LDFLAGS));
//  string libdir = getenv("GNUCAP_ADMS_LIBDIR", string(ADMS_LIBDIR));
  string libs = getenv("GNUCAP_ADMS_LIBS", string(ADMS_LIBS));

  error(bDEBUG, string("setting CPPFLAGS=") + cppf + "\n");
  error(bDEBUG, string("setting LDFLAGS=") + ldfl + "\n");
  error(bDEBUG, string("setting LIBS=") + libs + "\n");

  // we should append, but probably not necessary.
  setenv("CPPFLAGS", cppf.c_str(), 1);
  setenv("LIBS", libs.c_str(), 1);
  setenv("LDFLAGS", ldfl.c_str(), 1);

  CMD::command("load " + path + "/" + file + ".cc", &CARD_LIST::card_list);

  if (cppf_save) { untested();
    setenv("CPPFLAGS", cppf_save, 1);
  } else {
    unsetenv("CPPFLAGS");
  }
  if (ldfl_save) { untested();
    setenv("LDFLAGS", ldfl_save, 1);
  } else {
    unsetenv("LDFLAGS");
  }
  if (libs_save) { untested();
    setenv("LIBS", libs_save, 1);
  } else {
    unsetenv("LIBS");
  }
}
/*--------------------------------------------------------------------------*/
static const char hexchars[] = "0123456789abcdef";
static string hexfrommd5(unsigned char* md5)
{
  char result[33];
  char* seek = result;

  for (unsigned i=0; i<16; i++) {
    *(seek++) = hexchars[md5[i] >> 4];
    *(seek++) = hexchars[md5[i] & 15];
  }
  *seek = 0;
  seek = result;
  return string(seek);
}
/*--------------------------------------------------------------------------*/
void LANG_ADMS::parse_top_item(CS& cmd, CARD_LIST* Scope)
{
  cmd.get_line("gnucap-adms>");
  if( cmd >> "endadms" ) {
    unsigned char md5[16];
    __md5_finish_ctx (&_md5_ctx, md5);
    attachadms(_lines, _adms_tmpdir, hexfrommd5(md5));
    CMD::command("options lang=verilog", Scope); // switch to previous lang?
    init();
  } else {
    string line = cmd.fullstring(); //  + '\n';
    _lines.push_back(line);
    __md5_process_bytes (line.c_str(), line.size(), &_md5_ctx);
  }
  // hmm maybe later.
  //new__instance(cmd, NULL, Scope);
}
/*--------------------------------------------------------------------------*/
class CMD_ADMS : public CMD { //
public:
  void do_it(CS&, CARD_LIST* Scope)
  {
    command("options lang=adms", Scope);
  }
} p8;
DISPATCHER<CMD>::INSTALL d8(&command_dispatcher, "adms", &p8);
/*--------------------------------------------------------------------------*/
}
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
// vim:ts=8:sw=2:et:
