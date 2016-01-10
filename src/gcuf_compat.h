// gnucap-uf transition
#ifndef GCUF_COMPAT_H
#define GCUF_COMPAT_H

#ifndef USE
# define USE(x)
#endif
/*--------------------------------------------------------------------------*/
#ifndef HAVE_UINT_T
typedef int uint_t;
#endif
/*--------------------------------------------------------------------------*/
#ifndef HAVE_DOUBLE_TYPES
typedef double voltage_t;
typedef double current_t;
typedef double charge_t;
typedef double conductance_t;
#endif
/*--------------------------------------------------------------------------*/
#if !(defined(HAVE_METHOD)) && !(defined(E_STORAGE_H))
enum METHOD {mTRAPGEAR, mEULER, mTRAP, mGEAR, mTRAPEULER};
#endif
/*--------------------------------------------------------------------------*/
#ifndef HAVE_IS_NUMBER
inline bool is_number(long double x){
	const double inf = std::numeric_limits<float>::infinity( );
	return (( x != inf ) && (x != -inf ) && (x == x)) ;
}
#endif
#ifndef trace6
# define trace6(a,b,c,d,e,f,g)
#endif
/*--------------------------------------------------------------------------*/
#ifndef HAVE_GETENV_DEFAULTS
  template<class T>
  T getenv(const std::string&, T def);

  template<class T>
  inline T getenv(const std::string& s, T def) {
    char* ev = ::getenv(s.c_str());
    if(!ev){
      return def;
    }else{
      return T(ev);
    }
  }

  template<>
  inline bool getenv(const std::string& s, bool def) {
    char* ev = ::getenv(s.c_str());
    if(!ev){ untested();
      return def;
    } else { untested();
      return strcmp(ev,"no") && strcmp(ev,"0");
    }
  }
#else
using OS::getenv;
#endif
/*--------------------------------------------------------------------------*/
#ifndef PATCHLEVEL
#define PATCHLEVEL "maybe-upstream"
#endif
/*--------------------------------------------------------------------------*/
#ifndef GCUF_CONST
#define GCUF_CONST
#endif

#endif
