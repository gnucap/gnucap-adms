#ifndef IO_MATRIX_H_
#define IO_MATRIX_H_

#include "m_matrix.h"
#include "io_.h"
#include <complex>

//using namespace std;

/* template <>

*/

//inline OMSTREAM& operator<<( OMSTREAM& o, const std::complex<double> &c){
//	return (o << c.real() << "+" << c.imag() << "I"  );
//}

template <class T, class S>
inline S& operator<<( S& o, const BSMATRIX<T> &m)
{
  unsigned size=m.size();
  unsigned i,j;
  T x;

  for(i = 1; i <=size ; i++ ){
    for(j = 1; j <=size ; j++ ){
      x = ((m).s(i,j));
      o << x << " ";
    }
    if (i<size ) o << "\n";
  }
  return o;
}

#if 1 // why are these needed??

template <class T>
inline OMSTREAM& operator<<( OMSTREAM& o, const BSMATRIX<T> &m)
{
  std::string s;
  unsigned size=m.size();
  unsigned i,j;
  T x;

  for(i = 1; i <=size ; i++ ){
    for(j = 1; j <=size ; j++ ){
      x = ((m).s(i,j));
      o << x << " ";
    }
    if (i<size ) o << "\n";
  }
  return o;
}

template <class T>
inline std::ostream& operator<<( std::ostream& o, const BSMATRIX<T> &m)
{
  std::string s;
  unsigned size=m.size();
  unsigned i,j;
  T x;

  for(i = 1; i <=size ; i++ ){
    for(j = 1; j <=size ; j++ ){
      x = ((m).s(i,j));
      o << x << " ";
    }
    if (i<size ) o << "\n";
  }
  return o;
}

#if 0
template <class T>
inline T& operator<<( T& o, const needed_t x)
{
// enum needed_t  {nYES, nNO, nFILL};
	static std::string needed_t_names[] = {".", "*", "0"};
	assert(x>=0);
	assert(x<3);
	o << needed_t_names[x];
	return o;
}
#endif

#endif
#endif
