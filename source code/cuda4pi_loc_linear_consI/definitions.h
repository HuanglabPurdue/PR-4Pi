/*!
 * \file definitions.h
 * \author <your name>
 * \date April 20, 2012
 * \brief  If you are defining any constants they better be in this file and only this file.
 */
//#include "PSF_Interp.hpp"

#ifndef DEFINITIONS_H
#define DEFINITIONS_H
#ifndef max
//! not defined in the C standard used by visual studio
#define max(a,b) (((a) > (b)) ? (a) : (b))
#endif
#ifndef min
//! not defined in the C standard used by visual studio
#define min(a,b) (((a) < (b)) ? (a) : (b))
#endif

#define PI 3.141592f
#define SPD 32
#define NP 5
#define NPL 5
#define BSZ 1024;



#define index(m,n) (m*16)*PSFSize+(n)*PSFSize+i
#define Vadex(m,n) (m*10)*PSFSize+(n)*PSFSize+i
#define IndF(m,n) m+(n)*psfNum
#define IndA(m,n)   (m*2)*PSFSize+(n)*PSFSize+i
//#define index(m,n) (m)+(n)*Nfit+8*Nfit*i
//#define Vadex(m, n) (m)+(n)*Nfit + 5 * Nfit*i
//#define IndF(m,n) (m)*channelNum+n
#endif


