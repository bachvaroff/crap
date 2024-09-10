#ifndef _COMPLEX_H
#define _COMPLEX_H

#include <stdio.h>
#include <math.h>

#include "cconfig.h"
#ifdef __FASM__
#include "fasm.h"
#endif

typedef struct _complex_t {
	double re;
	double im;
} complex_t;

#define Re(Z) ((Z).re)
#define Im(Z) ((Z).im)
#define mkZ(Z, RE, IM) do { \
	Re(Z) = (RE); \
	Im(Z) = (IM); \
} while (0)
#define mkZ0(Z) mkZ(Z, 0.0, 0.0)
#define mkZ1(Z) mkZ(Z, 1.0, 0.0)
#define mkZi(Z) mkZ(Z, 0.0, 1.0)

#define printZfp(F, PFX, Z, SFX) do { \
	(void)fprintf((F), PFX "[ %.16lf %.16lf ]" SFX, Re(Z), Im(Z)); \
	(void)fflush((F)); \
} while (0)
#define printZ(PFX, Z, SFX) printZfp(stdout, PFX, Z, SFX)
#define printZerr(PFX, Z, SFX) printZfp(stderr, PFX, Z, SFX)

#define distZ(U, V) \
	sqrt((Re(U) - Re(V)) * (Re(U) - Re(V)) + (Im(U) - Im(V)) * (Im(U) - Im(V)))
#define distsqZ(U, V) \
	((Re(U) - Re(V)) * (Re(U) - Re(V)) + (Im(U) - Im(V)) * (Im(U) - Im(V)))
#define pseudoscZ(U, V) (Re(U) * Re(V) + Im(U) * Im(V))
#define magsqZ(Z) pseudoscZ(Z, Z)
#define magZ(Z) sqrt(pseudoscZ(Z, Z))
#define argZ(Z, BR) (atan2(Im(Z), Re(Z)) + 2.0 * (double)(BR) * M_PI)

#ifdef __BRANCH__
#define ArgZ(Z) argZ(Z, 0)
#else
#define ArgZ(Z) atan2(Im(Z), Re(Z))
#endif

#ifdef __FASM__
#define FSINCOS(PHI, C, S) fsincos((PHI), &(C), &(S))
#else
#define FSINCOS(PHI, C, S) do { \
	double __phi__ = (PHI); \
	(C) = cos(__phi__); \
	(S) = sin(__phi__); \
} while (0)
#endif

#define xchgZ(U, V) do { \
	complex_t __txchg__; \
	__txchg__ = (U); \
	(U) = (V); \
	(V) = __txchg__; \
} while (0)

#define swapZ2(Z, U) do { \
	double __ReU__ = Re(U); \
	Re(Z) = Im(U); \
	Im(Z) = __ReU__; \
} while (0)
#define swapZ(Z) swapZ2(Z, Z)

#define scalemZ3(Z, U, R) do { \
	Re(Z) = Re(U) * (R); \
	Im(Z) = Im(U) * (R); \
} while (0)
#define scalemZ2(Z, R) scalemZ3(Z, Z, R)

#define scaledZ3(Z, U, R) do { \
	Re(Z) = Re(U) / (R); \
	Im(Z) = Im(U) / (R); \
} while (0)
#define scaledZ2(Z, R) scaledZ3(Z, Z, R)

#define conjZ2(Z, U) do { \
	Re(Z) = Re(U); \
	Im(Z) = -Im(U); \
} while (0)
#define conjZ(Z) conjZ2(Z, Z)

#define negZ2(Z, U) do { \
	Re(Z) = -Re(U); \
	Im(Z) = -Im(U); \
} while (0)
#define negZ(Z) negZ2(Z, Z)

#define invZ2(Z, U) do { \
	double __magsqU__ = magsqZ(U); \
	conjZ2(Z, U); \
	scaledZ2(Z, __magsqU__); \
} while (0)
#define invZ(Z) invZ2(Z, Z)

#define addZ3(Z, U, V) do { \
	Re(Z) = Re(U) + Re(V); \
	Im(Z) = Im(U) + Im(V); \
} while (0)
#define addZ2(Z, U) addZ3(Z, Z, U)

#define subZ3(Z, U, V) do { \
	Re(Z) = Re(U) - Re(V); \
	Im(Z) = Im(U) - Im(V); \
} while (0)
#define subZ2(Z, U) subZ3(Z, Z, U)

#define mulZ3(Z, U, V) do { \
	double __ReU__ = Re(U); \
	double __ReV__ = Re(V); \
	Re(Z) = __ReU__ * __ReV__ - Im(U) * Im(V); \
	Im(Z) = Im(U) * __ReV__ + __ReU__ * Im(V); \
} while (0)
#define mulZ2(Z, U) mulZ3(Z, Z, U)

#define mulpZ3(Z, U, V) do { \
	Re(Z) = Re(U) * Re(V); \
	Im(Z) = Im(U) * Im(V); \
} while (0)
#define mulpZ2(Z, U) mulpZ3(Z, Z, U)

#define innerZ3(Z, U, V) do { \
	double __ReU__ = Re(U); \
	double __ReV__ = Re(V); \
	Re(Z) = __ReU__ * __ReV__ + Im(U) * Im(V); \
	Im(Z) = Im(U) * __ReV__ - __ReU__ * Im(V); \
} while (0)
#define innerZ2(Z, U) innerZ3(Z, Z, U)

#define maddZ(Z, U, V) do { \
	double __ReU__ = Re(U); \
	double __ReV__ = Re(V); \
	Re(Z) += __ReU__ * __ReV__ - Im(U) * Im(V); \
	Im(Z) += Im(U) * __ReV__ + __ReU__ * Im(V); \
} while (0)

#define divZ3(Z, U, V) do { \
	double __ReU__ = Re(U); \
	double __ReV__ = Re(V); \
	double __magsqV__ = magsqZ(V); \
	Re(Z) = (__ReU__ * __ReV__ + Im(U) * Im(V)) / __magsqV__; \
	Im(Z) = (Im(U) * __ReV__ - __ReU__ * Im(V)) / __magsqV__; \
} while (0)
#define divZ2(Z, U) divZ3(Z, Z, U)

#define expiphiZ(Z, PHI) FSINCOS((PHI), Re(Z), Im(Z))
#define normZ2(Z, U) expiphiZ(Z, ArgZ(U))
#define normZ(Z) normZ2(Z, Z)

#define expiZ2(Z, U) do { \
	double __expImU__ = exp(Im(U)); \
	FSINCOS(Re(U), Re(Z), Im(Z)); \
	scaledZ2(Z, __expImU__); \
} while (0)
#define expiZ(Z) expiZ2(Z, Z)

#define expZ2(Z, U) do { \
	double __expReU__ = exp(Re(U)); \
	FSINCOS(Im(U), Re(Z), Im(Z)); \
	scalemZ2(Z, __expReU__); \
} while (0)
#define expZ(Z) expZ2(Z, Z)

#define logZ2(Z, U, BR) do { \
	double __argU__ = argZ(U, BR); \
	Re(Z) = 0.5 * log(magsqZ(U)); \
	Im(Z) = __argU__; \
} while (0)
#define logZ(Z, BR) logZ2(Z, Z, BR)

#ifdef __BRANCH__
#define LogZ2(Z, U) logZ2(Z, U, 0)
#define LogZ(Z) logZ(Z, 0)
#else
#define LogZ2(Z, U) do { \
	double __ArgU__ = ArgZ(U); \
	Re(Z) = 0.5 * log(magsqZ(U)); \
	Im(Z) = __ArgU__; \
} while (0)
#define LogZ(Z) LogZ2(Z, Z)
#endif

#define SinZ2(Z, U) do { \
	double __expImU__ = exp(Im(U)); \
	double __expnImU__ = 1.0 / __expImU__; \
	FSINCOS(Re(U), Im(Z), Re(Z)); \
	Re(Z) *= 0.5 * (__expImU__ + __expnImU__); \
	Im(Z) *= 0.5 * (__expImU__ - __expnImU__); \
} while (0)
#define SinZ(Z) SinZ2(Z, Z)

#define CosZ2(Z, U) do { \
	double __expImU__ = exp(Im(U)); \
	double __expnImU__ = 1.0 / __expImU__; \
	FSINCOS(Re(U), Re(Z), Im(Z)); \
	Re(Z) *= 0.5 * (__expnImU__ + __expImU__); \
	Im(Z) *= 0.5 * (__expnImU__ - __expImU__); \
} while (0)
#define CosZ(Z) CosZ2(Z, Z)

#define SinhZ2(Z, U) do { \
	double __expReU__ = exp(Re(U)); \
	double __expnReU__ = 1.0 / __expReU__; \
	FSINCOS(Im(U), Re(Z), Im(Z)); \
	Re(Z) *= 0.5 * (__expReU__ - __expnReU__); \
	Im(Z) *= 0.5 * (__expReU__ + __expnReU__); \
} while (0)
#define SinhZ(Z) SinhZ2(Z, Z)

#define CoshZ2(Z, U) do { \
	double __expReU__ = exp(Re(U)); \
	double __expnReU__ = 1.0 / __expReU__; \
	FSINCOS(Im(U), Re(Z), Im(Z)); \
	Re(Z) *= 0.5 * (__expReU__ + __expnReU__); \
	Im(Z) *= 0.5 * (__expReU__ - __expnReU__); \
} while (0)
#define CoshZ(Z) CoshZ2(Z, Z)

#endif
