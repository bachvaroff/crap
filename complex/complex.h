#ifndef _COMPLEX_
#define _COMPLEX_

#include <stdio.h>
#include <math.h>

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
#define mk0(Z) mkZ(Z, 0.0, 0.0)
#define mk1(Z) mkZ(Z, 1.0, 0.0)
#define mki(Z) mkZ(Z, 0.0, 1.0)

#define printZ(PFX, Z, SFX) do { \
	(void)fprintf(stdout, PFX "[ %.16lf, %.16lf ]" SFX, Re(Z), Im(Z)); \
	(void)fflush(stdout); \
} while (0)

#define printZerr(PFX, Z, SFX) do { \
	(void)fprintf(stderr, PFX "[ %.16lf, %.16lf ]" SFX, Re(Z), Im(Z)); \
	(void)fflush(stderr); \
} while (0)

#define dist(U, V) \
	sqrt((Re(U) - Re(V)) * (Re(U) - Re(V)) + (Im(U) - Im(V)) * (Im(U) - Im(V)))
#define distsq(U, V) \
	((Re(U) - Re(V)) * (Re(U) - Re(V)) + (Im(U) - Im(V)) * (Im(U) - Im(V)))
#define pseudosc(U, V) (Re(U) * Re(V) + Im(U) * Im(V))
#define magsq(Z) pseudosc(Z, Z)
#define mag(Z) sqrt(pseudosc(Z, Z))
#define arg(Z, BR) (atan2(Im(Z), Re(Z)) + 2.0 * (double)(BR) * M_PI)

#ifdef __BRANCH
#define Arg(Z) arg(Z, 0)
#else
#define Arg(Z) atan2(Im(Z), Re(Z))
#endif

#define xchgZ(U, V) do { \
	complex_t __t__; \
	__t__ = (U); \
	(U) = (V); \
	(V) = __t__; \
} while (0)

#define swapZ2(Z, U) do { \
	register double __ReU__; \
	__ReU__ = Re(U); \
	Re(Z) = Im(U); \
	Im(Z) = __ReU__; \
} while (0)
#define swapZ(Z) swapZ2(Z, Z)

#define scale3(Z, U, R) do { \
	Re(Z) = Re(U) * (R); \
	Im(Z) = Im(U) * (R); \
} while (0)
#define scale2(Z, R) scale3(Z, Z, R)

#define con2(Z, U) do { \
	Re(Z) = Re(U); \
	Im(Z) = -Im(U); \
} while (0)
#define con(Z) con2(Z, Z)

#define neg2(Z, U) do { \
	Re(Z) = -Re(U); \
	Im(Z) = -Im(U); \
} while (0)
#define neg(Z) neg2(Z, Z)

#define add3(Z, U, V) do { \
	Re(Z) = Re(U) + Re(V); \
	Im(Z) = Im(U) + Im(V); \
} while (0)
#define add2(Z, U) add3(Z, Z, U)

#define sub3(Z, U, V) do { \
	Re(Z) = Re(U) - Re(V); \
	Im(Z) = Im(U) - Im(V); \
} while (0)
#define sub2(Z, U) sub3(Z, Z, U)

#define mul3(Z, U, V) do { \
	register double __ReU__ = Re(U); \
	register double __ReV__ = Re(V); \
	Re(Z) = __ReU__ * __ReV__ - Im(U) * Im(V); \
	Im(Z) = Im(U) * __ReV__ + __ReU__ * Im(V); \
} while (0)
#define mul2(Z, U) mul3(Z, Z, U)

#define inner3(Z, U, V) do { \
	register double __ReU__ = Re(U); \
	register double __ReV__ = Re(V); \
	register double __nImV__ = -Im(V); \
	Re(Z) = __ReU__ * __ReV__ - Im(U) * __nImV__; \
	Im(Z) = Im(U) * __ReV__ + __ReU__ * __nImV__; \
} while (0)
#define inner2(Z, U) inner3(Z, Z, U)

#define madd(Z, U, V) do { \
	register double __ReU__ = Re(U); \
	register double __ReV__ = Re(V); \
	Re(Z) += __ReU__ * __ReV__ - Im(U) * Im(V); \
	Im(Z) += Im(U) * __ReV__ + __ReU__ * Im(V); \
} while (0)

#define div3(Z, U, V) do { \
	register double __ReU__ = Re(U); \
	register double __ReV__ = Re(V); \
	register double __magsqV__ = magsq(V); \
	Re(Z) = (__ReU__ * __ReV__ + Im(U) * Im(V)) / __magsqV__; \
	Im(Z) = (Im(U) * __ReV__ - __ReU__ * Im(V)) / __magsqV__; \
} while (0)
#define div2(Z, U) div3(Z, Z, U)

#define expiphi(Z, PHI) do { \
	Re(Z) = cos(PHI); \
	Im(Z) = sin(PHI); \
} while (0)

#define expiZ2(Z, U) do { \
	register double __ReU__ = Re(U); \
	Re(Z) = exp(-Im(U)) * cos(__ReU__); \
	Im(Z) = exp(-Im(U)) * sin(__ReU__); \
} while (0)
#define expiZ(Z) expiZ(Z, Z)

#define expZ2(Z, U) do { \
	register double __ReU__ = Re(U); \
	Re(Z) = exp(__ReU__) * cos(Im(U)); \
	Im(Z) = exp(__ReU__) * sin(Im(U)); \
} while (0)
#define expZ(Z) expZ2(Z, Z)

#define logZ2(Z, U, BR) do { \
	register double __argU__ = arg(U, BR); \
	Re(Z) = 0.5 * log(magsq(U)); \
	Im(Z) = __argU__; \
} while (0)
#define logZ(Z, BR) logZ2(Z, Z, BR)

#ifdef __BRANCH
#define LogZ2(Z, U) logZ2(Z, U, 0)
#define LogZ(Z) logZ(Z, 0)
#else
#define LogZ2(Z, U) do { \
	register double __ArgU__ = Arg(U); \
	Re(Z) = 0.5 * log(magsq(U)); \
	Im(Z) = __ArgU__; \
} while (0)
#define LogZ(Z) LogZ2(Z, Z)
#endif

#define norm2(Z, U) do { \
	register double __ArgU__ = Arg(U); \
	expiphi(Z, __ArgU__); \
} while (0)
#define norm(Z) norm2(Z, Z)

#define SinZ2(Z, U) do { \
	register double __expImU__ = exp(Im(U)); \
	Im(Z) = -0.5 * (1.0 / __expImU__ - __expImU__) * cos(Re(U)); \
	Re(Z) = 0.5 * (1.0 / __expImU__ + __expImU__) * sin(Re(U)); \
} while (0)
#define SinZ(Z) SinZ2(Z, Z)

#define CosZ2(Z, U) do { \
	register double __expImU__ = exp(Im(U)); \
	Im(Z) = 0.5 * (1.0 / __expImU__ - __expImU__) * sin(Re(U)); \
	Re(Z) = 0.5 * (1.0 / __expImU__ + __expImU__) * cos(Re(U)); \
} while (0)
#define CosZ(Z) CosZ2(Z, Z)

#define SinhZ2(Z, U) do { \
	register double __expReU__ = exp(Re(U)); \
	Re(Z) = 0.5 * (__expReU__ - 1.0 / __expReU__) * cos(Im(U)); \
	Im(Z) = 0.5 * (__expReU__ + 1.0 / __expReU__) * sin(Im(U)); \
} while (0)
#define SinhZ(Z) SinhZ2(Z, Z)

#define CoshZ2(Z, U) do { \
	register double __expReU__ = exp(Re(U)); \
	Re(Z) = 0.5 * (__expReU__ + 1.0 / __expReU__) * cos(Im(U)); \
	Im(Z) = 0.5 * (__expReU__ - 1.0 / __expReU__) * sin(Im(U)); \
} while (0)
#define CoshZ(Z) CoshZ2(Z, Z)

#endif
