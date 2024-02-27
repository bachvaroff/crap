#ifndef _COMPLEX_
#define _COMPLEX_

typedef struct {
	double re;
	double im;
} complex_t;

#define Re(Z)		   ((Z).re)
#define Im(Z)		   ((Z).im)

#define printc(Z) do { \
	printf("[ %.16lf, %.16lf ]", Re(Z), Im(Z)); \
} while (0)

#define printcerr(Z) do { \
	fprintf(stderr, "[ %.16lf, %.16lf ]", Re(Z), Im(Z)); \
} while (0)

#define scale(Z, R) do { \
	Re(Z) *= (R); \
	Im(Z) *= (R); \
} while (0)

#define scale2(Z, U, R) { \
	Re(Z) = Re(U) * (R); \
	Im(Z) = Im(U) * (R); \
} while (0)

#define con(Z) do { \
	Im(Z) = -Im(Z); \
} while (0)

#define con2(Z, U) do { \
	Re(Z) = Re(U); \
	Im(Z) = -Im(U); \
} while (0)

#define mag(Z) (sqrt(Re(Z) * Re(Z) + Im(Z) * Im(Z)))

#define magsq(Z) (Re(Z) * Re(Z) + Im(Z) * Im(Z))

#define Arg(Z) (atan2(Im(Z), Re(Z)))

#define dist(U, V) \
	(sqrt((Re(U) - Re(V)) * (Re(U) - Re(V)) + \
	(Im(U) - Im(V)) * (Im(U) - Im(V))))

#define distsq(U, V) \
	((Re(U) - Re(V)) * (Re(U) - Re(V)) + \
	(Im(U) - Im(V)) * (Im(U) - Im(V)))

#define add2(Z, U) do { \
	Re(Z) += Re(U); \
	Im(Z) += Im(U); \
} while (0)

#define add3(Z, U, V) do { \
	Re(Z) = Re(U) + Re(V); \
	Im(Z) = Im(U) + Im(V); \
} while (0)

#define sub2(Z, U) do { \
	Re(Z) -= Re(U); \
	Im(Z) -= Im(U); \
} while (0)

#define sub3(Z, U, V) do { \
	Re(Z) = Re(U) - Re(V); \
	Im(Z) = Im(U) - Im(V); \
} while (0)

#define pseudosc(U, V) (Re(U) * Re(V) + Im(U) * Im(V))

#define mul2(Z, U) do { \
	register double __t__ = Re(Z); \
	Re(Z) = (__t__ * Re(U) - Im(Z) * Im(U)); \
	Im(Z) = (Im(Z) * Re(U) + __t__ * Im(U)); \
} while (0)

#define mul3(Z, U, V) do { \
	Re(Z) = (Re(U) * Re(V) - Im(U) * Im(V)); \
	Im(Z) = (Im(U) * Re(V) + Re(U) * Im(V)); \
} while (0)

#define madd(Z, U, V) do { \
	Re(Z) += (Re(U) * Re(V) - Im(U) * Im(V)); \
	Im(Z) += (Im(U) * Re(V) + Re(U) * Im(V)); \
} while (0)

#define div2(Z, U) do { \
	register double __t__ = Re(Z); \
	Re(Z) = (__t__ * Re(U) + Im(Z) * Im(U)) / magsq(U); \
	Im(Z) = (Im(Z) * Re(U) - __t__ * Im(U)) / magsq(U); \
} while (0)

#define div3(Z, U, V) do { \
	Re(Z) = (Re(U) * Re(V) + Im(U) * Im(V)) / magsq(V); \
	Im(Z) = (Im(U) * Re(V) - Re(U) * Im(V)) / magsq(V); \
} while (0)

#define expiphi(Z, PHI) do { \
	Re(Z) = cos(PHI); \
	Im(Z) = sin(PHI); \
} while (0)

#define expiZ(Z, U) do { \
	Re(Z) = exp(-Im(U)) * cos(Re(U)); \
	Im(Z) = exp(-Im(U)) * sin(Re(U)); \
} while (0)

#define expZ(Z) do { \
	register double __t__ = Re(Z); \
	Re(Z) = exp(__t__) * cos(Im(Z)); \
	Im(Z) = exp(__t__) * sin(Im(Z)); \
} while (0)

#define expZ2(Z, U) do { \
	Re(Z) = exp(Re(U)) * cos(Im(U)); \
	Im(Z) = exp(Re(U)) * sin(Im(U)); \
} while (0)

#define LogZ(Z) do { \
	register double __Arg__ = Arg(Z); \
	Re(Z) = 0.5 * log(magsq(Z)); \
	Im(Z) = __Arg__; \
} while (0)

#define LogZ2(Z, U) do { \
	Re(Z) = 0.5 * log(magsq(U)); \
	Im(Z) = Arg(U); \
} while (0)

#define norm(Z) do { \
	register double __Arg__ = Arg(Z); \
	expiphi(Z, __Arg__); \
} while (0)

#define _norm(Z) do { \
	register double __t__ = mag(Z); \
	Re(Z) /= __t__; \
	Im(Z) /= __t__; \
} while (0)

#define norm2(Z, U) do { \
	register double __Arg__ = Arg(U); \
	expiphi(Z, __Arg__); \
} while (0)

#define _norm2(Z, U) do { \
	register double __t__ = mag(U); \
	Re(Z) = Re(U) / __t__; \
	Im(Z) = Im(U) / __t__; \
} while (0)


#define SinZ(Z) do { \
	register double __t0__ = exp(Im(Z)); \
	register double __t1__ = Re(Z); \
	Re(Z) = 0.5 * (1.0 / __t0__ + __t0__) * sin(__t1__); \
	Im(Z) = -0.5 * (1.0 / __t0__ - __t0__) * cos(__t1__); \
} while (0)

#define SinZ2(Z, U) do { \
	register double __t0__ = exp(Im(U)); \
	register double __t1__ = Re(U); \
	Re(Z) = 0.5 * (1.0 / __t0__ + __t0__) * sin(__t1__); \
	Im(Z) = -0.5 * (1.0 / __t0__ - __t0__) * cos(__t1__); \
} while (0)

#define CosZ(Z) do { \
	register double __t0__ = exp(Im(Z)); \
	register double __t1__ = Re(Z); \
	Re(Z) = 0.5 * (1.0 / __t0__ + __t0__) * cos(__t1__); \
	Im(Z) = 0.5 * (1.0 / __t0__ - __t0__) * sin(__t1__); \
} while (0)

#define CosZ2(Z, U) do { \
	register double __t0__ = exp(Im(U)); \
	register double __t1__ = Re(U); \
	Re(Z) = 0.5 * (1.0 / __t0__ + __t0__) * cos(__t1__); \
	Im(Z) = 0.5 * (1.0 / __t0__ - __t0__) * sin(__t1__); \
} while (0)

#endif
