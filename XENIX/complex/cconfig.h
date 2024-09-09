#ifndef _CCONFIG_H
#define _CCONFIG_H

#ifdef _USE_ASM_
#define __FASM__
#endif

#ifdef _USE_BRANCH_
#define __BRANCH__
#endif

#ifdef _DEBUG_
#define __DEBUG_CI__
#undef __DEBUG_CI_D__
#endif

#ifdef __WATCOMC__
#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif
#ifndef M_PI_2
#define M_PI_2		1.57079632679489661923	/* pi/2 */
#endif
#ifndef M_PI_4
#define M_PI_4		0.78539816339744830962	/* pi/4 */
#endif
#ifndef M_LN2
#define M_LN2		0.69314718055994530942	/* log_e 2 */
#endif
#ifndef M_LN10
#define M_LN10		2.30258509299404568402	/* log_e 10 */
#define
#endif

#endif
