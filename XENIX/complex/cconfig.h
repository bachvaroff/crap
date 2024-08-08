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

#endif
