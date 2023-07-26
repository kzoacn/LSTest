#ifndef CONSTANT_H_
#define CONSTANT_H_

#define N 128

#if N==128
    #define LOGD 8
#else
    #define LOGD 16
#endif



#define T (2*LOGD)
#define FREE_VARS (2*LOGD)
#define QUAD_VARS (FREE_VARS+FREE_VARS*(FREE_VARS-1)/2)

#endif