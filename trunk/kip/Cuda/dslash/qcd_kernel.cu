
#ifndef _QCD_KERNEL_H_
#define _QCD_KERNEL_H_

#include "qcd.h"

texture<float4, 1, cudaReadModeElementType> spinor;
texture<float4, 1, cudaReadModeElementType> gauge;


__global__ void
testKernel(float4* g_out, int *index) {
    #include "qcd_core.cu"
}

#endif // #ifndef _QCD_KERNEL_H_
