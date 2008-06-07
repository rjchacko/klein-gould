#include <stdlib.h>
#include <stdio.h>
#include <cutil.h>
#include "qcd.h"

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (Nh/BLOCK_DIM) // there are Nh threads in total

#define SPINOR_SIZE (24) // spinors have 4*3*2 floats
#define PACKED_GAUGE_SIZE (4*20) // gauge matrices rounded up to fit float4 elements
#define SPINOR_BYTES (SPINOR_SIZE*sizeof(float))
#define PACKED_GAUGE_BYTES (PACKED_GAUGE_SIZE*sizeof(float))


// ----------------------------------------------------------------------
// Cuda code

texture<float4, 1, cudaReadModeElementType> gauge0Tex;
texture<float4, 1, cudaReadModeElementType> gauge1Tex;
texture<float4, 1, cudaReadModeElementType> spinorTex;


__global__ void
dslashKernel(float4* g_out, int oddBit) {
    #include "dslash_core.cu"
}

__global__ void
dslashDaggerKernel(float4* g_out, int oddBit) {
    #include "dslash_dagger_core.cu"
}


// ----------------------------------------------------------------------

void packGaugeField(float4 *res, float **gauge, int oddBit) {
    for (int dir = 0; dir < 4; dir++) {
        float *g = gauge[dir] + oddBit*(Nh*3*3*2);
        for (int i = 0; i < Nh; i++) {
            for (int j = 0; j < 5; j++) {
                float a1, a2, a3=0, a4=0;
                a1 = g[i*18 + j*4 + 0];
                a2 = g[i*18 + j*4 + 1];
                if (j < 4) {
                    a3 = g[i*18 + j*4 + 2];
                    a4 = g[i*18 + j*4 + 3];
                }
                float4 f4 = {a1, a2, a3, a4};
                res[(dir*5+j)*Nh + i] = f4;
            }
        }
    }
}

void packSpinorField(float4 *res, float *spinor, int oddBit) {
    float *s = spinor + oddBit*(Nh*4*3*2);
    for (int i = 0; i < Nh; i++) {
        for (int j = 0; j < 6; j++) {
            float a1 = s[i*(6*4) + j*(4) + 0];
            float a2 = s[i*(6*4) + j*(4) + 1];
            float a3 = s[i*(6*4) + j*(4) + 2];
            float a4 = s[i*(6*4) + j*(4) + 3];
            float4 f4 = {a1, a2, a3, a4};
            res[j*Nh + i] = f4;
        }
    }
}

void unpackSpinorField(float *res, float4 *spinorPacked, int oddBit) {
    float *r = res + oddBit*(Nh*4*3*2);
    
    for (int i = 0; i < Nh; i++) {
        for (int j = 0; j < 6; j++) {
            float4 f4 = spinorPacked[j*Nh + i];
            r[i*(6*4) + j*(4) + 0] = f4.x;
            r[i*(6*4) + j*(4) + 1] = f4.y;
            r[i*(6*4) + j*(4) + 2] = f4.z;
            r[i*(6*4) + j*(4) + 3] = f4.w;
        }
    }
}


CudaFullGauge loadGaugeField(float **gauge) {
    CudaFullGauge ret;
    CUDA_SAFE_CALL(cudaMalloc((void **)&ret.even, Nh*PACKED_GAUGE_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void **)&ret.odd,  Nh*PACKED_GAUGE_BYTES));

    float4 *packedEven = (float4*) malloc(Nh*PACKED_GAUGE_BYTES);
    float4 *packedOdd  = (float4*) malloc(Nh*PACKED_GAUGE_BYTES);
    packGaugeField(packedEven, gauge, 0);
    packGaugeField(packedOdd,  gauge, 1);
    CUDA_SAFE_CALL(cudaMemcpy(ret.even, packedEven, Nh*PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(ret.odd,  packedOdd,  Nh*PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice));    
    free(packedEven);
    free(packedOdd);
    
    return ret;
}

CudaFullSpinor loadSpinorField(float *spinor) {
    CudaFullSpinor ret;
    CUDA_SAFE_CALL(cudaMalloc((void**)&ret.even, Nh*SPINOR_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**)&ret.odd,  Nh*SPINOR_BYTES));
    
    float4 *packedEven = (float4*) malloc(Nh*SPINOR_BYTES);
    float4 *packedOdd  = (float4*) malloc(Nh*SPINOR_BYTES);
    packSpinorField(packedEven, spinor, 0);
    packSpinorField(packedOdd,  spinor, 1);
    CUDA_SAFE_CALL(cudaMemcpy(ret.even, packedEven, Nh*SPINOR_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(ret.odd,  packedOdd,  Nh*SPINOR_BYTES, cudaMemcpyHostToDevice));    
    free(packedEven);
    free(packedOdd);
    
    return ret;
}

void freeGaugeField(CudaFullGauge gauge) {
    CUDA_SAFE_CALL(cudaFree(gauge.even));
    CUDA_SAFE_CALL(cudaFree(gauge.odd));
}

void freeSpinorField(CudaFullSpinor spinor) {
    CUDA_SAFE_CALL(cudaFree(spinor.even));
    CUDA_SAFE_CALL(cudaFree(spinor.odd));
}

void retrieveParitySpinor(float *res, CudaPSpinor spinor) {
    float4 *packed = (float4*) malloc(Nh*SPINOR_BYTES);
    CUDA_SAFE_CALL(cudaMemcpy(packed,  spinor,  Nh*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(res, packed, 0);
    free(packed);
}

void retrieveSpinorField(float *res, CudaFullSpinor spinor) {
    float4 *packedEven = (float4*) malloc(Nh*SPINOR_BYTES);
    float4 *packedOdd = (float4*) malloc(Nh*SPINOR_BYTES);
    CUDA_SAFE_CALL(cudaMemcpy(packedEven, spinor.even, Nh*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    CUDA_SAFE_CALL(cudaMemcpy(packedOdd,  spinor.odd,  Nh*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(res, packedEven, 0);
    unpackSpinorField(res, packedEven, 1);
    free(packedEven);
    free(packedOdd);
}

void compareParitySpinors(float *spinor1, CudaPSpinor cudaSpinor2) {
    float *spinor2 = (float *) malloc(Nh*SPINOR_BYTES);
    retrieveParitySpinor(spinor2, cudaSpinor2);
    CUTBoolean res = cutComparefe(spinor1, spinor2, Nh*SPINOR_BYTES, 1e-4);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
}

void dslashCuda(CudaPSpinor res, CudaFullGauge gauge, CudaPSpinor spinor, int oddBit, int daggerBit) {
    if (oddBit) {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.odd, Nh*PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.even, Nh*PACKED_GAUGE_BYTES); 
    }
    else {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.even, Nh*PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.odd, Nh*PACKED_GAUGE_BYTES); 
    }
    cudaBindTexture(0 /*offset*/, spinorTex, spinor, Nh*SPINOR_BYTES); 

    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    
    if (!daggerBit) {
        dslashKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit);
    }
    else {
        dslashDaggerKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit);
    }
    
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
}

void printCudaDslashInfo() {
    printf("Spinors: %d\n", Nh);
    printf("Global kb: %f\n", N*(PACKED_GAUGE_BYTES+SPINOR_BYTES)/1024.);
    printf("Shared kb: %fkB\n", SHARED_BYTES/1024.);
}
