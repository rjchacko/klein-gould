#include <stdlib.h>
#include <stdio.h>
#include "qcd.h"

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (Nh/BLOCK_DIM) // there are Nh threads in total

#define SPINOR_BYTES (Nh*spinorSiteSize*sizeof(float))
#define PACKED_GAUGE_BYTES (4*Nh*packedGaugeSiteSize*sizeof(float))


// ----------------------------------------------------------------------
// Cuda code

texture<float4, 1, cudaReadModeElementType> gauge0Tex;
texture<float4, 1, cudaReadModeElementType> gauge1Tex;
texture<float4, 1, cudaReadModeElementType> spinorTex;
texture<float4, 1, cudaReadModeElementType> accumTex;

//#define WRITE_FLOAT4
//#define WRITE_FLOAT1_SMEM
#define WRITE_FLOAT1_STAGGERED

__global__ void
dslashKernel(float4* g_out, int oddBit) {
    #include "dslash_core.cu"
}

__global__ void
dslashDaggerKernel(float4* g_out, int oddBit) {
    #include "dslash_dagger_core.cu"
}


#define DSLASH_XPAY

__global__ void
dslashXpayKernel(float4* g_out, int oddBit, float a) {
    #include "dslash_core.cu"
}

__global__ void
dslashDaggerXpayKernel(float4* g_out, int oddBit, float a) {
    #include "dslash_dagger_core.cu"
}

#undef DSLASH_XPAY


// ----------------------------------------------------------------------

void packGaugeField(float4 *res, float **gauge, int oddBit) {
    for (int dir = 0; dir < 4; dir++) {
        float *g = gauge[dir] + oddBit*Nh*gaugeSiteSize;
        for (int i = 0; i < Nh; i++) {
            for (int j = 0; j < 3; j++) {
                float a1 = g[i*18 + j*4 + 0];
                float a2 = g[i*18 + j*4 + 1];
                float a3 = g[i*18 + j*4 + 2];
                float a4 = g[i*18 + j*4 + 3];
                float4 f4 = {a1, a2, a3, a4};
                res[(dir*3+j)*Nh + i] = f4;
            }
        }
    }
}

void packParitySpinor(float4 *res, float *spinor) {
    for (int i = 0; i < Nh; i++) {
        for (int j = 0; j < 6; j++) {
            float a1 = spinor[i*(6*4) + j*(4) + 0];
            float a2 = spinor[i*(6*4) + j*(4) + 1];
            float a3 = spinor[i*(6*4) + j*(4) + 2];
            float a4 = spinor[i*(6*4) + j*(4) + 3];
            float4 f4 = {a1, a2, a3, a4};
            res[j*Nh + i] = f4;
        }
    }
}

void unpackParitySpinor(float *res, float4 *spinorPacked) {
    for (int i = 0; i < Nh; i++) {
        for (int j = 0; j < 6; j++) {
            float4 f4 = spinorPacked[j*Nh + i];
            res[i*(6*4) + j*(4) + 0] = f4.x;
            res[i*(6*4) + j*(4) + 1] = f4.y;
            res[i*(6*4) + j*(4) + 2] = f4.z;
            res[i*(6*4) + j*(4) + 3] = f4.w;
        }
    }
}

CudaPSpinor allocateParitySpinor() {
    CudaPSpinor ret;
    cudaMalloc((void**)&ret, SPINOR_BYTES);
    return ret;
}

CudaFullGauge loadGaugeField(float **gauge) {
    CudaFullGauge ret;
    cudaMalloc((void **)&ret.even, PACKED_GAUGE_BYTES);
    cudaMalloc((void **)&ret.odd,  PACKED_GAUGE_BYTES);

    float4 *packedEven = (float4*) malloc(PACKED_GAUGE_BYTES);
    float4 *packedOdd  = (float4*) malloc(PACKED_GAUGE_BYTES);
    packGaugeField(packedEven, gauge, 0);
    packGaugeField(packedOdd,  gauge, 1);
    cudaMemcpy(ret.even, packedEven, PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice);
    cudaMemcpy(ret.odd,  packedOdd,  PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice);    
    free(packedEven);
    free(packedOdd);
    
    return ret;
}

CudaPSpinor loadParitySpinor(float *spinor) {
    CudaPSpinor ret;
    cudaMalloc((void**)&ret, SPINOR_BYTES);
    
    float4 *packed = (float4*) malloc(SPINOR_BYTES);
    packParitySpinor(packed, spinor);
    cudaMemcpy(ret, packed, SPINOR_BYTES, cudaMemcpyHostToDevice);
    free(packed);
    
    return ret;
}

CudaFullSpinor loadSpinorField(float *spinor) {
    CudaFullSpinor ret;
    ret.even = loadParitySpinor(spinor);
    ret.odd  = loadParitySpinor(spinor + Nh*spinorSiteSize);
    return ret;
}

void freeParitySpinor(CudaPSpinor spinor) {
    cudaFree(spinor);
}

void freeGaugeField(CudaFullGauge gauge) {
    cudaFree(gauge.even);
    cudaFree(gauge.odd);
}

void freeSpinorField(CudaFullSpinor spinor) {
    cudaFree(spinor.even);
    cudaFree(spinor.odd);
}

void retrieveParitySpinor(float *res, CudaPSpinor spinor) {
    float4 *packed = (float4*) malloc(SPINOR_BYTES);
    cudaMemcpy(packed, spinor, SPINOR_BYTES, cudaMemcpyDeviceToHost);
    unpackParitySpinor(res, packed);
    free(packed);
}

void retrieveSpinorField(float *res, CudaFullSpinor spinor) {
    retrieveParitySpinor(res, spinor.even);
    retrieveParitySpinor(res+Nh*spinorSiteSize, spinor.odd);
}

void dslashCuda(CudaPSpinor res, CudaFullGauge gauge, CudaPSpinor spinor, int oddBit, int daggerBit) {
    if (oddBit) {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.odd, PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.even, PACKED_GAUGE_BYTES); 
    }
    else {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.even, PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.odd, PACKED_GAUGE_BYTES); 
    }
    cudaBindTexture(0 /*offset*/, spinorTex, spinor, SPINOR_BYTES); 

    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    
    if (!daggerBit) {
        dslashKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit);
    }
    else {
        dslashDaggerKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit);
    }
}

void dslashXpayCuda(CudaPSpinor res, CudaFullGauge gauge, CudaPSpinor spinor, int oddBit, int daggerBit, CudaPSpinor x, float a) {
    if (oddBit) {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.odd, PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.even, PACKED_GAUGE_BYTES); 
    }
    else {
        cudaBindTexture(0 /*offset*/, gauge0Tex, gauge.even, PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, gauge.odd, PACKED_GAUGE_BYTES); 
    }
    cudaBindTexture(0 /*offset*/, spinorTex, spinor, SPINOR_BYTES); 
    cudaBindTexture(0 /*offset*/, accumTex, x, SPINOR_BYTES); 
    
    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    
    if (!daggerBit) {
        dslashXpayKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit, a);
    }
    else {
        dslashDaggerXpayKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)res, oddBit, a);
    }
}

int dslashCudaSharedBytes() {
    return SHARED_BYTES;
}

// Apply the even-odd preconditioned Dirac operator
void MatPCCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp) {
    float kappa2 = -kappa*kappa;
    dslashCuda(tmp, gauge, inEven, 1, 0);
    dslashXpayCuda(outEven, gauge, tmp, 0, 0, inEven, kappa2); 
}

// Apply the even-odd preconditioned Dirac operator
void MatPCDagCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp) {
    float kappa2 = -kappa*kappa;
    dslashCuda(tmp, gauge, inEven, 1, 1);
    dslashXpayCuda(outEven, gauge, tmp, 0, 1, inEven, kappa2);
}

void MatPCDagMatPCCuda(CudaPSpinor outEven, CudaFullGauge gauge, CudaPSpinor inEven, float kappa, CudaPSpinor tmp1, CudaPSpinor tmp2) {
    MatPCCuda(tmp2, gauge, inEven, kappa, tmp1);
    MatPCDagCuda(outEven, gauge, tmp2, kappa, tmp1);
}
