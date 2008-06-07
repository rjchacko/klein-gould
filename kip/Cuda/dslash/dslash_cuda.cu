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

float4 *d_gaugeEven, *d_gaugeOdd;
float4 *d_spinorEven, *d_spinorOdd;

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

void packGaugeField(float4 *res, float **gauge) {
    for (int dir = 0; dir < 4; dir++) {
        for (int i = 0; i < Nh; i++) {
            for (int j = 0; j < 5; j++) {
                float a1, a2, a3=0, a4=0;
                a1 = gauge[dir][i*18 + j*4 + 0];
                a2 = gauge[dir][i*18 + j*4 + 1];
                if (j < 4) {
                    a3 = gauge[dir][i*18 + j*4 + 2];
                    a4 = gauge[dir][i*18 + j*4 + 3];
                }
                float4 f4 = {a1, a2, a3, a4};
                res[(dir*5+j)*Nh + i] = f4;
            }
        }
    }
}

void packSpinorField(float4 *res, float *spinor) {
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

void unpackSpinorField(float *res, float4 *spinorPacked) {
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

void sendGaugeField(float **gaugeEven, float **gaugeOdd) {
    float4 *packed1 = (float4*) malloc(Nh*PACKED_GAUGE_BYTES);
    float4 *packed2 = (float4*) malloc(Nh*PACKED_GAUGE_BYTES);
    packGaugeField(packed1, gaugeEven);
    packGaugeField(packed2, gaugeOdd);
    CUDA_SAFE_CALL(cudaMemcpy(d_gaugeEven, packed1, Nh*PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_gaugeOdd, packed2, Nh*PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice));
    free(packed1);
    free(packed2);
}

void sendSpinorFieldEven(float *spinorEven) {
    float4 *packed = (float4*) malloc(Nh*SPINOR_BYTES);
    packSpinorField(packed, spinorEven);
    CUDA_SAFE_CALL(cudaMemcpy(d_spinorEven, packed, Nh*SPINOR_BYTES, cudaMemcpyHostToDevice));
    free(packed);
}

void sendSpinorFieldOdd(float *spinorOdd) {
    float4 *packed = (float4*) malloc(Nh*SPINOR_BYTES);
    packSpinorField(packed, spinorOdd);
    CUDA_SAFE_CALL(cudaMemcpy(d_spinorOdd, packed, Nh*SPINOR_BYTES, cudaMemcpyHostToDevice));
    free(packed);
}

void retrieveSpinorFieldEven(float *res) {
    float4 *packed = (float4*) malloc(Nh*SPINOR_BYTES);
    CUDA_SAFE_CALL(cudaMemcpy(packed, d_spinorEven, Nh*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(res, packed);
    free(packed);
}

void retrieveSpinorFieldOdd(float *res) {
    float4 *packed = (float4*) malloc(Nh*SPINOR_BYTES);
    CUDA_SAFE_CALL(cudaMemcpy(packed, d_spinorOdd, Nh*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(res, packed);
    free(packed);
}

void initializeCuda(int argc, char** argv) {
    CUT_DEVICE_INIT(argc, argv);
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_gaugeEven, Nh*PACKED_GAUGE_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_gaugeOdd, Nh*PACKED_GAUGE_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_spinorEven, Nh*SPINOR_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_spinorOdd, Nh*SPINOR_BYTES));

    printf("Spinors: %d\n", Nh);
    printf("Global kb: %f\n", 2*Nh*(PACKED_GAUGE_BYTES+SPINOR_BYTES)/1024.);
    printf("Shared kb: %fkB\n", SHARED_BYTES/1024.);
}

void releaseCuda() {
    CUDA_SAFE_CALL(cudaFree(d_gaugeEven));
    CUDA_SAFE_CALL(cudaFree(d_gaugeOdd));
    CUDA_SAFE_CALL(cudaFree(d_spinorEven));
    CUDA_SAFE_CALL(cudaFree(d_spinorOdd));
}

void dslashCuda(int oddBit, int daggerBit) {
    if (oddBit) {
        cudaBindTexture(0 /*offset*/, gauge0Tex, d_gaugeOdd, Nh*PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, d_gaugeEven, Nh*PACKED_GAUGE_BYTES); 
    }
    else {
        cudaBindTexture(0 /*offset*/, gauge0Tex, d_gaugeEven, Nh*PACKED_GAUGE_BYTES); 
        cudaBindTexture(0 /*offset*/, gauge1Tex, d_gaugeOdd, Nh*PACKED_GAUGE_BYTES); 
    }
    cudaBindTexture(0 /*offset*/, spinorTex, d_spinorEven, Nh*SPINOR_BYTES); 

    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    
    if (!daggerBit) {
        dslashKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)d_spinorOdd, oddBit);
    }
    else {
        dslashDaggerKernel <<<gridDim, blockDim, SHARED_BYTES>>> ((float4 *)d_spinorOdd, oddBit);
    }
    
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
}


void printSpinorHalfField(float *spinor) {
    printSpinor(&spinor[0*(4*3*2)]);
    printf("...\n");
    printSpinor(&spinor[(Nh-1)*(4*3*2)]);
    printf("\n");    
}

int main(int argc, char **argv) {
    initializeCuda(argc, argv);
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    // construct input fields
    float *gaugeEven[4], *gaugeOdd[4];
    for (int dir = 0; dir < 4; dir++) {
        gaugeEven[dir] = (float*)malloc(Nh*3*3*2*sizeof(float));
        gaugeOdd[dir]  = (float*)malloc(Nh*3*3*2*sizeof(float));
    }
    float *spinorEven = (float*)malloc(Nh*4*3*2*sizeof(float));
    float *spinorOdd  = (float*)malloc(Nh*4*3*2*sizeof(float));
    float *spinorRef  = (float*)malloc(Nh*4*3*2*sizeof(float));

    // copy inputs from host to device
    printf("Randomizing fields\n");
    constructGaugeField(gaugeEven, gaugeOdd);
    constructSpinorField(spinorEven);
    sendGaugeField(gaugeEven, gaugeOdd);
    sendSpinorFieldEven(spinorEven);
    
    int ODD_BIT = 0;
    int DAGGER_BIT = 0;
    
    // execute kernel
    printf("Beginning kernel execution\n");
    cutStartTimer(timer);
    const int LOOPS = 100;
    for (int i = 0; i < LOOPS; i++) {
        dslashCuda(ODD_BIT, DAGGER_BIT);
    }
    cutStopTimer(timer);

    // print timing information
    float millisecs = cutGetTimerValue(timer)/LOOPS;
    float secs = millisecs / 1000.;
    printf("Elapsed time %fms\n", millisecs);
    printf("GFLOPS = %f\n", 1e-9*1320*Nh/secs);
    printf("GiB = %f\n\n", Nh*(8*7+4)*3*2*sizeof(float)/(secs*(1<<30)));

    // compare to dslash reference implementation
    retrieveSpinorFieldOdd(spinorOdd);
    dslashReference(spinorRef, gaugeEven, gaugeOdd, spinorEven, ODD_BIT, DAGGER_BIT);
    printf("Reference:\n");
    printSpinorHalfField(spinorRef);
    printf("\nCUDA:\n");
    printSpinorHalfField(spinorOdd);
    CUTBoolean res = cutComparefe(spinorOdd, spinorRef, Nh*4*3*2, 1e-4);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    
    // release memory
    for (int dir = 0; dir < 4; dir++) {
        free(gaugeEven[dir]);
        free(gaugeOdd[dir]);
    }
    free(spinorEven);
    free(spinorOdd);
    free(spinorRef);
    cutDeleteTimer(timer);
    releaseCuda();
}
