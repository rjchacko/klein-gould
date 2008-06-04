
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cutil.h>

// includes, kernels
#include "qcd.h"


texture<float4, 1, cudaReadModeElementType> gaugeTex;
texture<float4, 1, cudaReadModeElementType> spinorTex;


__global__ void
testKernel(float4* g_out) {
    #include "qcd_core.cu"
}




void packGaugeField(float4 *res, float **gauge) {
    for (int dir = 0; dir < 4; dir++) {
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < 5; j++) {
                float a1, a2, a3=0, a4=0;
                a1 = gauge[dir][i*18 + j*4 + 0];
                a2 = gauge[dir][i*18 + j*4 + 1];
                if (j < 4) {
                    a3 = gauge[dir][i*18 + j*4 + 2];
                    a4 = gauge[dir][i*18 + j*4 + 3];
                }
                float4 f4 = {a1, a2, a3, a4};
                res[(dir*5+j)*L + i] = f4;
            }
        }
    }
}

void packSpinorField(float4 *res, float *spinor) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < 6; j++) {
            float a1 = spinor[i*(6*4) + j*(4) + 0];
            float a2 = spinor[i*(6*4) + j*(4) + 1];
            float a3 = spinor[i*(6*4) + j*(4) + 2];
            float a4 = spinor[i*(6*4) + j*(4) + 3];
            float4 f4 = {a1, a2, a3, a4};
            res[j*L + i] = f4;
        }
    }
}

void unpackSpinorField(float *res, float4 *spinorPacked) {
    for (int i = 0; i < L; i++) {
        if (0) {
            for (int j = 0; j < 6; j++) {
                float4 f4 = spinorPacked[j*L + i];
                res[i*(6*4) + j*(4) + 0] = f4.x;
                res[i*(6*4) + j*(4) + 1] = f4.y;
                res[i*(6*4) + j*(4) + 2] = f4.z;
                res[i*(6*4) + j*(4) + 3] = f4.w;
            }
        }
        else {
            for (int j = 0; j < 24; j++) {
                res[i*24 + j] = ((float *)spinorPacked)[j*L + i];
            }
        }
    }
}


void runTest( int argc, char** argv)  {
    CUT_DEVICE_INIT(argc, argv);
    
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    int shared_bytes = SHARED_BYTES;
    printf("Spinors: %d\n", L);
    printf("Global kb: %f\n", L*(PACKED_GAUGE_BYTES+SPINOR_BYTES)/1024.);
    printf("Shared kb: %fkB\n", shared_bytes/1024.);
    
    // construct input fields
    float *gaugeIn[4];
    for (int dir = 0; dir < 4; dir++)
        gaugeIn[dir] = (float *)malloc(L*3*3*2*sizeof(float));
    float* spinorIn = (float *)malloc(L*SPINOR_BYTES);
    float* spinorOut = (float*)malloc(L*SPINOR_BYTES);
    float* spinorRef = (float*)malloc(L*SPINOR_BYTES);

    // allocate host memory
    float4* packedGaugeIn = (float4*) malloc(L*PACKED_GAUGE_BYTES);
    float4* packedSpinorIn = (float4*) malloc(L*SPINOR_BYTES);
    float4* packedSpinorOut = (float4*) malloc(L*SPINOR_BYTES);
    
    // allocate device memory
    float4* d_gaugeIn;
    float4* d_spinorIn;
    float4* d_spinorOut;
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_gaugeIn, L*PACKED_GAUGE_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_spinorIn, L*SPINOR_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_spinorOut, L*SPINOR_BYTES));
    
    // copy inputs from host to device
    printf("Randomizing fields\n");
    constructGaugeField(gaugeIn);
    constructSpinorField(spinorIn);
    packGaugeField(packedGaugeIn, gaugeIn);
    packSpinorField(packedSpinorIn, spinorIn);

    CUDA_SAFE_CALL(cudaMemcpy(d_gaugeIn, packedGaugeIn, L*PACKED_GAUGE_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_spinorIn, packedSpinorIn, L*SPINOR_BYTES, cudaMemcpyHostToDevice));
    
    // bind textures
    cudaBindTexture(0 /*offset*/, gaugeTex, d_gaugeIn, L*PACKED_GAUGE_BYTES); 
    cudaBindTexture(0 /*offset*/, spinorTex, d_spinorIn, L*SPINOR_BYTES); 
    
    // execute kernel
    const int LOOPS = 100;
    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    printf("Beginning kernel execution\n");
    cutStartTimer(timer);
    for (int i = 0; i < LOOPS; i++) {
        testKernel <<<gridDim, blockDim, shared_bytes>>> ((float4 *)d_spinorOut);
        CUT_CHECK_ERROR("Kernel execution failed");
        cudaThreadSynchronize();
    }
    cutStopTimer(timer);
    
    // print timing information
    float millisecs = cutGetTimerValue(timer)/LOOPS;
    float secs = millisecs / 1000.;
    printf("Elapsed time %fms\n", millisecs);
    printf("GFLOPS = %f\n", 1e-9*1320*L/secs);
    printf("GiB = %f\n\n", L*(8*7+4)*3*2*sizeof(float)/(secs*(1<<30)));
    
    // retrieve output spinor from device
    CUDA_SAFE_CALL(cudaMemcpy(packedSpinorOut, d_spinorOut, L*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(spinorOut, packedSpinorOut);
    
    // compare to dslash reference implementation
    computeGold(spinorRef, gaugeIn, spinorIn);
    printf("Reference:\n");
    printSpinorField(spinorRef);
    printf("\nCUDA:\n");
    printSpinorField(spinorOut);
    CUTBoolean res = cutComparefe(spinorOut, spinorRef, L*SPINOR_SIZE, 1e-4);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    
    // release memory
    for (int dir = 0; dir < 4; dir++)
        free(gaugeIn[dir]);
    free(spinorIn);
    free(spinorOut);
    free(spinorRef);
    free(packedGaugeIn);
    free(packedSpinorIn);
    free(packedSpinorOut);
    CUDA_SAFE_CALL(cudaFree(d_gaugeIn));
    CUDA_SAFE_CALL(cudaFree(d_spinorIn));
    CUDA_SAFE_CALL(cudaFree(d_spinorOut));
    cutDeleteTimer(timer);
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) {
    runTest( argc, argv);
    // CUT_EXIT(argc, argv);
}
