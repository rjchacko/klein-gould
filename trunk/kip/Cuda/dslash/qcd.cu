
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


void runTest( int argc, char** argv)  {
    CUT_DEVICE_INIT(argc, argv);
    
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    int shared_bytes = SHARED_BYTES;
    printf("Spinors: %d\n", L);
    printf("Global kb: %f\n", L*(PACKED_GAUGE_BYTES+SPINOR_BYTES)/1024.);
    printf("Shared kb: %fkB\n", shared_bytes/1024.);
    
    // construct input fields
    float* gaugeIn = (float *)malloc(L*GAUGE_BYTES);
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
    cutStartTimer(timer);
    for (int i = 0; i < LOOPS; i++) {
        testKernel <<<gridDim, blockDim, shared_bytes>>> ((float4 *)d_spinorOut);
        CUT_CHECK_ERROR("Kernel execution failed");
        cudaThreadSynchronize();
    }
    cutStopTimer(timer);
    
    float millisecs = cutGetTimerValue(timer)/LOOPS;
    float secs = millisecs / 1000.;
    printf("Elapsed time %fms\n", millisecs);
    printf("GFLOPS = %f\n", 1e-9*1320*L/secs);
    printf("GiB = %f\n", 2*(8*7+4)*3*4*L/(secs*(1<<30)));
    
    // compare to host computation
    CUDA_SAFE_CALL(cudaMemcpy(packedSpinorOut, d_spinorOut, L*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    unpackSpinorField(spinorOut, packedSpinorOut);
    computeGold(spinorRef, gaugeIn, spinorIn);
    CUTBoolean res = cutComparefe(spinorOut, spinorRef, L*SPINOR_SIZE, 1e-4);
    
    printSpinorField(spinorRef);
    printf("\n");
    printSpinorField(spinorOut);

    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    
    
    testSpinorField(spinorOut);
    
    
    // cleanup memory
    free(gaugeIn);
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
