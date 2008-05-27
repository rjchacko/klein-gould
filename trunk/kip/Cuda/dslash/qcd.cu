
// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// includes, project
#include <cutil.h>

// includes, kernels
#include "qcd.h"
#include "qcd_kernel.cu"


extern "C"
void computeGold(float* reference, float* vdata, float *mdata, const unsigned int len);

void printVector(int i, float *v) {
    printf("%d : {(%f %f) (%f %f) (%f %f)}\n", i, v[0*L], v[1*L], v[2*L], v[3*L], v[4*L], v[5*L]);
}


void
runTest( int argc, char** argv)  {
    CUT_DEVICE_INIT(argc, argv);
    
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    int shared_bytes = SHARED_BYTES;

    printf("Spinors: %d\n", L);
    printf("Global kb: %f\n", L*(GAUGE_BYTES+SPINOR_BYTES)/1024.);
    printf("Shared kb: %fkB\n", shared_bytes/1024.);
    
    // allocate host memory
    float* h_vdata = (float*) malloc(L*SPINOR_BYTES);
    float* h_mdata = (float*) malloc(L*GAUGE_BYTES);
    float* h_odata = (float*) malloc(L*SPINOR_BYTES);
    int *h_index   =   (int*) malloc(8*L*sizeof(int));
    float* reference = (float*) malloc(L*SPINOR_BYTES);
    
    // initialize indices
    for (int d = 0; d < 8; d++) {
    
        for (int i = 0; i < L; i++) {
            h_index[d*L+i] = i;
        }
    
    }
    
    // initialize gauge field
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            for(int i = 0; i < L; i++) {
                h_mdata[L*(2*(3*m + n) + 0) + i] = 1;
                h_mdata[L*(2*(3*m + n) + 1) + i] = -1;
            }
        }
    }
    // initialize spinor field
    for (int s = 0; s < 4; s++) {
        for (int m = 0; m < 3; m++) {
            for(int i = 0; i < L; i++) {
                h_vdata[L*(2*(3*s + m) + 0) + i] = m;
                h_vdata[L*(2*(3*s + m) + 1) + i] = m;
            }
        }
    }
    
    // allocate and load device memory
    float* d_vdata; // input vectors
    float* d_mdata; // matrices
    float* d_odata; // output vectors
    int* d_index;
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_vdata, L*SPINOR_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_mdata, L*GAUGE_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_odata, L*SPINOR_BYTES));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_index, 8*L*sizeof(int)));
    CUDA_SAFE_CALL(cudaMemcpy(d_vdata, h_vdata, L*SPINOR_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_mdata, h_mdata, L*GAUGE_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_mdata, h_mdata, L*GAUGE_BYTES, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_index, h_index, 8*L*sizeof(int), cudaMemcpyHostToDevice));
    
    cudaBindTexture(0 /*offset*/, spinor, d_vdata, L*SPINOR_BYTES); 
    cudaBindTexture(0 /*offset*/, gauge, d_mdata, L*GAUGE_BYTES); 

    // execute kernel
    const int LOOPS = 100;
    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    cutStartTimer(timer);
    for (int i = 0; i < LOOPS; i++) {
        testKernel <<<gridDim, blockDim, shared_bytes>>> ((float4 *)d_odata, d_index);
        CUT_CHECK_ERROR("Kernel execution failed");
        cudaThreadSynchronize();
    }
    cutStopTimer(timer);
    
    float millisecs = cutGetTimerValue(timer)/LOOPS;
    float secs = millisecs / 1000.;
    printf("Elapsed time %fms\n", millisecs);
    printf("GFLOPS = %f\n", 1e-9*1320*L/secs);
    printf("GiB = %f\n", 2*(1*7+4)*3*4*L/(secs*(1<<30)));
    
    // compare to host computation
    CUDA_SAFE_CALL(cudaMemcpy(h_odata, d_odata, L*SPINOR_BYTES, cudaMemcpyDeviceToHost));
    
/*
    cutResetTimer(timer);
    cutStartTimer(timer);
    for (int i = 0; i < 100; i++) {
        computeGold(reference, h_vdata, h_mdata, NUM_MATS);
    }
    cutStopTimer(timer);
    printf( "Host processing time: %f (ms)\n", cutGetTimerValue( timer));
*/

/*
    printf("comparing %f kbytes\n", (L*SPINOR_BYTES)/1024.);
    CUTBoolean res = cutComparefe(reference, h_odata, L*SPINOR_BYTES, 1e-6);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
*/


    for (int i = 0; i < 10; i++) {
//        printVector(i, &reference[i]);
        printVector(i, &h_odata[i]);
    }

    // cleanup memory
    free(h_mdata);
    free(h_vdata);
    free(h_odata);
    free(h_index);
    free(reference);
    CUDA_SAFE_CALL(cudaFree(d_mdata));
    CUDA_SAFE_CALL(cudaFree(d_vdata));
    CUDA_SAFE_CALL(cudaFree(d_odata));
    CUDA_SAFE_CALL(cudaFree(d_index));
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
