
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

#define BLOCK_DIM (32*dim)
#define GRID_DIM (1<<15)
#define NUM_MATS (GRID_DIM*BLOCK_DIM/dim)

extern "C"
void computeGold(float* reference, float* vdata, float *mdata, const unsigned int len);

void printVector(int i, float *v) {
    printf("%d : {(%f %f) (%f %f) (%f %f)}\n", i, v[0], v[1], v[2], v[3], v[4], v[5]);
}


////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
void
runTest( int argc, char** argv)  {
    CUT_DEVICE_INIT();
    
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    unsigned int vec_bytes = vecsize*NUM_MATS*sizeof(float);
    unsigned int mat_bytes = matsize*NUM_MATS*sizeof(float);
    
    printf("Matrix bytes %d\n", mat_bytes+2*vec_bytes);
    
    // allocate host memory
    float* h_vdata = (float*) malloc(vec_bytes);
    float* h_mdata = (float*) malloc(mat_bytes);
    float* h_odata = (float*) malloc(vec_bytes);
    float* reference = (float*) malloc(vec_bytes);

    // initialize matrices and vectors
    for(int i = 0; i < NUM_MATS; i++) {
        for (int m = 0; m < dim; m++) {
            for (int n = 0; n < dim; n++) {
                h_mdata[i*matsize + m*vecsize + 2*n + 0] = i%3 - 3.2;
                h_mdata[i*matsize + m*vecsize + 2*n + 1] = i%4;
            }
            h_vdata[i*vecsize + 2*m + 0] = i%2 + 0.2;
            h_vdata[i*vecsize + 2*m + 1] = -0.2;
        }
    }
    
    // allocate and load device memory
    float* d_vdata; // input vectors
    float* d_mdata; // matrices
    float* d_odata; // output vectors
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_vdata, vec_bytes));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_mdata, mat_bytes));
    CUDA_SAFE_CALL(cudaMalloc((void**) &d_odata, vec_bytes));
    CUDA_SAFE_CALL(cudaMemcpy(d_vdata, h_vdata, vec_bytes, cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(d_mdata, h_mdata, mat_bytes, cudaMemcpyHostToDevice));
    
    // execute kernel
    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    int shared_bytes = sizeof(float) * (BLOCK_DIM/dim) * (2*vecsize + matsize);
    cutStartTimer(timer);
    for (int i = 0; i < 1000; i++) {
        testKernel <<<gridDim, blockDim, shared_bytes>>> (d_vdata, d_mdata, d_odata);
        CUT_CHECK_ERROR("Kernel execution failed");
    }
    cutStopTimer(timer);
    printf( "Cuda processing time: %f (ms)\n", cutGetTimerValue(timer));
    
    // compare to host computation
    CUDA_SAFE_CALL(cudaMemcpy(h_odata, d_odata, vec_bytes, cudaMemcpyDeviceToHost));
    
    cutResetTimer(timer);
    cutStartTimer(timer);
    for (int i = 0; i < 100; i++) {
        computeGold(reference, h_vdata, h_mdata, NUM_MATS);
    }
    cutStopTimer(timer);
    printf( "Host processing time: %f (ms)\n", cutGetTimerValue( timer));

//    for (int i = 0; i < 10; i++) {
//        printVector(i, &reference[i*vecsize]);
//        printVector(i, &h_odata[i*vecsize]);
//    }
    
    if(cutCheckCmdLineFlag(argc, (const char**) argv, "write")) {
        CUT_SAFE_CALL(cutWriteFilef("./data/regression.dat", h_odata, dim*NUM_MATS, 0.0));
    }
    else {
        printf("comparing %d bytes\n", vecsize*NUM_MATS*4);
        CUTBoolean res = cutComparefe(reference, h_odata, vecsize*NUM_MATS, 1e-6);
        printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    }

    
    // cleanup memory
    free(h_mdata);
    free(h_vdata);
    free(h_odata);
    free(reference);
    CUDA_SAFE_CALL(cudaFree(d_mdata));
    CUDA_SAFE_CALL(cudaFree(d_vdata));
    CUDA_SAFE_CALL(cudaFree(d_odata));
//    cutDeleteTimer(timer);
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) {
    runTest( argc, argv);
    // CUT_EXIT(argc, argv);
}
