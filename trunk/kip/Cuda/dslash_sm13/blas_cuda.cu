#include <stdlib.h>
#include <stdio.h>

#include "qcd.h"


//#define REDUCE_PRECISION float
#define REDUCE_PRECISION double
#define REDUCE_THREADS 128
#define REDUCE_MAX_BLOCKS 128


void zeroCuda(float* dst, int len) {
    // cuda's floating point format, IEEE-754, represents the floating point
    // zero as 4 zero bytes
    cudaMemset(dst, 0, len*sizeof(float));
}

void copyCuda(float* dst, float *src, int len) {
    cudaMemcpy(dst, src, len*sizeof(float), cudaMemcpyDeviceToDevice);
}


__global__ void axpbyKernel(float a, float *x, float b, float *y, int len) {
    unsigned int i = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;
    while (i < len) {
        y[i] = a*x[i] + b*y[i];
        i += gridSize;
    } 
}

// performs the operation y[i] = a*x[i] + b*y[i]
void axpbyCuda(float a, float *x, float b, float *y, int len) {
    int blocks = min(REDUCE_MAX_BLOCKS, max(len/REDUCE_THREADS, 1));
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    axpbyKernel<<<dimGrid, dimBlock>>>(a, x, b, y, len);
}


// performs the operation y[i] = a*x[i] + y[i]
void axpyCuda(float a, float *x, float *y, int len) {
    axpbyCuda(a, x, 1.0, y, len);
}

// performs the operation y[i] = x[i] + a*y[i]
void xpayCuda(float *x, float a, float *y, int len) {
    axpbyCuda(1.0, x, a, y, len);
}

// performs the operation y[i] -= x[i] (minus x plus y)
void mxpyCuda(float *x, float *y, int len) {
    axpbyCuda(-1.0, x, 1.0, y, len);
}

__global__ void axpyZpbxKernel(float a, float *x, float *y, float *z, float b, int len) {
    unsigned int i = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;
    while (i < len) {
        float x_i = x[i];
        y[i] = a*x_i + y[i];
        x[i] = z[i] + b*x_i;
        i += gridSize;
    }
}

// performs the operations: {y[i] = a x[i] + y[i]; x[i] = z[i] + b x[i]}
void axpyZpbxCuda(float a, float *x, float *y, float *z, float b, int len) {
    int blocks = min(REDUCE_MAX_BLOCKS, max(len/REDUCE_THREADS, 1));
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    axpyZpbxKernel<<<dimGrid, dimBlock>>>(a, x, y, z, b, len);
}

// performs the operation y[i] = a*x[i] + y[i], and returns norm(y)
// float axpyNormCuda(float a, float *x, float *y, int len);


// Computes c = a + b in "double single" precision.
__device__ void dsadd(float &c0, float &c1, const float a0, const float a1, const float b0, const float b1) {
    // Compute dsa + dsb using Knuth's trick.
    float t1 = a0 + b0;
    float e = t1 - a0;
    float t2 = ((b0 - e) + (a0 - (t1 - e))) + a1 + b1;
    // The result is t1 + t2, after normalization.
    c0 = e = t1 + t2;
    c1 = t2 - (e - t1);
}

//
// float sumCuda(float *a, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) sum##suffix
#define REDUCE_TYPES float *a
#define REDUCE_PARAMS a
#define REDUCE_AUXILIARY(i)
#define REDUCE_OPERATION(i) a[i]
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_AUXILIARY
#undef REDUCE_OPERATION

//
// float normCuda(float *a, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) norm##suffix
#define REDUCE_TYPES float *a
#define REDUCE_PARAMS a
#define REDUCE_AUXILIARY(i)
#define REDUCE_OPERATION(i) (a[i]*a[i])
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_AUXILIARY
#undef REDUCE_OPERATION

//
// float reDotProductCuda(float *a, float *b, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) reDotProduct##suffix
#define REDUCE_TYPES float *a, float *b
#define REDUCE_PARAMS a, b
#define REDUCE_AUXILIARY(i)
#define REDUCE_OPERATION(i) (a[i]*b[i])
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_AUXILIARY
#undef REDUCE_OPERATION


//
// float axpyNormCuda(float a, float *x, float *y, n){}
//
// First performs the operation y[i] = a*x[i] + y[i]
// Second returns the norm of y
//

#define REDUCE_FUNC_NAME(suffix) axpyNorm##suffix
#define REDUCE_TYPES float a, float *x, float *y
#define REDUCE_PARAMS a, x, y
#define REDUCE_AUXILIARY(i) y[i] = a*x[i] + y[i]
#define REDUCE_OPERATION(i) (y[i]*y[i])
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_AUXILIARY
#undef REDUCE_OPERATION



double cpuDouble(float *data, int size) {
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += data[i];
    return sum;
}

void blasTest() {
    int n = 1<<26;
    double mib = (double)n*sizeof(float) / (1 << 20);
    
    float *h_data = (float *)malloc(n*sizeof(float));
    float *d_data;
    cudaMalloc((void **)&d_data,  n*sizeof(float));
    
    for (int i = 0; i < n; i++) {
        h_data[i] = rand()/(float)RAND_MAX - 0.5; // n-1.0-i;
    }
    
    cudaMemcpy(d_data, h_data, n*sizeof(float), cudaMemcpyHostToDevice);
    
    cudaThreadSynchronize();
    stopwatchStart();
    int LOOPS = 20;
    for (int i = 0; i < LOOPS; i++) {
        sumCuda(d_data, n);
    }
    cudaThreadSynchronize();
    float secs = stopwatchReadSeconds();
    
    printf("%f GiB/s\n", (mib/1024) * LOOPS / secs);
    printf("Device: %f MiB\n", mib);
    printf("Shared: %f KiB\n", (float)REDUCE_THREADS*sizeof(float) / (1 << 10));

    float correctDouble = cpuDouble(h_data, n);
    printf("CPU: %f\n", correctDouble);
    printf("CUDA: %f\n", sumCuda(d_data, n));
    printf("Error: %f\n", fabs(correctDouble-sumCuda(d_data, n)));
    
    cudaFree(d_data) ;
    free(h_data);
}

void axpbyTest() {
    int n = 3 * 1 << 20;
    float *h_x = (float *)malloc(n*sizeof(float));
    float *h_y = (float *)malloc(n*sizeof(float));
    float *h_res = (float *)malloc(n*sizeof(float));
    
    float *d_x, *d_y;
    cudaMalloc((void **)&d_x,  n*sizeof(float));
    cudaMalloc((void **)&d_y,  n*sizeof(float));
    
    for (int i = 0; i < n; i++) {
        h_x[i] = 1;
        h_y[i] = 2;
    }
    
    cudaMemcpy(d_x, h_x, n*sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(d_y, h_y, n*sizeof(float), cudaMemcpyHostToDevice);
    
    axpbyCuda(4, d_x, 3, d_y, n/2);
    
    cudaMemcpy( h_res, d_y, n*sizeof(float), cudaMemcpyDeviceToHost);

    for (int i = 0; i < n; i++) {
        float expect = (i < n/2) ? 4*h_x[i] + 3*h_y[i] : h_y[i];
        if (h_res[i] != expect)
            printf("FAILED %d : %f != %f\n", i, h_res[i], h_y[i]);
    }
    
    cudaFree(d_y);
    cudaFree(d_x);
    free(h_x);
    free(h_y);
}
