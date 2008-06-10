#include <stdlib.h>
#include <stdio.h>

#include "qcd.h"


#define REDUCE_THREADS 128
#define REDUCE_MAX_BLOCKS 64


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


//
// float sumCuda(float *a, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) sum##suffix
#define REDUCE_TYPES float *a
#define REDUCE_PARAMS a
#define REDUCE_OPERATION(i) a[i]
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_OPERATION

//
// float normCuda(float *a, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) norm##suffix
#define REDUCE_TYPES float *a
#define REDUCE_PARAMS a
#define REDUCE_OPERATION(i) (a[i]*a[i])
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_OPERATION

//
// float reDotProductCuda(float *a, float *b, int n) {}
//
#define REDUCE_FUNC_NAME(suffix) reDotProduct##suffix
#define REDUCE_TYPES float *a, float *b
#define REDUCE_PARAMS a, b
#define REDUCE_OPERATION(i) (a[i]*b[i])
#include "reduce_core.cu"
#undef REDUCE_FUNC_NAME
#undef REDUCE_TYPES
#undef REDUCE_PARAMS
#undef REDUCE_OPERATION




void blasTest() {
    int n = 3*1<<8;
    float *h_data = (float *)malloc(n*sizeof(float));
    float *d_data;
    cudaMalloc((void **)&d_data,  n*sizeof(float));
    
    double acc = 0;
    for (int i = 0; i < n; i++) {
        h_data[i] = i;
        acc += i*i;
    }
    cudaMemcpy(d_data, h_data, n*sizeof(float), cudaMemcpyHostToDevice);
    
    printf("Size: %f MiB\n", (float)n*sizeof(float) / (1 << 20));
    printf("cuda: %f, expected: %f\n", reDotProductCuda(d_data, d_data, n), acc);
    
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
