
#include <stdio.h>
#include <cutil.h>


#define SHARED_BYTES (0)
#define ARRAY_BYTES (1<<29)
#define ARRAY_FLOATS (ARRAY_BYTES / sizeof(float))
#define BLOCK_DIM 128
#define GRID_DIM 128


texture<float, 1, cudaReadModeElementType> tex1;
texture<float2, 1, cudaReadModeElementType> tex2;
texture<float4, 1, cudaReadModeElementType> tex4;

extern __shared__ float s_data1[];
extern __shared__ float2 s_data2[];
extern __shared__ float4 s_data4[];


__global__ void testKernel1(float *data) {
    int tid = threadIdx.x;
    int i = blockIdx.x*BLOCK_DIM + tid;
    
    while (i < ARRAY_FLOATS) {
        s_data1[tid] = data[i];
        i += GRID_DIM*BLOCK_DIM;
    }
}

__global__ void testKernel1b (float *data) {
    int tid = threadIdx.x;
    int i = blockIdx.x*(2*BLOCK_DIM) + tid;
    
    while (i < ARRAY_FLOATS) {
        s_data1[tid] = data[i+0*BLOCK_DIM];
        s_data1[tid] = data[i+1*BLOCK_DIM];
        i += 2*GRID_DIM*BLOCK_DIM;
    }
}

__global__ void testKernel1T() {
    int tid = threadIdx.x;
    int i = blockIdx.x*BLOCK_DIM + tid;
    
    while (i < ARRAY_FLOATS) {
        s_data1[tid] = tex1Dfetch(tex1, i);
        i += GRID_DIM*BLOCK_DIM;
    }
}

__global__ void testKernel2(float2 *data) {
    int tid = threadIdx.x;
    int i = blockIdx.x*BLOCK_DIM + tid;

    while (i < ARRAY_FLOATS/2) {
        s_data2[tid] = data[i];
        i += GRID_DIM*BLOCK_DIM;
    }
}

__global__ void testKernel2T() {
    int tid = threadIdx.x;
    int i = blockIdx.x*BLOCK_DIM + tid;
    
    while (i < ARRAY_FLOATS/2) {
        s_data2[tid] = tex1Dfetch(tex2, i);
        i += GRID_DIM*BLOCK_DIM;
    }
}

__global__ void testKernel4T() {
    int tid = threadIdx.x;
    int i = blockIdx.x*BLOCK_DIM + tid;
    
    while (i < ARRAY_FLOATS/4) {
        s_data4[tid] = tex1Dfetch(tex4, i);
        i += GRID_DIM*BLOCK_DIM;
    }
}

__global__ void
testKernel4Tb() {
    int tid = threadIdx.x;
    int i = blockIdx.x*(2*BLOCK_DIM) + tid;
    
    while (i < ARRAY_FLOATS/4) {
        s_data4[tid] = tex1Dfetch(tex4, i+0*BLOCK_DIM);
        s_data4[tid] = tex1Dfetch(tex4, i+1*BLOCK_DIM);
        i += 2*GRID_DIM*BLOCK_DIM;
    }
}

void printKernel(int kernel) {
    int flts = kernel / 10;
    int tex = (kernel % 10) / 2;
    int reorder = kernel % 2;
    printf("Read float%d, ", flts);
    printf(tex ? "texture" : "array");
    if (reorder)
        printf(", reordered ");
    printf("\n");
}

double benchmark(int kernel, void *d_idata) {
    int nIters = 100;
    dim3 threads(BLOCK_DIM, 1, 1);
    dim3 grid(GRID_DIM, 1, 1);
    
    cudaEvent_t start, end;
    cudaEventCreate(&start);
    cudaEventCreate(&end);
    cudaEventRecord(start, 0);
    for (int i=0; i < nIters; ++i) {
        switch (kernel) {
        case 10:
            testKernel1<<< grid, threads, SHARED_BYTES >>>((float*)d_idata);
            break;
        case 11:
            testKernel1b<<< grid, threads, SHARED_BYTES >>>((float*)d_idata);
            break;
        case 12:
            testKernel1T<<< grid, threads, SHARED_BYTES >>>();
            break;
            
        case 20:
            testKernel2<<< grid, threads, SHARED_BYTES >>>((float2*)d_idata);
            break;
        case 22:
            testKernel2T<<< grid, threads, SHARED_BYTES >>>();
            break;
            
        case 42:
            testKernel4T<<< grid, threads, SHARED_BYTES >>>();
            break;
        case 43:
            testKernel4Tb<<< grid, threads, SHARED_BYTES >>>();
            break;
        
        default:
            printf("Undefined kernel %d\n", kernel);
            exit(1);
        }
    }
    cudaEventRecord(end, 0);
    cudaEventSynchronize(end);
    float runTime;
    cudaEventElapsedTime(&runTime, start, end);
    double secs = runTime / 1000.;
    return secs / nIters;
}


int main(int argc, char** argv) {
    double MiB = (double)ARRAY_BYTES / (1<<20);
    double GiB = (double)ARRAY_BYTES / (1<<30);

    CUT_DEVICE_INIT(argc, argv);

    void *d_idata;
    cudaMalloc(&d_idata, ARRAY_BYTES);
    cudaBindTexture(0 /*offset*/, tex1, d_idata, ARRAY_BYTES); 
    cudaBindTexture(0 /*offset*/, tex2, d_idata, ARRAY_BYTES); 
    cudaBindTexture(0 /*offset*/, tex4, d_idata, ARRAY_BYTES); 

    printf("Block dim = %d, Grid dim = %d\n", BLOCK_DIM, GRID_DIM);
    printf("Array size = %f MiB\n\n", MiB);
    
    int kernels[] = {10, 11, 12, 20, 22, 42, 43};
    for (int i = 0; i < 7; i++) {
        printKernel(kernels[i]);
        double secs = benchmark(kernels[i], d_idata);
        printf("Average time: %f s\n", secs);
        printf("Bandwidth:    %f GiB/s\n\n", GiB / secs);
    }
}
