
#include <stdio.h>
#include <cutil.h>


#define SHARED_BYTES (0)

#define ARRAY_BYTES (1<<29)
#define ARRAY_FLOATS (ARRAY_BYTES / sizeof(float))
#define BLOCK_DIM 128
#define GRID_DIM 128

#define THREADS (BLOCK_DIM*GRID_DIM)
#define LOOPS (ARRAY_FLOATS/(8*THREADS))


texture<float2, 1, cudaReadModeElementType> tex2;
texture<float4, 1, cudaReadModeElementType> tex4;

extern __shared__ float s_data1[];
extern __shared__ float2 s_data2[];
extern __shared__ float4 s_data4[];


__global__ void testKernel1(float *data) {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    float acc = 0;
    for (int i = 0; i < ARRAY_FLOATS/(2*THREADS); i++) {
        acc += data[sp_idx]; sp_idx += THREADS;
        acc += data[sp_idx]; sp_idx += THREADS;
    }
    s_data1[tid] = acc;
}


__global__ void testKernel1b (float *data) {
    unsigned int tid = threadIdx.x;
    unsigned int bid = blockIdx.x;
    unsigned int gridSize = GRID_DIM*(2*BLOCK_DIM);
    unsigned int i = bid*(2*BLOCK_DIM) + tid;
    
    float acc = 0;
    while (i < ARRAY_FLOATS) {
        acc += data[i] + data[BLOCK_DIM+i];
        i += gridSize;
    }
    
    s_data1[tid] = acc;
}


__global__ void
testKernel2(float2 *data) {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < LOOPS; i++) {
        s_data2[tid] = data[sp_idx]; sp_idx += THREADS;
        s_data2[tid] = data[sp_idx]; sp_idx += THREADS;
        s_data2[tid] = data[sp_idx]; sp_idx += THREADS;
        s_data2[tid] = data[sp_idx]; sp_idx += THREADS;
    }
}

__global__ void
testKernel2T() {
    const int bid = blockIdx.x;
    int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < LOOPS; i++) {
        s_data2[tid] = tex1Dfetch(tex2, sp_idx+=THREADS);
        s_data2[tid] = tex1Dfetch(tex2, sp_idx+=THREADS);
        s_data2[tid] = tex1Dfetch(tex2, sp_idx+=THREADS);
        s_data2[tid] = tex1Dfetch(tex2, sp_idx+=THREADS);
    }
}

__global__ void
testKernel4T() {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < LOOPS/2; i++) {
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=THREADS);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=THREADS);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=THREADS);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=THREADS);
    }
}


int main(int argc, char** argv) {
    int nIters = 100;
    double MiB = (double)ARRAY_BYTES / (1<<20);
    double GiB = (double)ARRAY_BYTES / (1<<30);

    CUT_DEVICE_INIT(argc, argv);
    // CUDA_SAFE_CALL( cudaSetDevice(1) );

    void *d_idata;
    CUDA_SAFE_CALL( cudaMalloc(&d_idata, ARRAY_BYTES) );
    CUDA_SAFE_CALL( cudaBindTexture(0 /*offset*/, tex2, d_idata, ARRAY_BYTES) ); 
    CUDA_SAFE_CALL( cudaBindTexture(0 /*offset*/, tex4, d_idata, ARRAY_BYTES) ); 

    printf("Global MiB: %f\n", MiB);
    printf("Shared kb: %fkB\n", SHARED_BYTES/1024.);
    
    dim3 threads(BLOCK_DIM, 1, 1);
    dim3 grid(GRID_DIM, 1, 1);
    
    cudaEvent_t start, end;
    CUDA_SAFE_CALL( cudaEventCreate(&start) );
    CUDA_SAFE_CALL( cudaEventCreate(&end) );
    CUDA_SAFE_CALL( cudaEventRecord(start, 0) );
    for (int i=0; i < nIters; ++i) {
        switch (-2) {
        case -2:
            if (i == 0) printf("Reduction float\n");
            testKernel1b<<< grid, threads, SHARED_BYTES >>>((float*)d_idata);
            break;
        case -1:
            if (i == 0) printf("Straight float\n");
            testKernel1<<< grid, threads, SHARED_BYTES >>>((float*)d_idata);
            break;
        case 0:
            if (i == 0) printf("Straight float2\n");
            testKernel2<<< grid, threads, SHARED_BYTES >>>((float2*)d_idata);
            break;
        case 1:
            if (i == 0) printf("Texture float2\n");
            testKernel2T<<< grid, threads, SHARED_BYTES >>>();
            break;
        case 2:
            if (i == 0) printf("Texture float4\n");
            testKernel4T<<< grid, threads, SHARED_BYTES >>>();
            break;
        }
    }
    CUDA_SAFE_CALL( cudaEventRecord(end, 0) );
    CUDA_SAFE_CALL( cudaEventSynchronize(end) );

    float runTime;
    CUDA_SAFE_CALL( cudaEventElapsedTime(&runTime, start, end) );
    double secs = runTime / 1000.;
    
    printf("Average time: %f ms\n", runTime);
    printf("Bandwidth:    %f GiB/s\n\n", GiB*nIters / secs);
}
