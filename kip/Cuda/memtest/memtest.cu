
#include <stdio.h>

#include <cutil.h>


#define L (1<<14)
#define BLOCK_DIM 64
#define GRID_DIM (L/BLOCK_DIM)
#define ROWS 200
#define ARRAY_BYTES (L*ROWS*8*sizeof(float))



texture<float2, 1, cudaReadModeElementType> tex2;
texture<float4, 1, cudaReadModeElementType> tex4;

extern __shared__ float s_data1[];
extern __shared__ float2 s_data[];
extern __shared__ float4 s_data4[];

__global__ void
testKernel1(float *data) {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < ROWS; i++) {
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
        s_data1[tid] = data[sp_idx]; sp_idx += L;
    }
}

__global__ void
testKernel2(float2 *data) {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < ROWS; i++) {
        s_data[tid] = data[sp_idx]; sp_idx += L;
        s_data[tid] = data[sp_idx]; sp_idx += L;
        s_data[tid] = data[sp_idx]; sp_idx += L;
        s_data[tid] = data[sp_idx]; sp_idx += L;
    }
}

__global__ void
testKernel2T() {
    const int bid = blockIdx.x;
    int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < ROWS; i++) {
        s_data[tid] = tex1Dfetch(tex2, sp_idx+=L);
        s_data[tid] = tex1Dfetch(tex2, sp_idx+=L);
        s_data[tid] = tex1Dfetch(tex2, sp_idx+=L);
        s_data[tid] = tex1Dfetch(tex2, sp_idx+=L);
    }
}

__global__ void
testKernel4T() {
    const int bid = blockIdx.x;
    const int tid = threadIdx.x;
    const int sid = BLOCK_DIM*bid + tid;
    int sp_idx = sid;
    
    for (int i = 0; i < ROWS/2; i++) {
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=L);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=L);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=L);
        s_data4[tid] = tex1Dfetch(tex4, sp_idx+=L);
    }
}


int main() {
    int nIters = 100;
    int sharedBytes = 4500;
    
    CUDA_SAFE_CALL( cudaSetDevice(1) );

    float2 *d_idata;
    CUDA_SAFE_CALL( cudaMalloc((void**)&d_idata, ARRAY_BYTES) );
    CUDA_SAFE_CALL( cudaBindTexture(0 /*offset*/, tex2, d_idata, ARRAY_BYTES) ); 
    CUDA_SAFE_CALL( cudaBindTexture(0 /*offset*/, tex4, (float4*)d_idata, ARRAY_BYTES) ); 

    printf("Global kb: %f\n", ARRAY_BYTES/1024.);
    printf("Shared kb: %fkB\n", sharedBytes/1024.);


    dim3 threads(BLOCK_DIM, 1, 1);
    dim3 grid(GRID_DIM, 1, 1);
    
    cudaEvent_t start, end;
    CUDA_SAFE_CALL( cudaEventCreate(&start) );
    CUDA_SAFE_CALL( cudaEventCreate(&end) );
    CUDA_SAFE_CALL( cudaEventRecord(start, 0) );
    for (int i=0; i < nIters; ++i) {
        switch (-1) {
        case -1:
            if (i == 0) printf("Straight float\n");
            testKernel1<<< grid, threads, sharedBytes >>>((float*)d_idata);
            break;
        case 0:
            if (i == 0) printf("Straight float2\n");
            testKernel2<<< grid, threads, sharedBytes >>>(d_idata);
            break;
        case 1:
            if (i == 0) printf("Texture float2\n");
            testKernel2T<<< grid, threads, sharedBytes >>>();
            break;
        case 2:
            if (i == 0) printf("Texture float4\n");
            testKernel4T<<< grid, threads, sharedBytes >>>();
            break;
        }
    }
    CUDA_SAFE_CALL( cudaEventRecord(end, 0) );
    CUDA_SAFE_CALL( cudaEventSynchronize(end) );

    float runTime;
    CUDA_SAFE_CALL( cudaEventElapsedTime(&runTime, start, end) );
    runTime /= float(nIters);
    printf("Average time: %f ms\n", runTime);
    printf("Bandwidth:    %f GiB/s\n\n", ARRAY_BYTES / (runTime * 1.0e-3 * 1024*1024*1024));
}
