#include <stdlib.h>
#include <stdio.h>

#define CUDA_DEVICE (0)
#define NUM_THREADS (1<<13)
#define BLOCK_DIM (64)
#define GRID_DIM (NUM_THREADS/BLOCK_DIM)
#define NUM_BYTES (NUM_THREADS*4*sizeof(float))

// Compile and run with the commands:
//   nvcc float4_test.cu
//   ./a.out
//
// Failure occurs on my Tesla C870 card when KERNEL_INVOCATIONS is a large
// number (e.g. 100), and TEST_KERNEL is 1 or 3. Kernels 1 and 3 are those
// which write float4 values to device memory. Failure does not occur on
// my Quadro NVS 290 card for any of the kernels.
//
#define KERNEL_INVOCATIONS (100)
#define TEST_KERNEL (1)


__global__ void testKernel1(float4* g_out, float4* g_in) {
    const int idx = blockDim.x*blockIdx.x + threadIdx.x;
    g_out[idx] = g_in[idx];
}

__global__ void testKernel2(float* g_out, float4* g_in) {
    const int idx = BLOCK_DIM*blockIdx.x + threadIdx.x;
    float4 f4 = g_in[idx];
    g_out[4*idx+0] = f4.x;
    g_out[4*idx+1] = f4.y;
    g_out[4*idx+2] = f4.z;
    g_out[4*idx+3] = f4.w;
}

__global__ void testKernel3(float4* g_out, float* g_in) {
    const int idx = BLOCK_DIM*blockIdx.x + threadIdx.x;
    float x = g_in[4*idx+0];
    float y = g_in[4*idx+1];
    float z = g_in[4*idx+2];
    float w = g_in[4*idx+3];
    g_out[idx] = make_float4(x, y, z, w);
}

__global__ void testKernel4(float* g_out, float* g_in) {
    const int idx = BLOCK_DIM*blockIdx.x + threadIdx.x;
    g_out[NUM_THREADS*0 + idx] = g_in[NUM_THREADS*0 + idx];
    g_out[NUM_THREADS*1 + idx] = g_in[NUM_THREADS*1 + idx];
    g_out[NUM_THREADS*2 + idx] = g_in[NUM_THREADS*2 + idx];
    g_out[NUM_THREADS*3 + idx] = g_in[NUM_THREADS*3 + idx];
}

int main( int argc, char** argv)  {
    cudaSetDevice(CUDA_DEVICE);
    
    float *input = (float *)malloc(NUM_BYTES);
    float *output = (float *)malloc(NUM_BYTES);
    void* d_input;
    void* d_output;
    cudaMalloc(&d_input, NUM_BYTES);
    cudaMalloc(&d_output, NUM_BYTES);
    
    for (int i = 0; i < NUM_THREADS*4; i++) {
        input[i] = i;
    }
    cudaMemcpy(d_input, input, NUM_BYTES, cudaMemcpyHostToDevice);
    
    dim3 gridDim(GRID_DIM, 1, 1);
    dim3 blockDim(BLOCK_DIM, 1, 1);
    for (int i = 0; i < KERNEL_INVOCATIONS; i++) {
        switch (TEST_KERNEL) {
        case 1:
            testKernel1 <<<gridDim, blockDim>>> ((float4 *)d_output, (float4 *)d_input);
            break;
        case 2:
            testKernel2 <<<gridDim, blockDim>>> ((float *)d_output, (float4 *)d_input);
            break;
        case 3:
            testKernel3 <<<gridDim, blockDim>>> ((float4 *)d_output, (float *)d_input);
            break;
        case 4:
            testKernel4 <<<gridDim, blockDim>>> ((float *)d_output, (float *)d_input);
            break;
        }
        cudaThreadSynchronize();
    }
    
    cudaMemcpy(output, d_output, NUM_BYTES, cudaMemcpyDeviceToHost);
    for (int i = 0; i < NUM_THREADS*4; i++) {
        if (output[i] != i) {
            printf("KERNEL=%d FAILED: elem #%d = %f\n", TEST_KERNEL, i, output[i]);
        }
    }
    
    free(input);
    free(output);
    cudaFree(d_input);
    cudaFree(d_output);
}
