#define RP REDUCE_PRECISION


__global__ void REDUCE_FUNC_NAME(Kernel) (REDUCE_TYPES, RP *g_odata, unsigned int n) {
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(2*REDUCE_THREADS) + threadIdx.x;
    unsigned int gridSize = gridDim.x*(2*REDUCE_THREADS);
    
    RP acc = 0;
    while (i < n) {
        REDUCE_AUXILIARY(i);
        REDUCE_AUXILIARY(i+REDUCE_THREADS);        
        acc += (RP)REDUCE_OPERATION(i) + (RP)REDUCE_OPERATION(i+REDUCE_THREADS);
        i += gridSize;
    }

    extern __shared__ RP s[];
    s[tid] = acc;
    __syncthreads();
    
    // do reduction in shared mem
    if (REDUCE_THREADS >= 256) { if (tid < 128) { s[tid] += s[128+tid]; } __syncthreads(); }
    if (REDUCE_THREADS >= 128) { if (tid <  64) { s[tid] += s[ 64+tid]; } __syncthreads(); }
    
    if (tid < 32) {
        s[tid] += s[32+tid];
        s[tid] += s[16+tid];
        s[tid] += s[ 8+tid];
        s[tid] += s[ 4+tid];
        s[tid] += s[ 2+tid];
        s[tid] += s[ 1+tid];
    }
    
    if (tid == 0) g_odata[blockIdx.x] = s[0];
}


float REDUCE_FUNC_NAME(Cuda) (REDUCE_TYPES, int n) {
    if (n % (REDUCE_THREADS) != 0) {
        printf("ERROR reduceCuda(): length must be a multiple of %d\n", (REDUCE_THREADS));
        return 0.;
    }
    
    // allocate arrays on device and host to store one float for each block
    int blocks = min(REDUCE_MAX_BLOCKS, n / (2*REDUCE_THREADS));
    RP h_odata[REDUCE_MAX_BLOCKS];
    RP *d_odata;
    cudaMalloc((void**) &d_odata, blocks*sizeof(RP));
    
    // partial reduction; each block generates one number
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = REDUCE_THREADS * sizeof(RP);
    REDUCE_FUNC_NAME(Kernel)<<< dimGrid, dimBlock, smemSize >>>(REDUCE_PARAMS, d_odata, n);
    
    // copy result from device to host, and perform final reduction on CPU
    cudaMemcpy(h_odata, d_odata, blocks*sizeof(RP), cudaMemcpyDeviceToHost);
    RP gpu_result = 0;
    for (int i = 0; i < blocks; i++)  {
        gpu_result += h_odata[i];
    }
    
    cudaFree(d_odata);    
    return (float)gpu_result;
}
