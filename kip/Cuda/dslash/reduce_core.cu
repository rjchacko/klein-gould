

template <unsigned int blockSize>
__global__ void REDUCE_FUNC_NAME(Kernel) (REDUCE_TYPES, float *g_odata, unsigned int n) {
    extern __shared__ float sdata[];

    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*(blockSize*2) + threadIdx.x;
    unsigned int gridSize = blockSize*2*gridDim.x;
    sdata[tid] = 0;

    // we reduce multiple elements per thread.  The number is determined by the 
    // number of active thread blocks (via gridSize).  More blocks will result
    // in a larger gridSize and therefore fewer elements per thread
    while (i < n)
    {
        sdata[tid] += REDUCE_OPERATION(i) + REDUCE_OPERATION(i+blockSize);  
        i += gridSize;
    } 
    __syncthreads();

    // do reduction in shared mem
    if (blockSize >= 512) { if (tid < 256) { sdata[tid] += sdata[tid + 256]; } __syncthreads(); }
    if (blockSize >= 256) { if (tid < 128) { sdata[tid] += sdata[tid + 128]; } __syncthreads(); }
    if (blockSize >= 128) { if (tid <  64) { sdata[tid] += sdata[tid +  64]; } __syncthreads(); }
    
    if (tid < 32) {
        if (blockSize >=  64) { sdata[tid] += sdata[tid + 32]; }
        if (blockSize >=  32) { sdata[tid] += sdata[tid + 16]; }
        if (blockSize >=  16) { sdata[tid] += sdata[tid +  8]; }
        if (blockSize >=   8) { sdata[tid] += sdata[tid +  4]; }
        if (blockSize >=   4) { sdata[tid] += sdata[tid +  2]; }
        if (blockSize >=   2) { sdata[tid] += sdata[tid +  1]; }
    }
    
    // write result for this block to global mem 
    if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


float REDUCE_FUNC_NAME(Cuda) (REDUCE_TYPES, int n) {
    if (n % (2*REDUCE_THREADS) != 0) {
        printf("ERROR reduceCuda(): length must be a multiple of %d\n", (2*REDUCE_THREADS));
        return 0.;
    }
    
    // allocate arrays on device and host to store one float for each block
    int blocks = min(REDUCE_MAX_BLOCKS, n / (2*REDUCE_THREADS));
    float h_odata[REDUCE_MAX_BLOCKS];
    float *d_odata;
    CUDA_SAFE_CALL( cudaMalloc((void**) &d_odata, blocks*sizeof(float)) );
    
    // partial reduction; each block generates one number
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = REDUCE_THREADS * sizeof(float);
    REDUCE_FUNC_NAME(Kernel)<REDUCE_THREADS><<< dimGrid, dimBlock, smemSize >>>(REDUCE_PARAMS, d_odata, n);
    CUT_CHECK_ERROR("Kernel execution failed");

    // copy result from device to host, and perform final reduction on CPU
    CUDA_SAFE_CALL( cudaMemcpy( h_odata, d_odata, blocks*sizeof(float), cudaMemcpyDeviceToHost) );
    float gpu_result = 0;
    for (int i = 0; i < blocks; i++) 
        gpu_result += h_odata[i];
    
    CUDA_SAFE_CALL( cudaFree(d_odata) );    
    return gpu_result;
}
