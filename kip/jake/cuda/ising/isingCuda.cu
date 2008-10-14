#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "ising.h"


// ----------------------------------------------------------------------------
// Cuda implementation
//

#define GRID_DIM 128
#define BLOCK_DIM 128

#define CUDA_INCLUDE
#include "bits.cpp" // include CUDA device functions
#include "rand48.cu" // random numbers


// ----------------------------------------------------------------------------
// Cuda magnetization calculation

#define REDUCE_THREADS 128
#define REDUCE_MAX_BLOCKS 128

__global__ void isingCuda_magnetizationKernel (unsigned int *g_idata, unsigned int *g_odata, unsigned int n) {
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*REDUCE_THREADS + threadIdx.x;
    unsigned int gridSize = gridDim.x*REDUCE_THREADS;
    
    unsigned int acc = 0;
    while (i < n) {
        acc += (unsigned int) bitCount(g_idata[i]);
        i += gridSize;
    }

    extern __shared__ unsigned int s[];
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
    if (tid == 0) g_odata[blockIdx.x] = s[tid];
}

int divideCeil(int x, int y) {
    return (x + y - 1) / y;
}

double isingCuda_bitCount(unsigned int *d_idata, int n) {
    // allocate arrays on device and host to store one float for each block
    int blocks = min(REDUCE_MAX_BLOCKS, divideCeil(n, REDUCE_THREADS));
    unsigned int h_odata[REDUCE_MAX_BLOCKS];
    unsigned int *d_odata;
    cudaMalloc((void**) &d_odata, blocks*sizeof(unsigned int));
    
    // partial reduction; each block generates one number
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = REDUCE_THREADS * sizeof(unsigned int);
    isingCuda_magnetizationKernel<<< blocks, REDUCE_THREADS, smemSize >>>(d_idata, d_odata, n);
    
    // copy result from device to host, and perform final reduction on CPU
    cudaMemcpy(h_odata, d_odata, blocks*sizeof(unsigned int), cudaMemcpyDeviceToHost);
    double gpu_result = 0;
    for (int i = 0; i < blocks; i++)  {
        gpu_result += h_odata[i];
    }
    cudaFree(d_odata);    
    return gpu_result;
}


// ----------------------------------------------------------------------------
// Cuda update implementation

typedef struct {
    int len, dim, nblocks;
    unsigned int *blocks;
    float h, T; // external field and temperature
    int parityTarget;
} IsingCudaParams;


__device__ inline int shouldFlipSpin(IsingCudaParams p, Rand48 &rng, int s, int m) {
    float dE = 2*s*(m + p.h);
    if (dE < 0)
        return 1;
    else {
#ifdef DETERMINISTIC
        float r = 0.1;
#else
        float r = (float)rand48_nextInt(rng) / (unsigned int)(1<<31);
#endif
        return exp(- dE / p.T) > r;
    }
}

__device__ inline void isingCuda_updateSite(IsingCudaParams p, Rand48 &rng, int ip) {
    int parity = 0;
    Bits128 acc = {0, 0, 0, 0};
    int lenp_d = 1;
    Bits128 n1 = bitsExpand(p.blocks[ip]);
    
    for (int d = 0; d < p.dim; d++) {
        int lenp = (d < 5) ? p.len/2 : p.len;
        int xp = (ip / lenp_d) % lenp;
        parity += (d < 5) ? 0 : xp;
        
        int dx2 = (xp+1+lenp)%lenp - xp;
        int dx0 = (xp-1+lenp)%lenp - xp;
        Bits128 n2 = bitsExpand(p.blocks[ip+dx2*lenp_d]);
        Bits128 n0 = bitsExpand(p.blocks[ip+dx0*lenp_d]);
        
        if (d < 5) {
            int shift = 4 << d; // 4, 8, 16, 32, 64
            acc = bitsAdd(acc, bitsMaskShiftL(bitsAdd(n1,n2), shift));
            acc = bitsAdd(acc, bitsMaskShiftR(bitsAdd(n1,n0), shift));
        }
        else {
            acc = bitsAdd(bitsAdd(n0,n2), acc);
        }
        
        lenp_d *= lenp;
    }
    
    int deltaMax = p.dim < 5 ? (1 << p.dim) : 32;
    int cube = p.blocks[ip];
    for (int delta = 0; delta < deltaMax; delta++) {
        if ((parity + bitCount(delta)) % 2 == p.parityTarget) {
            // m = total magnetization of neighbors; in range [-2 dim, +2 dim]
            int m = 2*(bitsPick4(acc, delta) - p.dim);
            // s = spin at this site (+/- 1)
            int s = 2*((cube >> delta) & 1) - 1;
            if (shouldFlipSpin(p, rng, s, m)) {
                cube ^= (1 << delta);
            }
        }
    }
    p.blocks[ip] = cube;
}

__global__ void isingCuda_update(IsingCudaParams p, Rand48 rng) {
    rand48_loadState(rng);

    unsigned int ip = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;
    while (ip < p.nblocks) {
        isingCuda_updateSite(p, rng, ip);
        ip += gridSize;
    }
    
    rand48_storeState(rng);
}


// ----------------------------------------------------------------------------
// Ising class interface
//

IsingCuda::IsingCuda(int len, int dim, float h, float T) : Ising(len, dim, h, T) {
    assert(len % 2 == 0);
    assert(dim <= 7);
    len = len;
    dim = dim;
    n = (int)powl(len, dim);
    nblocks = n >> min(5,dim);
    
    int nbytes = nblocks*sizeof(unsigned int);
    blocks = (unsigned int *)malloc(nbytes);
    cudaMalloc((void**)&d_blocks, nbytes);
    
    for (int i = 0; i < nblocks; i++) {
        blocks[i] = 0;
    }
    transferHostToDevice();
    
    rng = new Rand48();
    rng->init(GRID_DIM*BLOCK_DIM, 0); // initialize random numbers
}

IsingCuda::~IsingCuda() {
    free(blocks);
    cudaFree(d_blocks);
    rng->destroy();
    delete rng;
}


void IsingCuda::completeNeighborSum(int *sum) {
    assert(0==1);
}

void IsingCuda::update(int parityTarget) {
    IsingCudaParams p;
    p.len = len;
    p.dim = dim;
    p.nblocks = nblocks;
    p.blocks = d_blocks;
    p.h = h;
    p.T = T;
    p.parityTarget = parityTarget;

    int sharedBytes = 0;
    
    isingCuda_update <<<GRID_DIM, BLOCK_DIM, sharedBytes>>> (p, *rng);
}

double IsingCuda::magnetization() {
    return 2.0*isingCuda_bitCount(d_blocks, nblocks) - n;
}

void IsingCuda::transferHostToDevice() {
    cudaMemcpy(d_blocks, blocks, nblocks*sizeof(unsigned int), cudaMemcpyHostToDevice);
}

void IsingCuda::transferDeviceToHost() {
    cudaMemcpy(blocks, d_blocks, nblocks*sizeof(unsigned int), cudaMemcpyDeviceToHost);
}


// given index 'i' into the full lattice, return compressed index 'ip'
// and bit index 'delta'.
void IsingCuda::index(int i, int *ip, int *delta) {
    int len_d = 1;
    int lenp_d = 1;
    *ip = 0;
    *delta = 0;
    
    for (int d = 0; d < dim; d++) {
        int x = (i / len_d) % len;
        int xp = (d < 5) ? x/2 : x;
        int del = (d < 5) ? x%2 : 0;
        
        *delta += (del << d) ;
        *ip += xp*lenp_d;
        
        int lenp = (d < 5) ? len/2 : len;
        len_d *= len;
        lenp_d *= lenp;
    }
    
    assert(*ip < nblocks);
    assert(*delta < 32);
}

void IsingCuda::set(int i, int s) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    int mask = ~(1 << delta);
    blocks[ip] = (blocks[ip] & mask) | (s << delta);
}

int IsingCuda::get(int i) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    return (blocks[ip]>>delta) & 1;
}

// ----------------------------------------------------------------------------
// Cuda utility methods
//

void initCuda(int argc, char *argv[]) {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        fprintf(stderr, "No devices supporting CUDA.\n");
        exit(EXIT_FAILURE);
    }
    int dev = deviceCount - 1;
    if (argc > 1) {
        sscanf(argv[1], "%d", &dev);
    }
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major < 1) {
        fprintf(stderr, "Device %d does not support CUDA.\n", dev);
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Using device %d: %s\n", dev, deviceProp.name);
    cudaSetDevice(dev);
}
