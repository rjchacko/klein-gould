#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "ising.h"


// ----------------------------------------------------------------------------
// Cuda implementation
//

#define GRID_DIM 128
#define BLOCK_DIM 128

typedef struct {
    int len, dim, nblocks;
    unsigned int *blocks;
    float h, T; // external field and temperature
    int parityTarget;
} IsingCudaParams;


// include functions from "bits.cpp" as CUDA device functions
#define CUDA_INCLUDE
#include "bits.cpp"


__device__ int shouldFlipSpin(IsingCudaParams p, int s, int m) {
    float dE = 2*s*(m + p.h);
    if (dE < 0)
        return 1;
    else {
        float r = 0.5; // rand() / (float)RAND_MAX;
        return exp(- dE / p.T) > r;
    }
}

__device__ void isingCuda_updateSite(IsingCudaParams p, int ip) {
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
            if (shouldFlipSpin(p, s, m)) {
                cube ^= (1 << delta);
            }
        }
    }
    p.blocks[ip] = cube;
}

__global__ void isingCuda_update(IsingCudaParams p) {
    unsigned int ip = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;
    while (ip < p.nblocks) {
        isingCuda_updateSite(p, ip);
        ip += gridSize;
    }
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
}

IsingCuda::~IsingCuda() {
    free(blocks);
    cudaFree(d_blocks);
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
    isingCuda_update <<<GRID_DIM, BLOCK_DIM, sharedBytes>>> (p);
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