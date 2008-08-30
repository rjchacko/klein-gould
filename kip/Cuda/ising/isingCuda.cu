#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "ising.h"


// ----------------------------------------------------------------------------
// Cuda implementation
//

typedef struct {
    // a lattice containing (n = len^dim) spins
    int len, dim, n;
    unsigned int *spins; // n spins stored using n/32 unsigned ints
    float h, T; // external field and temperature
} NVIsing;


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
