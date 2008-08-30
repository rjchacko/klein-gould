#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "ising.h"


// ------------------------------------------------------------------------------------------
// Utilities
//

int min(int a, int b) {
    return (a < b) ? a : b;
}



// ------------------------------------------------------------------------------------------
// Ising2
//

Ising2::Ising2(int len, int dim, float h, float T) : Ising(len, dim, h, T) {
    assert(len % 2 == 0);
    assert(dim <= 7);
    len = len;
    dim = dim;
    n = (int)powl(len, dim);
    nblocks = n >> min(5,dim);
    blocks = (unsigned int *)malloc(nblocks*sizeof(unsigned int));
    
    for (int i = 0; i < nblocks; i++) {
        blocks[i] = 0;
    }
}

Ising2::~Ising2() {
    free(blocks);
}


// given index 'i' into the full lattice, return compressed index 'ip'
// and bit index 'delta'.
void Ising2::index(int i, int *ip, int *delta) {
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

// given compressed index 'ip' and bit index 'delta', return index 'i' into the
// full lattice
int Ising2::reverseIndex(int ip, int delta) {
    int len_d = 1;
    int lenp_d = 1;
    int i = 0;
    
    for (int d = 0; d < dim; d++) {
      int lenp = (d < 5) ? len/2 : len;
      int xp = (ip / lenp_d) % lenp;
      int x = (d < 5) ? 2*xp + ((unsigned int)delta >> d)%2 : xp;
      i += x*len_d;
      
      len_d *= len;
      lenp_d *= lenp;
    }
    
    assert(i < n);
    return i;
}


void Ising2::set(int i, int s) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    int mask = ~(1 << delta);
    blocks[ip] = (blocks[ip] & mask) | (s << delta);
}

int Ising2::get(int i) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    return (blocks[ip]>>delta) & 1;
}

void Ising2::completeNeighborSum(int *sum) {
    for (int ip = 0; ip < nblocks; ip++) {
        Bits128 acc = {0, 0, 0, 0};
        int lenp_d = 1;
        Bits128 n1 = bitsExpand(blocks[ip]);
        
        for (int d = 0; d < dim; d++) {
            int lenp = (d < 5) ? len/2 : len;
            int xp = (ip / lenp_d) % lenp;
            
            int dx2 = (xp+1+lenp)%lenp - xp;
            int dx0 = (xp-1+lenp)%lenp - xp;
            Bits128 n2 = bitsExpand(blocks[ip+dx2*lenp_d]);
            Bits128 n0 = bitsExpand(blocks[ip+dx0*lenp_d]);
            
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
        
        int deltaMax = min(32, (int)powl(2, dim));
        for (int delta = 0; delta < deltaMax; delta++) {
            int i = reverseIndex(ip, delta);
            sum[i] = bitsPick4(acc, delta);
        }
    }
}

void Ising2::update(int parityTarget) {
    for (int ip = 0; ip < nblocks; ip++) {
        int parity = 0;
        Bits128 acc = {0, 0, 0, 0};
        int lenp_d = 1;
        Bits128 n1 = bitsExpand(blocks[ip]);
        
        for (int d = 0; d < dim; d++) {
            int lenp = (d < 5) ? len/2 : len;
            int xp = (ip / lenp_d) % lenp;
            parity += (d < 5) ? 0 : xp;
            
            int dx2 = (xp+1+lenp)%lenp - xp;
            int dx0 = (xp-1+lenp)%lenp - xp;
            Bits128 n2 = bitsExpand(blocks[ip+dx2*lenp_d]);
            Bits128 n0 = bitsExpand(blocks[ip+dx0*lenp_d]);
            
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
        
        int deltaMax = min(32, (int)powl(2, dim));
        int cube = blocks[ip];
        for (int delta = 0; delta < deltaMax; delta++) {
            if ((parity + bitCount(delta)) % 2 == parityTarget) {
                int m = 2*(bitsPick4(acc, delta) - dim);
                int s = 2*((cube >> delta) & 1) - 1;
                if (shouldFlipSpin(s, m)) {
                    cube ^= (1 << delta);
                }
            }
        }
        blocks[ip] = cube;
    }
}
