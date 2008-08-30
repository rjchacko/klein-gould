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


// count the number of one-bits in the binary representation of i.
int bitCount(unsigned int v) {
    unsigned int c;
    c = v - ((v >> 1) & 0x55555555);
    c = ((c >> 2) & 0x33333333) + (c & 0x33333333);
    c = ((c >> 4) + c) & 0x0F0F0F0F;
    c = ((c >> 8) + c) & 0x00FF00FF;
    c = ((c >> 16) + c) & 0x0000FFFF;
    return c;
}


// ------------------------------------------------------------------------------------------
// Bits128
//

// Bits128 should be thought of as a "128 bit unsigned int". the highest (127th) bit
// is the highest (31st) bit of a3, which explains why a3 appears on the left.
typedef struct {
    unsigned int a3, a2, a1, a0;
} Bits128;


// expand 8 low bits to be equally spaced among the 32 bits of an int
unsigned int expandFF(unsigned int x) {
    return
        x&1 |
        ((x>>1)&1)<<4 |
        ((x>>2)&1)<<8 |
        ((x>>3)&1)<<12 |
        ((x>>4)&1)<<16 |
        ((x>>5)&1)<<20 |
        ((x>>6)&1)<<24 |
        ((x>>7)&1)<<28;
}

// expand 32 bits of an int to be equally spaced among 128 bits of a Bits128
// structure. i.e., each input bit gets 4 bits of spacing in the return value.
Bits128 bitsExpand(int a) {
    Bits128 ret = {
        expandFF((a>>24) & 0xFF),
        expandFF((a>>16) & 0xFF),
        expandFF((a>>8) & 0xFF),
        expandFF((a>>0) & 0xFF)
    };
    return ret;
}

// add two Bits128 structs. note: overflow between 32 components not supported
Bits128 bitsAdd(Bits128 x, Bits128 y) {
    Bits128 ret = {
        x.a3+y.a3, x.a2+y.a2, x.a1+y.a1, x.a0+y.a0
    };
    return ret;
}

// shift input 'n' bits to the left, and then apply a 'mask'. specifically:
//  ret = (x << n) & (A B A B ... A B)
//
// where A and B represent strings of 'n' 1 and 0 bits, respectively
//
Bits128 bitsMaskShiftL(Bits128 x, int n) {
    Bits128 ret;
    unsigned int m;
    switch (n) {
        case 64:
            ret.a3 = x.a1;
            ret.a2 = x.a0;
            ret.a1 = 0;
            ret.a0 = 0;
            return ret;
        case 32:
            ret.a3 = x.a2;
            ret.a2 = 0;
            ret.a1 = x.a0;
            ret.a0 = 0;
            return ret;
        case 16:
            m = 0xffff0000;
            break;
        case 8:
            m = 0xff00ff00;
            break;
        case 4:
            m = 0xf0f0f0f0;
            break;
        default:
            printf("Illegal argument n=%d, bitsMaskShiftL\n", n);
            exit(1);
    }
    ret.a3 = (x.a3<<n)&m;
    ret.a2 = (x.a2<<n)&m;
    ret.a1 = (x.a1<<n)&m;
    ret.a0 = (x.a0<<n)&m;
    return ret;
}

// shift input 'n' bits to the right, and then apply a 'mask'. specifically:
//  ret = (x >> n) & (B A B A ... B A)
//
// where A and B represent strings of 'n' 1 and 0 bits, respectively
//
Bits128 bitsMaskShiftR(Bits128 x, int n) {
    Bits128 ret;
    unsigned int m;
    switch (n) {
        case 64:
            ret.a3 = 0;
            ret.a2 = 0;
            ret.a1 = x.a3;
            ret.a0 = x.a2;
            return ret;
        case 32:
            ret.a3 = 0;
            ret.a2 = x.a3;
            ret.a1 = 0;
            ret.a0 = x.a1;
            return ret;
        case 16:
            m = 0x0000ffff;
            break;
        case 8:
            m = 0x00ff00ff;
            break;
        case 4:
            m = 0x0f0f0f0f;
            break;
        default:
            printf("Illegal argument n=%d, bitsMaskShiftR\n", n);
            exit(1);
    }
    ret.a3 = (x.a3>>n)&m;
    ret.a2 = (x.a2>>n)&m;
    ret.a1 = (x.a1>>n)&m;
    ret.a0 = (x.a0>>n)&m;
    return ret;

}

// shift input 'n' bits to the right, and then pick the low 4 bits. specifically,
//  ret = (x >> n) & 0xf
//
int bitsPick4(Bits128 x, int n) {
    assert(n >= 0 && n < 32);
    int rshift = 4 * (n%8);
    switch (n / 8) {
    case 0: return (x.a0 >> rshift) & 0xf;
    case 1: return (x.a1 >> rshift) & 0xf;
    case 2: return (x.a2 >> rshift) & 0xf;
    case 3: return (x.a3 >> rshift) & 0xf;
    }
    printf("Illegal state, bitsPick\n");
    exit(1);
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
