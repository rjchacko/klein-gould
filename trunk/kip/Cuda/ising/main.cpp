#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

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

int pentacubeParity[32];

void initPentacubeParity() {
    pentacubeParity[0] = 0;
    int m = 1;
    for (int i = 1; i < 32; i++) {
        if (i == 2*m)
            m = i;
        pentacubeParity[i] = 1 - pentacubeParity[i-m];
    }
}


double T = 0;
double h = 0.1;

int flipSpin(int s, int m) {
    float dE = 2*s*(m + h);
    if (dE < 0)
        return 1;
    else {
        float r = 0.5; // rand() / (float)RAND_MAX;
        return exp(- dE / T) > r;
    }
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
// Ising1
//

typedef struct {
    int len, dim, n;
    int *spins;
    int *sum;
} Ising1;


Ising1 createIsing1(int len, int dim) {
    Ising1 ret;
    
    ret.len = len;
    ret.dim = dim;
    ret.n = (int)powl(len, dim);
    ret.spins = (int *)malloc(ret.n*sizeof(int));
    ret.sum = (int *)malloc(ret.n*sizeof(int));
    
    for (int i = 0; i < ret.n; i++) {
        ret.spins[i] = ret.sum[i] = 0;
    }
    return ret;
}

void set1(Ising1 self, int i, int s) {
    self.spins[i] = s;
}

int get1(Ising1 self, int i) {
    return self.spins[i];
}

void init1(Ising1 self) {
    for (int i = 0; i < self.n; i++) {
        set1(self, i, (i%5+i%4)%2);
    }
}

void print1(Ising1 self) {
    for (int y = 0; y < self.len; y++) {
        for (int x = 0; x < self.len; x++) {
            printf("%d ", get1(self, y*self.len+x));
        }
        printf("\n");
    }
    printf("\n");
}

void printNeighbors1(Ising1 self) {
    for (int y = 0; y < self.len; y++) {
        for (int x = 0; x < self.len; x++) {
            printf("%d ", self.sum[y*self.len+x]);
        }
        printf("\n");
    }
    printf("\n");
}

int neighborSum1(Ising1 self, int i) {
    int acc = 0;
    int len = self.len;
    int len_d = 1;
    for (int d = 0; d < self.dim; d++) {
        int x = (i / len_d) % len;
        for (int dir = -1; dir <= 1; dir += 2) {
            int xp = (x + dir + len) % len;
            int dx = xp - x;
            acc += get1(self, i + dx * len_d);
        }
        len_d *= len;
    }
    return acc;
}

void completeNeighborSum1(Ising1 self) {
    for (int i = 0; i < self.n; i++) {
        self.sum[i] = neighborSum1(self, i);
    }
}


int indexParity1(Ising1 self, int i) {
    int acc = 0;
    int len = self.len;
    int len_d = 1;
    for (int d = 0; d < self.dim; d++) {
        int x = (i / len_d) % len;
        acc += x;
        len_d *= len;
    }
    return acc % 2;
}


// update even or odd sites if parityTarget is 0 or 1.
void update1(Ising1 self, int parityTarget) {
    for (int i = 0; i < self.n; i++) {
        if (indexParity1(self, i) == parityTarget) {
            int m = 2 * (neighborSum1(self, i) - self.dim);
            int s = 2 * get1(self, i) - 1;
            if (flipSpin(s, m)) {
                set1(self, i, (1-s)/2);
            }
        }
    }
}

// ------------------------------------------------------------------------------------------
// Ising2
//


typedef struct {
    // a lattice containing (n = len^dim) spins
    int len, dim, n;
    // stored in groups of 32 bits
    unsigned int *spins; // n/32 elements
    int *sum;
} Ising2;


Ising2 createIsing2(int len, int dim) {
    Ising2 ret;
    
    ret.len = len;
    ret.dim = dim;
    ret.n = (int)powl(len, dim);
    ret.spins = (unsigned int *)malloc((ret.n/32)*sizeof(unsigned int));
    ret.sum = (int *)malloc(ret.n*sizeof(int));
    
    for (int ip = 0; ip < ret.n/32; ip++)
        ret.spins[ip] = 0;
    for (int i = 0; i < ret.n; i++) {
        ret.sum[i] = 0;
    }
    return ret;
}

// given index 'i' into the full lattice, return compressed index 'ip'
// and bit index 'delta'.
void index2(Ising2 self, int i, int *ip, int *delta) {
    int len = self.len;
    int len_d = 1;
    int lenp_d = 1;
    *ip = 0;
    *delta = 0;
    
    for (int d = 0; d < self.dim; d++) {
        int x = (i / len_d) % len;
        int xp = (d < 5) ? x/2 : x;
        int del = (d < 5) ? x%2 : 0;
        
        *delta += (del << d) ;
        *ip += xp*lenp_d;
      
        int lenp = (d < 5) ? len/2 : len;
        len_d *= len;
        lenp_d *= lenp;
    }
}

// given compressed index 'ip' and bit index 'delta', return index 'i' into the
// full lattice
int reverseIndex2(Ising2 self, int ip, int delta) {
    int len = self.len;
    int len_d = 1;
    int lenp_d = 1;
    int i = 0;
    
    for (int d = 0; d < self.dim; d++) {
      int lenp = (d < 5) ? len/2 : len;
      int xp = (ip / lenp_d) % lenp;
      int x = (d < 5) ? 2*xp + ((unsigned int)delta >> d)%2 : xp;
      i += x*len_d;
      
      len_d *= len;
      lenp_d *= lenp;
    }
    return i;
}


void set2(Ising2 self, int i, int s) {
    int ip, delta;
    index2(self, i, &ip, &delta);
    int mask = ~(1 << delta);
    self.spins[ip] = (self.spins[ip] & mask) | (s << delta);
}

int get2(Ising2 self, int i) {
    int ip, delta;
    index2(self, i, &ip, &delta);
    return (self.spins[ip]>>delta) & 1;
}

void init2(Ising2 self) {
    for (int i = 0; i < self.n; i++) {
        set2(self, i, (i%5+i%4)%2);
    }
}

void print2(Ising2 self) {
    for (int y = 0; y < self.len; y++) {
        for (int x = 0; x < self.len; x++) {
            printf("%d ", get2(self, y*self.len+x));
        }
        printf("\n");
    }
    printf("\n");
}

void printNeighbors2(Ising2 self) {
    for (int y = 0; y < self.len; y++) {
        for (int x = 0; x < self.len; x++) {
            printf("%d ", self.sum[y*self.len+x]);
        }
        printf("\n");
    }
    printf("\n");
}


void completeNeighborSum2(Ising2 self) {
    for (int ip = 0; ip < self.n/32; ip++) {
        Bits128 acc = {0, 0, 0, 0};
        int len = self.len;
        int lenp_d = 1;
        Bits128 n1 = bitsExpand(self.spins[ip]);
        
        for (int d = 0; d < self.dim; d++) {
            int lenp = (d < 5) ? len/2 : len;
            int xp = (ip / lenp_d) % lenp;
            
            int dx2 = (xp+1+lenp)%lenp - xp;
            int dx0 = (xp-1+lenp)%lenp - xp;
            Bits128 n2 = bitsExpand(self.spins[ip+dx2*lenp_d]);
            Bits128 n0 = bitsExpand(self.spins[ip+dx0*lenp_d]);
            
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
        
        int deltaMax = min(32, (int)powl(2, self.dim));
        for (int delta = 0; delta < deltaMax; delta++) {
            int i = reverseIndex2(self, ip, delta);
            self.sum[i] = bitsPick4(acc, delta);
        }
    }
}



void update2(Ising2 self, int parityTarget) {
    for (int ip = 0; ip < self.n/32; ip++) {
        int parity = 0;
        Bits128 acc = {0, 0, 0, 0};
        int len = self.len;
        int lenp_d = 1;
        Bits128 n1 = bitsExpand(self.spins[ip]);
        
        for (int d = 0; d < self.dim; d++) {
            int lenp = (d < 5) ? len/2 : len;
            int xp = (ip / lenp_d) % lenp;
            parity += (d < 5) ? 0 : xp;
            
            int dx2 = (xp+1+lenp)%lenp - xp;
            int dx0 = (xp-1+lenp)%lenp - xp;
            Bits128 n2 = bitsExpand(self.spins[ip+dx2*lenp_d]);
            Bits128 n0 = bitsExpand(self.spins[ip+dx0*lenp_d]);
            
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
        
        int deltaMax = min(32, (int)powl(2, self.dim));
        int cube = self.spins[ip];
        for (int delta = 0; delta < deltaMax; delta++) {
            if ((parity + pentacubeParity[delta]) % 2 == parityTarget) {
                int m = 2*(bitsPick4(acc, delta) - self.dim);
                int s = 2*((cube >> delta) & 1) - 1;
                if (flipSpin(s, m)) {
                    cube ^= (1 << delta);
                }
            }
        }
        self.spins[ip] = cube;
    }
}



// ------------------------------------------------------------------------------------------
// Main
//


void isingEqual(Ising1 ising1, Ising2 ising2) {
    assert(ising1.n == ising2.n);
    for (int i = 0; i < ising1.n; i++) {
        if (get1(ising1, i) != get2(ising2, i)) {
            printf("Lattices are not equal at %d\n", i);
            return;
        }
    }
    printf("Lattices equal\n");
}


int main (int argc, char * const argv[]) {
    initPentacubeParity();
    
    Ising1 ising1 = createIsing1(8, 7);
    init1(ising1);
//    srand(0);
    update1(ising1, 0);
    update1(ising1, 1);
    print1(ising1);
    
//    completeNeighborSum1(ising1);
//    printNeighbors1(ising1);

    Ising2 ising2 = createIsing2(8, 7);
    init2(ising2);
//    srand(0);
    update2(ising2, 0);
    update2(ising2, 1);
    print2(ising2);
    
//    completeNeighborSum2(ising2);
//    printNeighbors2(ising2);
    
    isingEqual(ising1, ising2);
    return 0;
}
