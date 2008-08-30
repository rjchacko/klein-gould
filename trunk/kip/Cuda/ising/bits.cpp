#include "ising.h"


// if CUDA_INCLUDE is defined, then compile the following routines for execution
// on device

#ifdef CUDA_INCLUDE
#define PREFIX __device__
#define MAYBE_ASSERT(x)
#else
#include <assert.h>
#define PREFIX
#define MAYBE_ASSERT(x) assert(x)
#endif




// count the number of one-bits in the binary representation of i.
PREFIX int bitCount(unsigned int v) {
    unsigned int c;
    c = v - ((v >> 1) & 0x55555555);
    c = ((c >> 2) & 0x33333333) + (c & 0x33333333);
    c = ((c >> 4) + c) & 0x0F0F0F0F;
    c = ((c >> 8) + c) & 0x00FF00FF;
    c = ((c >> 16) + c) & 0x0000FFFF;
    return c;
}


// expand 8 low bits to be equally spaced among the 32 bits of an int
PREFIX unsigned int expandFF(unsigned int x) {
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
PREFIX Bits128 bitsExpand(int a) {
    Bits128 ret = {
        expandFF((a>>24) & 0xFF),
        expandFF((a>>16) & 0xFF),
        expandFF((a>>8) & 0xFF),
        expandFF((a>>0) & 0xFF)
    };
    return ret;
}

// add two Bits128 structs. note: overflow between 32 components not supported
PREFIX Bits128 bitsAdd(Bits128 x, Bits128 y) {
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
PREFIX Bits128 bitsMaskShiftL(Bits128 x, int n) {
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
            MAYBE_ASSERT(0);
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
PREFIX Bits128 bitsMaskShiftR(Bits128 x, int n) {
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
            MAYBE_ASSERT(0);
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
PREFIX int bitsPick4(Bits128 x, int n) {
    MAYBE_ASSERT(n >= 0 && n < 32);
    int rshift = 4 * (n%8);
    switch (n / 8) {
    case 0: return (x.a0 >> rshift) & 0xf;
    case 1: return (x.a1 >> rshift) & 0xf;
    case 2: return (x.a2 >> rshift) & 0xf;
    case 3: return (x.a3 >> rshift) & 0xf;
    }
    MAYBE_ASSERT(0);
    return -1;
}
