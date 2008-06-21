#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <ctime>
#include <sys/time.h>

#include "qcd.h"


int compareFloats(float *a, float *b, int len, float epsilon) {
    for (int i = 0; i < len; i++) {
        float diff = a[i] - b[i];
        if (diff < 0)
            diff = -diff;
        if (diff > epsilon)
            return 0;
    }
    return 1;
}


struct timeval startTime;

void stopwatchStart() {
    gettimeofday(&startTime, NULL);
}

double stopwatchReadSeconds() {
    struct timeval endTime;
    gettimeofday( &endTime, 0);
    
    long ds = endTime.tv_sec - startTime.tv_sec;
    long dus = endTime.tv_usec - startTime.tv_usec;
    return ds + 0.000001*dus;
}


void printVector(float *v) {
    printf("{(%f %f) (%f %f) (%f %f)}\n", v[0], v[1], v[2], v[3], v[4], v[5]);
}

void printSpinor(float *spinor) {
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[s*(3*2)]);
    }
}

// X indexes the full lattice
void printSpinorElement(float *spinor, int X) {
    if (getOddBit(X) == 0)
        printSpinor(&spinor[(X/2)*(4*3*2)]);
    else
        printSpinor(&spinor[(X/2)*(4*3*2)+Nh*spinorSiteSize]);
}

void printGauge(float *gauge) {
    for (int m = 0; m < 3; m++) {
        printVector(&gauge[m*(3*2)]);
    }
}

// X indexes the full lattice
void printGaugeElement(float *gauge, int X) {
    if (getOddBit(X) == 0)
        printGauge(&gauge[(X/2)*gaugeSiteSize]);
    else
        printGauge(&gauge[(X/2+Nh)*gaugeSiteSize]);
}


// given a "half index" i into either an even or odd half lattice (corresponding
// to oddBit = {0, 1}), returns the corresponding full lattice index.
int fullLatticeIndex(int i, int oddBit) {
    int boundaryCrossings = i/L1h + i/(L2*L1h) + i/(L3*L2*L1h);
    return 2*i + (boundaryCrossings + oddBit) % 2;
}

// returns 0 or 1 if the full lattice index X is even or odd
int getOddBit(int X) {
    int x4 = X/(L3*L2*L1);
    int x3 = (X/(L2*L1)) % L3;
    int x2 = (X/L1) % L2;
    int x1 = X % L1;
    return (x4+x3+x2+x1) % 2;
}

void accumulateConjugateProduct(float *a, float *b, float *c, float sign) {
    a[0] += sign * (b[0]*c[0] - b[1]*c[1]);
    a[1] -= sign * (b[0]*c[1] + b[1]*c[0]);
}

// given first two rows (u,v) of SU(3) matrix mat, reconstruct the third row
// as the cross product of the conjugate vectors: w = u* x v*
// 
void su3_reconstruct(float *mat) {
    float *u = &mat[0*(3*2)];
    float *v = &mat[1*(3*2)];
    float *w = &mat[2*(3*2)];
    accumulateConjugateProduct(w+0*(2), u+1*(2), v+2*(2), +1);
    accumulateConjugateProduct(w+0*(2), u+2*(2), v+1*(2), -1);
    accumulateConjugateProduct(w+1*(2), u+2*(2), v+0*(2), +1);
    accumulateConjugateProduct(w+1*(2), u+0*(2), v+2*(2), -1);
    accumulateConjugateProduct(w+2*(2), u+0*(2), v+1*(2), +1);
    accumulateConjugateProduct(w+2*(2), u+1*(2), v+0*(2), -1);
}


void applyGaugeFieldScaling(float **gauge) {
    // Apply spatial scaling factor (u0) to spatial links
    for (int d=0; d<3; d++) {
        for (int i=0; i<gaugeSiteSize*N; i++) {
            gauge[d][i] /= SPATIAL_SCALING;
        }
    }
    
    // Apply boundary conditions to temporal links
    int surface = L1h*L2*L3;
    for (int j=Nh-surface; j<Nh; j++) {
        for (int i=0; i<gaugeSiteSize; i++) {
            gauge[3][j*gaugeSiteSize+i] *= TIME_SYMMETRY;
            gauge[3][(Nh+j)*gaugeSiteSize+i] *= TIME_SYMMETRY;
        }
    }
}

void constructUnitGaugeField(float **res) {
    float *resOdd[4], *resEven[4];
    for (int dir = 0; dir < 4; dir++) {  
        resEven[dir] = res[dir];
        resOdd[dir]  = res[dir]+Nh*gaugeSiteSize;
    }
    
    for (int dir = 0; dir < 4; dir++) {
        for(int i = 0; i < Nh; i++) {
            for (int m = 0; m < 3; m++) {
                for (int n = 0; n < 3; n++) {
                    resEven[dir][i*(3*3*2) + m*(3*2) + n*(2) + 0] = (m==n) ? 1 : 0;
                    resEven[dir][i*(3*3*2) + m*(3*2) + n*(2) + 1] = 0.0;
                    resOdd[dir][i*(3*3*2) + m*(3*2) + n*(2) + 0] = (m==n) ? 1 : 0;
                    resOdd[dir][i*(3*3*2) + m*(3*2) + n*(2) + 1] = 0.0;
                }
            }
        }
    }
    
    applyGaugeFieldScaling(res);
}

void constructGaugeField(float **res) {
    float *resOdd[4], *resEven[4];
    for (int dir = 0; dir < 4; dir++) {  
        resEven[dir] = res[dir];
        resOdd[dir]  = res[dir]+Nh*gaugeSiteSize;
    }
    
    for (int dir = 0; dir < 4; dir++) {
        for (int i = 0; i < Nh; i++) {
            for (int m = 0; m < 2; m++) {
                for (int n = 0; n < 3; n++) {
                    resEven[dir][i*(3*3*2) + m*(3*2) + n*(2) + 0] = rand() / (float)RAND_MAX;
                    resEven[dir][i*(3*3*2) + m*(3*2) + n*(2) + 1] = rand() / (float)RAND_MAX;
                    resOdd[dir][i*(3*3*2) + m*(3*2) + n*(2) + 0] = rand() / (float)RAND_MAX;
                    resOdd[dir][i*(3*3*2) + m*(3*2) + n*(2) + 1] = rand() / (float)RAND_MAX;                    
                }
            }
            su3_reconstruct(&resEven[dir][i*(3*3*2)]);
            su3_reconstruct(&resOdd[dir][i*(3*3*2)]);
        }
    }
    
    applyGaugeFieldScaling(res);
}

void constructPointSpinorField(float *res, int i0, int s0, int c0) {
    float *resEven = res;
    float *resOdd = res + Nh*spinorSiteSize;
    
    for(int i = 0; i < Nh; i++) {
        for (int s = 0; s < 4; s++) {
            for (int m = 0; m < 3; m++) {
                resEven[i*(4*3*2) + s*(3*2) + m*(2) + 0] = 0;
                resEven[i*(4*3*2) + s*(3*2) + m*(2) + 1] = 0;
                resOdd[i*(4*3*2) + s*(3*2) + m*(2) + 0] = 0;
                resOdd[i*(4*3*2) + s*(3*2) + m*(2) + 1] = 0;
                if (s == s0 && m == c0) {
                    if (fullLatticeIndex(i, 0) == i0)
                        resEven[i*(4*3*2) + s*(3*2) + m*(2) + 0] = 1;
                    if (fullLatticeIndex(i, 1) == i0)
                        resOdd[i*(4*3*2) + s*(3*2) + m*(2) + 0] = 1;
                }
            }
        }
    }
}

void constructSpinorField(float *res) {
    for(int i = 0; i < N; i++) {
        for (int s = 0; s < 4; s++) {
            for (int m = 0; m < 3; m++) {
                res[i*(4*3*2) + s*(3*2) + m*(2) + 0] = rand() / (float)RAND_MAX;
                res[i*(4*3*2) + s*(3*2) + m*(2) + 1] = rand() / (float)RAND_MAX;
            }
        }
    }
}


void applyGamma5(float *out, float *in, int sites) {
    for (int i=0; i<sites*spinorSiteSize; i+=spinorSiteSize) {
        for (int j=0; j<spinorSiteSize/2; j++) 
            out[i+j] = in[i+j];
        for (int j=0; j<spinorSiteSize/2; j++) 
            out[i+j+spinorSiteSize/2] = -in[i+j+spinorSiteSize/2];
    }
}


int main(int argc, char **argv) {
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

//    blasTest();
//    axpbyTest();
    dslashTest();
//    cgTest();
}
