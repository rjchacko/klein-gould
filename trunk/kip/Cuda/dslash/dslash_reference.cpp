#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qcd.h"


void zero(float* a, int cnt) {
    for (int i = 0; i < cnt; i++)
        a[i] = 0;
}

void sum(float *dst, float *a, float *b, int cnt) {
    for (int i = 0; i < cnt; i++)
        dst[i] = a[i] + b[i];
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

// i represents a "half index" into an even or odd "half lattice".
// when oddBit={0,1} the half lattice is {even,odd}.
// 
// the displacements, such as dx, refer to the full lattice coordinates. 
//
// neighborIndex() takes a "half index", displaces it, and returns the
// new "half index", which can be an index into either the even or odd lattices.
// displacements of magnitude one always interchange odd and even lattices.
//
int neighborIndex(int i, int oddBit, int dx4, int dx3, int dx2, int dx1) {
    int X = fullLatticeIndex(i, oddBit);
    int x4 = X/(L3*L2*L1);
    int x3 = (X/(L2*L1)) % L3;
    int x2 = (X/L1) % L2;
    int x1 = X % L1;

    // assert (oddBit == (x+y+z+t)%2);
    
    x4 = (x4+dx4+L4) % L4;
    x3 = (x3+dx3+L3) % L3;
    x2 = (x2+dx2+L2) % L2;
    x1 = (x1+dx1+L1) % L1;
    
    return (x4*(L3*L2*L1) + x3*(L2*L1) + x2*(L1) + x1) / 2;
}

float *gaugeLink(int i, int dir, int oddBit, float **gaugeEven, float **gaugeOdd) {
    float **gaugeField;
    int j;
    
    if (dir % 2 == 0) {
        j = i;
        gaugeField = (oddBit ? gaugeOdd : gaugeEven);
    }
    else {
        switch (dir) {
            case 1: j = neighborIndex(i, oddBit, 0, 0, 0, -1); break;
            case 3: j = neighborIndex(i, oddBit, 0, 0, -1, 0); break;
            case 5: j = neighborIndex(i, oddBit, 0, -1, 0, 0); break;
            case 7: j = neighborIndex(i, oddBit, -1, 0, 0, 0); break;
            default: j = -1; break;
        }
        gaugeField = (oddBit ? gaugeEven : gaugeOdd);
    }

    return &gaugeField[dir/2][j*(3*3*2)];
}

float *spinorNeighbor(int i, int dir, int oddBit, float *spinorField) {
    int j;
    switch (dir) {
        case 0: j = neighborIndex(i, oddBit, 0, 0, 0, +1); break;
        case 1: j = neighborIndex(i, oddBit, 0, 0, 0, -1); break;
        case 2: j = neighborIndex(i, oddBit, 0, 0, +1, 0); break;
        case 3: j = neighborIndex(i, oddBit, 0, 0, -1, 0); break;
        case 4: j = neighborIndex(i, oddBit, 0, +1, 0, 0); break;
        case 5: j = neighborIndex(i, oddBit, 0, -1, 0, 0); break;
        case 6: j = neighborIndex(i, oddBit, +1, 0, 0, 0); break;
        case 7: j = neighborIndex(i, oddBit, -1, 0, 0, 0); break;
        default: j = -1; break;
    }
    
    return &spinorField[j*(4*3*2)];
}

void dot(float* res, float* a, float* b) {
    res[0] = res[1] = 0;
    for (int m = 0; m < 3; m++) {
        float a_re = a[2*m+0];
        float a_im = a[2*m+1];
        float b_re = b[2*m+0];
        float b_im = b[2*m+1];
        res[0] += a_re * b_re - a_im * b_im;
        res[1] += a_re * b_im + a_im * b_re;
    }
}

void su3_transpose(float *res, float *mat) {
    for (int m = 0; m < 3; m++) {
        for (int n = 0; n < 3; n++) {
            res[m*(3*2) + n*(2) + 0] = + mat[n*(3*2) + m*(2) + 0];
            res[m*(3*2) + n*(2) + 1] = - mat[n*(3*2) + m*(2) + 1];
        }
    }
}

void su3_mul(float *res, float *mat, float *vec) {
    for (int n = 0; n < 3; n++) {
        dot(&res[n*(2)], &mat[n*(3*2)], vec);
    }
}

void su3_Tmul(float *res, float *mat, float *vec) {
    float matT[3*3*2];
    su3_transpose(matT, mat);
    su3_mul(res, matT, vec);
}

const float projector[8][4][4][2] = {
    {
      {{1,0}, {0,0}, {0,0}, {0,-1}},
      {{0,0}, {1,0}, {0,-1}, {0,0}},
      {{0,0}, {0,1}, {1,0}, {0,0}},
      {{0,1}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {0,1}},
      {{0,0}, {1,0}, {0,1}, {0,0}},
      {{0,0}, {0,-1}, {1,0}, {0,0}},
      {{0,-1}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {1,0}},
      {{0,0}, {1,0}, {-1,0}, {0,0}},
      {{0,0}, {-1,0}, {1,0}, {0,0}},
      {{1,0}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {-1,0}},
      {{0,0}, {1,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {1,0}, {0,0}},
      {{-1,0}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,-1}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {0,1}},
      {{0,1}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {0,-1}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,1}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {0,-1}},
      {{0,-1}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {0,1}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {-1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {-1,0}},
      {{-1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {-1,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {1,0}},
      {{1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {1,0}}
    }
};


void multiplySpinorByDiracProjector(float *res, float *spinorIn, int dir, int daggerBit) {
    zero(res, 4*3*2);
    
    int projIdx = !daggerBit ? dir : (dir + (1 - 2*(dir%2)));
    
    for (int s = 0; s < 4; s++) {
        for (int t = 0; t < 4; t++) {
            float projRe = projector[projIdx][s][t][0];
            float projIm = projector[projIdx][s][t][1];
            
            for (int m = 0; m < 3; m++) {
                float spinorRe = spinorIn[t*(3*2) + m*(2) + 0];
                float spinorIm = spinorIn[t*(3*2) + m*(2) + 1];
                res[s*(3*2) + m*(2) + 0] += projRe*spinorRe - projIm*spinorIm;
                res[s*(3*2) + m*(2) + 1] += projRe*spinorIm + projIm*spinorRe;
            }
        }
    }
}


// ---------------------------------------------------------------------------------------
// dslashReference()
//
// calculates the forward "d-slash" operation, given gauge field 'gauge' and
// spinor field 'spinor'. this function is intended to be a reference implementation for
// the much faster CUDA kernel.
//
// indices, varying slowest to fastest, are,
//   spinor field:  (z, y, x, t, spinor idx, color idx, complex idx)
//   gauge field:   (dir) (z, y, x, t, color row, color column, complex idx)
//
// constants (such as lattice lengths) are given in the file 'qcd.h'
// 
// if oddBit is zero/one then the even/odd spinor sites will be updated.
//
// if daggerBit is zero/one then perform dslash without/with dagger operator
//
void dslashReference(float *res, float **gaugeEven, float **gaugeOdd, float *spinorField, int oddBit, int daggerBit) {
    zero(res, Nh*4*3*2);
    
    for (int i = 0; i < Nh; i++) {
        for (int dir = 0; dir < 8; dir++) {
            float *gauge = gaugeLink(i, dir, oddBit, gaugeEven, gaugeOdd);
            float *spinor = spinorNeighbor(i, dir, oddBit, spinorField);
            
            float projectedSpinor[4*3*2], gaugedSpinor[4*3*2];
            multiplySpinorByDiracProjector(projectedSpinor, spinor, dir, daggerBit);
            
            for (int s = 0; s < 4; s++) {
                if (dir % 2 == 0)
                    su3_mul(&gaugedSpinor[s*(3*2)], gauge, &projectedSpinor[s*(3*2)]);
                else
                    su3_Tmul(&gaugedSpinor[s*(3*2)], gauge, &projectedSpinor[s*(3*2)]);
            }
            
            sum(&res[i*(4*3*2)], &res[i*(4*3*2)], gaugedSpinor, 4*3*2);
        }
    }
}



// ---------------------------------------------------------------------------------------
// helper functions
//

void printVector(float *v) {
    printf("{(%f %f) (%f %f) (%f %f)}\n", v[0], v[1], v[2], v[3], v[4], v[5]);
}

void printSpinor(float *spinor) {
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[s*(3*2)]);
    }
}

// X indexes the full lattice
void printSpinorElement(float *spinorEven, float *spinorOdd, int X) {
    if (getOddBit(X) == 0)
        printSpinor(&spinorEven[(X/2)*(4*3*2)]);
    else
        printSpinor(&spinorOdd[(X/2)*(4*3*2)]);
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

void constructUnitGaugeField(float **resEven, float **resOdd) {
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
}

void constructGaugeField(float **resEven, float **resOdd) {
    for (int dir = 0; dir < 4; dir++) {
        for(int i = 0; i < Nh; i++) {
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
}

void constructPointSpinorField(float *resEven, float *resOdd, int i0, int s0, int c0) {
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
    for(int i = 0; i < Nh; i++) {
        for (int s = 0; s < 4; s++) {
            for (int m = 0; m < 3; m++) {
                res[i*(4*3*2) + s*(3*2) + m*(2) + 0] = rand() / (float)RAND_MAX;
                res[i*(4*3*2) + s*(3*2) + m*(2) + 1] = rand() / (float)RAND_MAX;
            }
        }
    }
}

