
typedef struct {
    float x, y, z, w;
} float4;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qcd.h"


int index(int z, int y, int x, int t) {
    z = (z + L3) % L3;
    y = (y + L2) % L2;
    x = (x + L1) % L1;
    t = (t + L0) % L0;
    return z*L2*L1*L0 + y*L1*L0 + x*L0 + t;
}

int coord_0(int idx) {
    return idx % L0;
}

int coord_1(int idx) {
    return (idx/L0) % L1;
}

int coord_2(int idx) {
    return (idx/(L0*L1)) % L2;
}

int coord_3(int idx) {
    return idx/(L0*L1*L2);
}

int nextIndex(int i, int dir) {
    int t = coord_0(i);
    int x = coord_1(i);
    int y = coord_2(i);
    int z = coord_3(i);
    switch (dir) {
        case 0: return index(z, y, x, t+1);
        case 1: return index(z, y, x+1, t-1);
        case 2: return index(z, y, x-2, t);
        case 3: return index(z, y+1, x+1, t);
        case 4: return index(z, y-2, x, t);
        case 5: return index(z+1, y+1, x, t);
        case 6: return index(z-2, y, x, t);
        default: return -1;
    }
}

void packGaugeField(float4 *res, float *gauge) {
    for (int i = 0; i < L; i++) {
        for (int dir = 0; dir < 4; dir++) {
            for (int j = 0; j < 5; j++) {
                float a1, a2, a3, a4;
                a1 = gauge[i*4*18 + dir*18 + 4*j + 0];
                a2 = gauge[i*4*18 + dir*18 + 4*j + 1];
                if (j < 4) {
                    a3 = gauge[i*4*18 + dir*18 + 4*j + 2];
                    a4 = gauge[i*4*18 + dir*18 + 4*j + 3];
                }
                else {
                    int i3 = nextIndex(i, 2*dir+0);
                    int i4 = nextIndex(i, 2*dir+1);
                    a3 = *(float *)&i3;
                    a4 = *(float *)&i4;
                }
                float4 f4 = {a1, a2, a3, a4};
                res[(dir*5+j)*L + i] = f4;
            }
        }
    }
}

void packSpinorField(float4 *res, float *spinor) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < 6; j++) {
            float a1 = spinor[i*24 + 4*j + 0];
            float a2 = spinor[i*24 + 4*j + 1];
            float a3 = spinor[i*24 + 4*j + 2];
            float a4 = spinor[i*24 + 4*j + 3];
            float4 f4 = {a1, a2, a3, a4};
            res[j*L + i] = f4;
        }
    }
}

void unpackSpinorField(float *res, float4 *spinorPacked) {
    for (int i = 0; i < L; i++) {
        if (0) {
            for (int j = 0; j < 6; j++) {
                float4 f4 = spinorPacked[j*L + i];
                res[i*(6*4) + j*(4) + 0] = f4.x;
                res[i*(6*4) + j*(4) + 1] = f4.y;
                res[i*(6*4) + j*(4) + 2] = f4.z;
                res[i*(6*4) + j*(4) + 3] = f4.w;
            }
        }
        else {
            for (int j = 0; j < 24; j++) {
                res[i*24 + j] = ((float *)spinorPacked)[j*L + i];
            }
        }
    }
}


void constructGaugeField(float *res) {
    for(int i = 0; i < L; i++) {
        for (int dir = 0; dir < 4; dir++) {
            for (int m = 0; m < 3; m++) {
                for (int n = 0; n < 3; n++) {
                    res[i*(4*3*3*2) + dir*(3*3*2) + m*(3*2) + n*(2) + 0] = 1;
                    res[i*(4*3*3*2) + dir*(3*3*2) + m*(3*2) + n*(2) + 1] = 0;
                }
            }
        }
    }
}


void constructSpinorField(float *res) {
    for(int i = 0; i < L; i++) {
        for (int s = 0; s < 4; s++) {
            for (int m = 0; m < 3; m++) {
                res[i*(4*3*2) + s*(3*2) + m*(2) + 0] = m;
                res[i*(4*3*2) + s*(3*2) + m*(2) + 1] = 0;
            }
        }
    }
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

void zero(float* a, int cnt) {
    for (int i = 0; i < cnt; i++)
        a[i] = 0;
}

void copy(float *dst, float *src, int cnt) {
    for (int i = 0; i < cnt; i++)
        dst[i] = src[i];
}

void sum(float *dst, float *a, float *b, int cnt) {
    for (int i = 0; i < cnt; i++)
        dst[i] = a[i] + b[i];
}

    
const float projector[8][4][4][2] =
{
    {
      {{1,0}, {0,0}, {0,0}, {0,1}},
      {{0,0}, {1,0}, {0,1}, {0,0}},
      {{0,0}, {0,-1}, {1,0}, {0,0}},
      {{0,-1}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {0,-1}},
      {{0,0}, {1,0}, {0,-1}, {0,0}},
      {{0,0}, {0,1}, {1,0}, {0,0}},
      {{0,1}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {-1,0}},
      {{0,0}, {1,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {1,0}, {0,0}},
      {{-1,0}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,0}, {1,0}},
      {{0,0}, {1,0}, {-1,0}, {0,0}},
      {{0,0}, {-1,0}, {1,0}, {0,0}},
      {{1,0}, {0,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,1}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {0,-1}},
      {{0,-1}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {0,1}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {0,-1}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {0,1}},
      {{0,1}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {0,-1}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {1,0}},
      {{1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {1,0}}
    },
    {
      {{1,0}, {0,0}, {-1,0}, {0,0}},
      {{0,0}, {1,0}, {0,0}, {-1,0}},
      {{-1,0}, {0,0}, {1,0}, {0,0}},
      {{0,0}, {-1,0}, {0,0}, {1,0}}
    }
};


void multiplySpinorByDiracProjector(float *res, int dir, float *spinorIn) {
    zero(res, 4*3*2);
    
    for (int s = 0; s < 4; s++) {
        for (int t = 0; t < 4; t++) {
            float projRe = projector[dir][s][t][0];
            float projIm = projector[dir][s][t][1];
            
            for (int m = 0; m < 3; m++) {
                float spinorRe = spinorIn[t*(3*2) + m*(2) + 0];
                float spinorIm = spinorIn[t*(3*2) + m*(2) + 1];
                res[s*(3*2) + m*(2) + 0] += projRe*spinorRe - projIm*spinorIm;
                res[s*(3*2) + m*(2) + 1] += projRe*spinorIm + projIm*spinorRe;
            }
        }
    }
}

void computeGold(float* res, float* gauge, float *spinor) {
    zero(res, L*4*3*2);
    
    for (int idxOut = 0; idxOut < L; idxOut++) {
        int idxIn = idxOut;
        for (int dir = 0; dir < 1; dir++) {
            float projectedSpinor[4*3*2], gaugedSpinor[4*3*2];
            float *gaugeMatrix = &gauge[idxIn*(4*3*3*2) + (dir/2)*(3*3*2)];
            multiplySpinorByDiracProjector(projectedSpinor, dir, &spinor[idxIn*(4*3*2)]);
            
            for (int s = 0; s < 4; s++) {
                if (dir % 2 == 0)
                    su3_mul(&gaugedSpinor[s*(3*2)], gaugeMatrix, &projectedSpinor[s*(3*2)]);
                else
                    su3_Tmul(&gaugedSpinor[s*(3*2)], gaugeMatrix, &projectedSpinor[s*(3*2)]);
            }
            
            sum(&res[idxIn*(4*3*2)], &res[idxIn*(4*3*2)], gaugedSpinor, 4*3*2);
            idxIn = nextIndex(idxIn, dir);
        }
    }
}


float dist(float x, float y) {
    return fabs(x - y);
}

void testSpinorField(float *spinor) {
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < 12; j++) {
            if (dist(spinor[i*(4*3*2) + j], 3) > 1e-6) {
                printf("doh1 %d %f %f\n", i, spinor[i*(4*3*2) + j], 3.0);
            }
        }
        for (int j = 12; j < 24; j++) {
            float target = (j%2 == 0) ? 3 : -3;
            if (dist(spinor[i*(4*3*2) + j], target) > 1e-6) {
                printf("doh2 %d %f %f\n", i, spinor[i*(4*3*2) + j], target);
            }
        }
    }
}

void printVector(float *v) {
    printf("{(%f %f) (%f %f) (%f %f)}\n", v[0], v[1], v[2], v[3], v[4], v[5]);
}


void printSpinorField(float *spinor) {
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[0*(4*3*2) + s*(3*2)]);
    }
    printf("...\n");
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[(L-1)*(4*3*2) + s*(3*2)]);
    }
    printf("\n");    
}