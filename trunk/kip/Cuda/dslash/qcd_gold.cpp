
typedef struct {
    float x, y, z, w;
} float4;

#include <stdio.h>
#include <stdlib.h>
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
        for (int j = 0; j < 6; j++) {
            float4 f4 = spinorPacked[j*L + i];
            res[i*24 + 4*j + 0] = f4.x;
            res[i*24 + 4*j + 1] = f4.y;
            res[i*24 + 4*j + 2] = f4.z;
            res[i*24 + 4*j + 3] = f4.w;
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

void su3_mul(float *res, float *mat, float *vec) {
    for (int n = 0; n < 3; n++) {
        dot(&res[n*(2)], &mat[n*(3*2)], vec);
    }
}

void multiplyGaugeBySpinor(float *res, int dir, float *gauge, float *spinor) {
    for (int s = 0; s < 4; s++) {
        su3_mul(&res[s*(3*2)], gauge, &spinor[s*(3*2)]);
    }
}

void computeGold(float* res, float* gauge, float *spinor) {
    for (int i = 0; i < L; i++) {
        multiplyGaugeBySpinor(&res[i*(4*3*2)], 0, &gauge[i*(4*3*3*2)], &spinor[i*(4*3*2)]);
    }
}


void printVector(float *v) {
    printf("{(%f %f) (%f %f) (%f %f)}\n", v[0], v[1], v[2], v[3], v[4], v[5]);
}


void printSpinorField(float *spinor) {
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[s*(3*2)]);
    }
    printf("...\n");
    for (int s = 0; s < 4; s++) {
        printVector(&spinor[(L-1)*(4*3*2) + s*(3*2)]);
    }
    printf("\n");    
}
