
#include <stdlib.h>
#include <stdio.h>
#include "qcd.h"

void floatTimesSpinor(float a, float *x, int len) {
    for (int i=0; i<len; i++) x[i] *= a;
}

void axpy(float a, float *x, float *y, int len) {
    for (int i=0; i<len; i++) y[i] += a*x[i];
}

int main() {
    int gaugeSiteSize = 3*3*2;
    int spinorSiteSize = 4*3*2;
    
    float *gauge[4], *gaugeEven[4], *gaugeOdd[4];
    for (int dir = 0; dir < 4; dir++) {
        gauge[dir] = (float*)malloc(N*gaugeSiteSize*sizeof(float));
        gaugeEven[dir] = gauge[dir];
        gaugeOdd[dir]  = gauge[dir]+(N/2)*gaugeSiteSize;
    }
    
    float *spinorIn = (float*)malloc(N*spinorSiteSize*sizeof(float));
    float *spinorEvenIn = spinorIn;
    float *spinorOddIn  = spinorIn + (N/2)*spinorSiteSize;
    float *spinorOut = (float*)malloc(N*spinorSiteSize*sizeof(float));
    float *spinorEvenOut = spinorOut;
    float *spinorOddOut  = spinorOut + (N/2)*spinorSiteSize;
    
    constructUnitGaugeField(gaugeEven, gaugeOdd);
    
    int i0 = 0;
    int s0 = 0;
    int c0 = 0;
    constructPointSpinorField(spinorEvenIn, spinorOddIn, i0, s0, c0);

    // full dslash operator
    int daggerBit = 0;
    dslashReference(spinorEvenOut, gaugeEven, gaugeOdd, spinorOddIn, 0, daggerBit);
    dslashReference(spinorOddOut, gaugeEven, gaugeOdd, spinorEvenIn, 1, daggerBit);

    // factor of -0.5
    floatTimesSpinor(-0.5, spinorOut, N*spinorSiteSize);
    
    // lastly apply the mass term
    float mass = 0.01;
    float d = 4.0+mass;
    axpy(d, spinorIn, spinorOut, N*spinorSiteSize);
    
    for (int i = 0; i < N; i++) {
        printSpinorElement(spinorEvenOut, spinorOddOut, i);
        printf("\n");
    }
    
    // release memory
    for (int dir = 0; dir < 4; dir++) {
        free(gauge[dir]);
    }
    free(spinorIn);
    free(spinorOut);
}
