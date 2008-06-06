
#include <stdlib.h>
#include <stdio.h>
#include "qcd.h"

int main() {
    // construct input fields
    float *gaugeEven[4], *gaugeOdd[4];
    for (int dir = 0; dir < 4; dir++) {
        gaugeEven[dir] = (float*)malloc(Nh*3*3*2*sizeof(float));
        gaugeOdd[dir]  = (float*)malloc(Nh*3*3*2*sizeof(float));
    }
    float *spinorEvenIn = (float*)malloc(Nh*4*3*2*sizeof(float));
    float *spinorOddIn  = (float*)malloc(Nh*4*3*2*sizeof(float));
    float *spinorEvenOut = (float*)malloc(Nh*4*3*2*sizeof(float));
    float *spinorOddOut  = (float*)malloc(Nh*4*3*2*sizeof(float));
    
    constructUnitGaugeField(gaugeEven, gaugeOdd);
    
    int i0 = 0;
    int s0 = 0;
    int c0 = 0;
    constructPointSpinorField(spinorEvenIn, spinorOddIn, i0, s0, c0);

    // full dslash operator; final bit indicates "dagger" modifier
    dslashReference(spinorEvenOut, gaugeEven, gaugeOdd, spinorOddIn, 0, 1);
    dslashReference(spinorOddOut, gaugeEven, gaugeOdd, spinorEvenIn, 1, 1);

    // lastly apply the mass term
    float kappa = 0.125;
    float kappa2 = kappa*kappa;
    // apxy(spinorIn, -kappa, spinorOut);
    
    printf("Spinor (t=0,z=0,y=0,x=1):\n");
    printSpinorElement(spinorEvenOut, spinorOddOut, 1);
    printf("Spinor (t=0,z=0,y=1,x=0):\n");
    printSpinorElement(spinorEvenOut, spinorOddOut, L1);
    
    // release memory
    for (int dir = 0; dir < 4; dir++) {
        free(gaugeEven[dir]);
        free(gaugeOdd[dir]);
    }
    free(spinorEvenIn);
    free(spinorOddIn);
    free(spinorEvenOut);
    free(spinorOddOut);
}
