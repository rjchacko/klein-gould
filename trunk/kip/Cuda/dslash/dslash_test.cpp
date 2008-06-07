#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cutil.h>
#include "qcd.h"


void printSpinorHalfField(float *spinor) {
    printSpinor(&spinor[0*(4*3*2)]);
    printf("...\n");
    printSpinor(&spinor[(Nh-1)*(4*3*2)]);
    printf("\n");    
}


void dslashTest(int argc, char **argv) {
    CUT_DEVICE_INIT(argc, argv);
    printCudaDslashInfo();
    
    unsigned int timer = 0;
    cutCreateTimer(&timer);
    
    // construct input fields
    float *gauge[4];
    for (int dir = 0; dir < 4; dir++) {
        gauge[dir] = (float*)malloc(N*gaugeSiteSize*sizeof(float));
    }
    float *spinor = (float*)malloc(N*4*3*2*sizeof(float));
    float *spinorRef = (float*)malloc(N*4*3*2*sizeof(float));
    
    float *spinorEven = spinor;
    float *spinorOdd = spinor + Nh*spinorSiteSize;
    
    printf("Randomizing fields\n");
    constructUnitGaugeField(gauge);
    constructSpinorField(spinor);

    // copy inputs from host to device
    CudaFullGauge cudaGauge = loadGaugeField(gauge);
    CudaFullSpinor cudaSpinor = loadSpinorField(spinor);

    int ODD_BIT = 1;
    int DAGGER_BIT = 1;
    
    // execute kernel
    printf("Beginning kernel execution\n");
    cutStartTimer(timer);
    const int LOOPS = 100;
    for (int i = 0; i < LOOPS; i++) {
        dslashCuda(cudaSpinor.odd, cudaGauge, cudaSpinor.even, ODD_BIT, DAGGER_BIT);
    }
    cutStopTimer(timer);

    // print timing information
    float millisecs = cutGetTimerValue(timer)/LOOPS;
    float secs = millisecs / 1000.;
    printf("Elapsed time %fms\n", millisecs);
    printf("GFLOPS = %f\n", 1e-9*1320*Nh/secs);
    printf("GiB/s = %f\n\n", Nh*(8*7+4)*3*2*sizeof(float)/(secs*(1<<30)));

    // compare to dslash reference implementation
    retrieveParitySpinor(spinorOdd, cudaSpinor.odd);
    dslashReference(spinorRef, gauge, spinorEven, ODD_BIT, DAGGER_BIT);
    
    printf("Reference:\n");
    printSpinorHalfField(spinorRef);
    
    printf("\nCUDA:\n");
    printSpinorHalfField(spinorOdd);
    
    CUTBoolean res = cutComparefe(spinorOdd, spinorRef, Nh*4*3*2, 1e-4);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    
    // release memory
    for (int dir = 0; dir < 4; dir++) {
        free(gauge[dir]);
    }
    free(spinor);
    free(spinorRef);
    freeGaugeField(cudaGauge);
    freeSpinorField(cudaSpinor);
    cutDeleteTimer(timer);
}
