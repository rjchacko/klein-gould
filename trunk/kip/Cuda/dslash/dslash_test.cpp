#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include "qcd.h"

#define FULL_WILSON 1


void printSpinorHalfField(float *spinor) {
    printSpinor(&spinor[0*(4*3*2)]);
    printf("...\n");
    printSpinor(&spinor[(Nh-1)*(4*3*2)]);
    printf("\n");    
}


void dslashTest() {
    float spinorGiB = (float)Nh*spinorSiteSize*sizeof(float) / (1 << 30);
    float gaugeGiB = (float)4*N*packedGaugeSiteSize*sizeof(float) / (1 << 30);
    float sharedKB = (float)dslashCudaSharedBytes() / (1 << 10);
    printf("Spinor mem: %f GiB\n", spinorGiB);
    printf(" Gauge mem: %f GiB\n", gaugeGiB);
    printf("Shared mem: %f KB\n", sharedKB);
    
    // construct input fields
    float *gauge[4];
    for (int dir = 0; dir < 4; dir++) {
        gauge[dir] = (float*)malloc(N*gaugeSiteSize*sizeof(float));
    }
    float *spinor = (float*)malloc(N*4*3*2*sizeof(float));
    float *spinorRef = (float*)malloc(N*4*3*2*sizeof(float));
    
    float *spinorEven = spinor;
    float *spinorOdd = spinor + Nh*spinorSiteSize;
    
    printf("Randomizing fields...");
    fflush(stdout);
    constructGaugeField(gauge);
    constructSpinorField(spinor);
    printf("done.\n");
    
    // copy inputs from host to device
    printf("Sending fields to GPU...");
    fflush(stdout);
    CudaFullGauge cudaGauge = loadGaugeField(gauge);
    CudaFullSpinor cudaSpinor = loadSpinorField(spinor);
    CudaPSpinor tmp = allocateParitySpinor();
    printf("done.\n");
    
    float kappa = 1.0;
    int ODD_BIT = 0;
    int DAGGER_BIT = 0;
    
    // execute kernel
    const int LOOPS = 200;
    printf("Executing %d kernel loops...", LOOPS);
    fflush(stdout);
    stopwatchStart();
    for (int i = 0; i < LOOPS; i++) {
    	if (FULL_WILSON)
    		MatPCCuda(cudaSpinor.odd, cudaGauge, cudaSpinor.even, kappa, tmp);
		else
        	dslashCuda(cudaSpinor.odd, cudaGauge, cudaSpinor.even, ODD_BIT, DAGGER_BIT);
    }
    cudaThreadSynchronize();
    double secs = stopwatchReadSeconds() / LOOPS;
    printf("done.\n\n");
    
    // print timing information
    printf("%fms per loop\n", 1000*secs);
    int flops = FULL_WILSON ? 1320*2 + 48 : 1320;
    int floats = FULL_WILSON ? 2*(8*(24+12)+24)+24 : 8*(24+12)+24;
    printf("GFLOPS = %f\n", 1.0e-9*flops*Nh/secs);
    printf("GiB/s = %f\n\n", Nh*floats*sizeof(float)/(secs*(1<<30)));

    // compare to dslash reference implementation
    printf("Comparing to reference implementation...");
    fflush(stdout);
    retrieveParitySpinor(spinorOdd, cudaSpinor.odd);
    if (FULL_WILSON)
	    MatPC(spinorRef, gauge, spinorEven, kappa);
    else
    	dslashReference(spinorRef, gauge, spinorEven, ODD_BIT, DAGGER_BIT);
    printf("done.\n");
    
    printf("Reference:\n");
    printSpinorHalfField(spinorRef);
    
    printf("\nCUDA:\n");
    printSpinorHalfField(spinorOdd);
    
    int res = compareFloats(spinorOdd, spinorRef, Nh*4*3*2, 1e-4);
    printf("Test %s\n", (1 == res) ? "PASSED" : "FAILED");
    
    // release memory
    for (int dir = 0; dir < 4; dir++) {
        free(gauge[dir]);
    }
    free(spinor);
    free(spinorRef);
    freeGaugeField(cudaGauge);
    freeSpinorField(cudaSpinor);
}
