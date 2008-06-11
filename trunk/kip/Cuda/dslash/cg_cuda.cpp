#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qcd.h"


void cgCuda(float *h_x, float **h_gauge, float *h_b, float kappa, float tol) {
#ifndef EVEN_ODD
    printf("Cuda CG only does preconditioned inverse.\n");
    return;
#endif
    
    float spinorGiB = (float)Nh*spinorSiteSize*sizeof(float) / (1 << 30);
    float gaugeGiB = (float)4*N*packedGaugeSiteSize*sizeof(float) / (1 << 30);
    printf("Cuda Space Required. Spinor:%f + Gauge:%f GiB\n", 7*spinorGiB, gaugeGiB);
    
    int len = Nh*spinorSiteSize;
    
    CudaFullGauge gauge = loadGaugeField(h_gauge);
    CudaPSpinor b = loadParitySpinor(h_b);
    
    CudaPSpinor x = allocateParitySpinor();
    CudaPSpinor r = allocateParitySpinor();
    CudaPSpinor p = allocateParitySpinor();
    CudaPSpinor Ap = allocateParitySpinor();
    CudaPSpinor tmp1 = allocateParitySpinor();
    CudaPSpinor tmp2 = allocateParitySpinor();
    
    float b2 = normCuda((float *)b, len);
    float r2 = b2;
    float r2_old;
    float stop = r2*tol*tol; // stopping condition of solver
    
    float alpha, beta, pAp;
    
    
    copyCuda((float *)r, (float *)b, len);
    
    copyCuda((float *)p, (float *)r, len);
    zeroCuda((float *)x, len);
    
    int k=0;
    // printf("%d iterations, r2 = %e\n", k, r2);
    stopwatchStart();
    while (r2 > stop) {
        MatPCDagMatPCCuda(Ap, gauge, p, kappa, tmp1, tmp2);
        
        pAp = reDotProductCuda((float *)p, (float *)Ap, len);
        
        alpha = r2 / pAp;        
        r2_old = r2;
        
        axpyCuda(-alpha,(float *)Ap,(float *)r,len);
        r2 = normCuda((float *)r, len);
        
        beta = r2 / r2_old;
        
        axpyZpbxCuda(alpha, (float *)p, (float *)x, (float *)r, beta, len);
        
        k++;
        // printf("%d iterations, r2 = %e %e\n", k, r2, normCuda((float *)x, len));
    }
    float gflops = (1.0e-9*Nh)*(4*1320 + 14*spinorSiteSize);
    printf("%f gflops\n", k*gflops / stopwatchReadSeconds());
    
    // Calculate the true residual
    MatPCDagMatPCCuda(Ap, gauge, x, kappa, tmp1, tmp2);
    
    copyCuda((float *)r, (float *)b, len);
    mxpyCuda((float *)Ap, (float *)r, len);
    double true_res = normCuda((float *)r, len);
    
    
    printf("Converged after %d iterations, r2 = %e, true_res^2 = %e\n", k, r2, true_res);

    retrieveParitySpinor(h_x, x);
    
    freeGaugeField(gauge);
    freeParitySpinor(x);
    freeParitySpinor(b);
    freeParitySpinor(r);
    freeParitySpinor(p);
    freeParitySpinor(Ap);
    freeParitySpinor(tmp1);
    freeParitySpinor(tmp2);
    
    return;
}
