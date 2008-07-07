#include <stdlib.h>
#include <stdio.h>
#include "qcd.h"


void cgTest() {
  float mass = 0.01;
  float kappa = 1.0 / (2.0*(4 + mass));

  float *gauge[4];
  for (int dir = 0; dir < 4; dir++) {
    gauge[dir] = (float*)malloc(N*gaugeSiteSize*sizeof(float));
  }
  //constructGaugeField(gauge);
  constructUnitGaugeField(gauge);
    
  float *spinorIn = (float*)malloc(N*spinorSiteSize*sizeof(float));
  float *spinorOut = (float*)malloc(N*spinorSiteSize*sizeof(float));

#ifdef EVEN_ODD
  float *source = (float *)malloc(Nh*spinorSiteSize*sizeof(float));
  float *tmp = (float *)malloc(Nh*spinorSiteSize*sizeof(float));
#else
  float *source = (float *)malloc(N*spinorSiteSize*sizeof(float));
#endif

  int i0 = 0;
  int s0 = 0;
  int c0 = 0;
  constructPointSpinorField(spinorIn, i0, s0, c0);
  //constructSpinorField(spinorIn);

  // Prepare the source term
  ax(2*kappa, spinorIn, N*spinorSiteSize);

  // see output element
  // Mat(source, gauge, spinorIn, kappa);
  // printSpinorElement(source, 0);

#ifdef EVEN_ODD
  float *spinorInOdd = spinorIn + Nh*spinorSiteSize;
  dslashReference(tmp, gauge, spinorInOdd, 0, 0);
  xpay(spinorIn, kappa, tmp, Nh*spinorSiteSize);
  MatPCDag(source, gauge, tmp, kappa);
#else
  MatDag(source, gauge, spinorIn, kappa);
#endif

  cgCuda(spinorOut, gauge, source, kappa, 1e-7);
//  cg_reference(spinorOut, gauge, source, kappa, 1e-7);

  // Reconstruct the full inverse
#ifdef EVEN_ODD
  float *spinorOutOdd = spinorOut + Nh*spinorSiteSize;
  dslashReference(spinorOutOdd, gauge, spinorOut, 1, 0);
  xpay(spinorInOdd, kappa, spinorOutOdd, Nh*spinorSiteSize);
#endif

  printf("Result norm = %e\n", norm(spinorOut, N*spinorSiteSize));

  // release memory
  for (int dir = 0; dir < 4; dir++) free(gauge[dir]);
  free(source);
#ifdef EVEN_ODD
  free(tmp);
#endif
  free(spinorIn);
  free(spinorOut);
}
