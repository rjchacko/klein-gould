#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "qcd.h"


void cgReference(float *x, float **gauge, float *b, float kappa, float tol) {
#ifdef EVEN_ODD
  int len = Nh*spinorSiteSize;
#else
  int len = N*spinorSiteSize;
#endif

  float *r = (float*)malloc(len*spinorSiteSize*sizeof(float));
  float *p = (float*)malloc(len*spinorSiteSize*sizeof(float));
  float *Ap = (float*)malloc(len*spinorSiteSize*sizeof(float));

  float b2 = norm(b, len);
  float r2 = b2;
  float r2_old;
  float stop = r2*tol*tol; // stopping condition of solver

  float alpha, beta, pAp;

  copy(r, b, len);

  copy(p, r, len);
  zero(x, len);

  int k=0;
  printf("%d iterations, r2 = %e\n", k, r2);
  while (r2 > stop) {

#ifdef EVEN_ODD
    MatPCDagMatPC(Ap, gauge, p, kappa);
#else
    MatDagMat(Ap, gauge, p, kappa);
#endif

    pAp = reDotProduct(p, Ap, len);
    
    alpha = r2 / pAp;

    axpy(alpha,p,x,len);
    axpy(-alpha,Ap,r,len);

    r2_old = r2;
    r2 = norm(r, len);

    beta = r2 / r2_old;

    xpay(r, beta, p, len);

    k++;
    printf("%d iterations, r2 = %e %e\n", k, r2, norm(x, len));
  }

  // Calculate the true residual
#ifdef EVEN_ODD
  MatPCDagMatPC(Ap, gauge, x, kappa);
#else
  MatDagMat(Ap, gauge, x, kappa);
#endif

  copy(r, b, len);
  mxpy(Ap, r, len);
  double true_res = norm(r, len);


  printf("Converged after %d iterations, r2 = %e, true_res^2 = %e\n", 
	 k, r2, true_res);

  free(Ap);
  free(p);
  free(r);
  return;
}
