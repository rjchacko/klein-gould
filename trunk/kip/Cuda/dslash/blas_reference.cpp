#include "qcd.h"


// sets all elements of the destination vector to zero
void zero(float* a, int cnt) {
    for (int i = 0; i < cnt; i++)
        a[i] = 0;
}

// copy one spinor to the other
void copy(float* a, float *b, int len) {
  for (int i = 0; i < len; i++) a[i] = b[i];
}

// performs the operation x[i] *= a
void ax(float a, float *x, int len) {
    for (int i=0; i<len; i++) x[i] *= a;
}

// performs the operation y[i] = a*x[i] + b*y[i]
void axpby(float a, float *x, float b, float *y, int len) {
    for (int i=0; i<len; i++) y[i] = a*x[i] + b*y[i];
}

// performs the operation y[i] = a*x[i] + y[i]
void axpy(float a, float *x, float *y, int len) {
    for (int i=0; i<len; i++) y[i] += a*x[i];
}

// performs the operation y[i] = x[i] + a*y[i]
void xpay(float *x, float a, float *y, int len) {
    for (int i=0; i<len; i++) y[i] = x[i] + a*y[i];
}

// performs the operation y[i] -= x[i] (minus x plus y)
void mxpy(float *x, float *y, int len) {
    for (int i=0; i<len; i++) y[i] -= x[i];
}


// returns the square of the L2 norm of the vector
float norm(float *v, int len) {

  float sum=0.0;
  for (int i=0; i<len; i++) {
    sum += v[i]*v[i];
  }

  return sum;
}

// returns the real part of the dot product of 2 complex valued vectors
float reDotProduct(float *v1, float *v2, int len) {

  float dot=0.0;
  for (int i=0; i<len; i++) {
    dot += v1[i]*v2[i];
  }

  return dot;
}

// returns the imaginary part of the dot product of 2 complex valued vectors
float imDotProduct(float *v1, float *v2, int len) {

  float dot=0.0;
  for (int i=0; i<len; i+=2) {
    dot += v1[i]*v2[i+1] - v1[i+1]*v2[i];
  }

  return dot;
}

// returns the square of the L2 norm of the vector
double normD(float *v, int len) {

  double sum=0.0;
  for (int i=0; i<len; i++) {
    sum += v[i]*v[i];
  }

  return sum;
}

// returns the real part of the dot product of 2 complex valued vectors
double reDotProductD(float *v1, float *v2, int len) {

  double dot=0.0;
  for (int i=0; i<len; i++) {
    dot += v1[i]*v2[i];
  }

  return dot;
}

// returns the imaginary part of the dot product of 2 complex valued vectors
double imDotProductD(float *v1, float *v2, int len) {

  double dot=0.0;
  for (int i=0; i<len; i+=2) {
    dot += v1[i]*v2[i+1] - v1[i+1]*v2[i];
  }

  return dot;
}
