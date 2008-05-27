#include "qcd.h"

extern "C" 
void computeGold(float* reference, float *vecs, float* mats, const int len);

/*
void
dot(float* a, float* b, float& re, float &im) {
    re = im = 0;
    for (int i = 0; i < dim; i++) {
        float a_re = a[2*i+0];
        float a_im = a[2*i+1];
        float b_re = b[2*i+0];
        float b_im = b[2*i+1];
        re += a_re * b_re - a_im * b_im;
        im += a_re * b_im + a_im * b_re;
    }
}

void
computeGold(float* reference, float *vecs, float* mats, const int len) {
    // loop over vectors
    for (int i = 0; i < len; i++) {
        // loop over elements of vector
        for (int j = 0; j < dim; j++) {
            float re = 0, im = 0;
            dot(&mats[i*matsize + j*vecsize], &vecs[i*vecsize], re, im);
            reference[i*vecsize + 2*j + 0] = re;
            reference[i*vecsize + 2*j + 1] = im;
        }
    }
}
*/
