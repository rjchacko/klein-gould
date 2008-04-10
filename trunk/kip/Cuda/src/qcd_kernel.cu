
#ifndef _QCD_KERNEL_H_
#define _QCD_KERNEL_H_

#include <stdio.h>
#include <qcd.h>

// #define SDATA( index)      CUT_BANK_CHECKER(sdata, index)

extern __shared__ float sdata[];


__device__ void
dot(float* vec1, float *vec2, float* output) {
    float out_re = 0.0f, out_im = 0.0f;
    for (int i = 0; i < dim; i++) {
        float v1_re = vec1[2*i + 0];
        float v1_im = vec1[2*i + 1];
        float v2_re = vec2[2*i + 0];
        float v2_im = vec2[2*i + 1];
        out_re += v1_re*v2_re - v1_im*v2_im;
        out_im += v1_re*v2_im + v1_im*v2_re;
    }
    output[0] = out_re;
    output[1] = out_im;
}


////////////////////////////////////////////////////////////////////////////////
//! @param g_vdata  input vectors
//! @param g_mdata  input matrices
//! @param g_odata  output vectors
////////////////////////////////////////////////////////////////////////////////
__global__ void
testKernel(float* g_vecs, float* g_mats, float* g_out) {
    const int tid = threadIdx.x;        // thread index
    const int bid = blockIdx.x;         // block index
    const int nts = blockDim.x;         // number of threads in one block
    const int matsPerBlock = nts / dim;
    
    float *vecs = &sdata[0*nts];
    float *out  = &sdata[2*nts];
    float *mats = &sdata[4*nts];
    
    // offset global array bases by the appropriate block index.
    g_mats = &g_mats[matsize*matsPerBlock*bid];
    g_vecs = &g_vecs[vecsize*matsPerBlock*bid];
    g_out  = &g_out [vecsize*matsPerBlock*bid];
    
    // load vector and matrix data for this block; loading a vector takes
    // two iterations since each thread maps to a single complex
    // vector element (which has two float components.)
    // n.b.: for efficiency it is crucial to load contiguous blocks of memory.
    for (int i = 0; i < 2; i++)
        vecs[i*nts+tid] = g_vecs[i*nts+tid];
    for (int i = 0; i < 2*dim; i++)
        mats[i*nts+tid] = g_mats[i*nts+tid];
    __syncthreads();
    
    // vector id
    int vid = tid/dim; // particularly costly operation; can eliminate by reshaping block
    int rem = tid - dim*vid;
    
    // perform matrix vector multiplies
    dot(&vecs[vecsize*vid], &mats[matsize*vid+vecsize*rem], &out[vecsize*vid+2*rem]);
    __syncthreads();
    
    // write data to global memory
    for (int i = 0; i < 2; i++)
        g_out[i*nts+tid] = out[i*nts+tid];
}

#endif // #ifndef _TEMPLATE_KERNEL_H_
