#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "ising.h"


// ----------------------------------------------------------------------------
// Cuda implementation
//

#define GRID_DIM 128
#define BLOCK_DIM 128

#define CUDA_INCLUDE
#include "bits.cpp" // include CUDA device functions
#include "rand48.cu" // random numbers

typedef struct {
    int len, dim, nblocks;
    unsigned int *blocks;
    float h, T; // external field and temperature
    int parityTarget;
} IsingCudaParams;

// ----------------------------------------------------------------------------
// Cuda lookup table (JEE)
// The indexing is as follows:
//      for m*s, look at table index (m + 2*dim)/2


// This variable stays in CUDA memory until the termination of the program
// Note that we must choose a sufficiently large array. Dynamic memory
// allocation does not work for type __constant__.
#define TABLESIZE 2*(2*7+1)
__constant__ unsigned int table[TABLESIZE];

void IsingCuda::buildLookupTable ()
{
    // Note tables are of size 2*len in order to take care of both spin
    // cases
    int len = 2*dim + 1;

    // Build table on host then copy it to device
    unsigned int * localTable = 
        (unsigned int *) malloc (2*len*sizeof (unsigned int));
    for (int i=0; i<=1; ++i)
    {
        int s = (i==0 ? -1 : 1);
        for (int m=-2*dim; m<=2*dim; m+=2)
        {
            int index = (m + 2*dim)/2 + i*len;
            double dE = 2*s*(m + h);
            if (dE < 0)
                localTable[index] = KIP_RAND_MAX;
            else
                localTable[index] = (unsigned int) round(exp( -dE / T )*KIP_RAND_MAX);
        }
    }

    cudaMemcpyToSymbol 
    (table,localTable, 2*len*sizeof (unsigned int));

    free (localTable);
}

void IsingCuda::flipH ()
{
    h = -h;
    buildLookupTable ();
}

__device__ inline int shouldFlipSpin_table
(IsingCudaParams p, Rand48 &rng, int s, int m) 
{
    int len = 2*p.dim + 1;
    //int spinIndexOffset = ( s==-1 ? 0 : 1 )*len;
    int spinIndexOffset = ((s+1)>>1)*len;
    //int index = (m + 2*p.dim)/2 + spinIndexOffset;
    int index = ((m + 2*p.dim)>>1) + spinIndexOffset;

    return rand48_nextInt (rng) < table[index];
}

// ----------------------------------------------------------------------------
// Cuda magnetization calculation

#define REDUCE_THREADS 128
#define REDUCE_MAX_BLOCKS 128

__global__ void isingCuda_magnetizationKernel (unsigned int *g_idata, unsigned int *g_odata, unsigned int n) {
    unsigned int tid = threadIdx.x;
    unsigned int i = blockIdx.x*REDUCE_THREADS + threadIdx.x;
    unsigned int gridSize = gridDim.x*REDUCE_THREADS;
    
    unsigned int acc = 0;
    while (i < n) {
        acc += (unsigned int) bitCount(g_idata[i]);
        i += gridSize;
    }

    extern __shared__ unsigned int s[];
    s[tid] = acc;
    __syncthreads();
    
    // do reduction in shared mem
    if (REDUCE_THREADS >= 256) { if (tid < 128) { s[tid] += s[128+tid]; } __syncthreads(); }
    if (REDUCE_THREADS >= 128) { if (tid <  64) { s[tid] += s[ 64+tid]; } __syncthreads(); }
    if (tid < 32) {
        s[tid] += s[32+tid];
        s[tid] += s[16+tid];
        s[tid] += s[ 8+tid];
        s[tid] += s[ 4+tid];
        s[tid] += s[ 2+tid];
        s[tid] += s[ 1+tid];
    }
    if (tid == 0) g_odata[blockIdx.x] = s[tid];
}

int divideCeil(int x, int y) {
    return (x + y - 1) / y;
}

double isingCuda_bitCount(unsigned int *d_idata, int n) {
    // allocate arrays on device and host to store one float for each block
    int blocks = min(REDUCE_MAX_BLOCKS, divideCeil(n, REDUCE_THREADS));
    unsigned int h_odata[REDUCE_MAX_BLOCKS];
    unsigned int *d_odata;
    cudaMalloc((void**) &d_odata, blocks*sizeof(unsigned int));
    
    // partial reduction; each block generates one number
    dim3 dimBlock(REDUCE_THREADS, 1, 1);
    dim3 dimGrid(blocks, 1, 1);
    int smemSize = REDUCE_THREADS * sizeof(unsigned int);
    isingCuda_magnetizationKernel<<< blocks, REDUCE_THREADS, smemSize >>>(d_idata, d_odata, n);
    
    // copy result from device to host, and perform final reduction on CPU
    cudaMemcpy(h_odata, d_odata, blocks*sizeof(unsigned int), cudaMemcpyDeviceToHost);
    double gpu_result = 0;
    for (int i = 0; i < blocks; i++)  {
        gpu_result += h_odata[i];
    }
    cudaFree(d_odata);    
    return gpu_result;
}

// ----------------------------------------------------------------------------
// Cuda Energy calculation JEE

__device__ inline void isingCuda_localInternalEnergy
(IsingCudaParams p, int ip, int & internal_energy) {

    int parity = 0;
    Bits128 acc = {0, 0, 0, 0};
    int lenp_d = 1;
    Bits128 n1 = bitsExpand(p.blocks[ip]);
    
    for (int d = 0; d < p.dim; d++) {
        int lenp = (d < 5) ? p.len>>2 : p.len;
        int xp = (ip / lenp_d) % lenp;
        parity += (d < 5) ? 0 : xp;
        
        int dx2 = (xp+1+lenp)%lenp - xp;
        int dx0 = (xp-1+lenp)%lenp - xp;
        Bits128 n2 = bitsExpand(p.blocks[ip+dx2*lenp_d]);
        Bits128 n0 = bitsExpand(p.blocks[ip+dx0*lenp_d]);
        
        if (d < 5) {
            int shift = 4 << d; // 4, 8, 16, 32, 64
            acc = bitsAdd(acc, bitsMaskShiftL(bitsAdd(n1,n2), shift));
            acc = bitsAdd(acc, bitsMaskShiftR(bitsAdd(n1,n0), shift));
        }
        else {
            acc = bitsAdd(bitsAdd(n0,n2), acc);
        }
        
        lenp_d *= lenp;
    }

    int deltaMax = p.dim < 5 ? (1 << p.dim) : 32;
    int cube = p.blocks[ip];
    for (int delta = 0; delta < deltaMax; delta++) {
        // Make sure we only test even parity, to avoid double bond
        // counting. This check is not done in the primary kernel.
        if (((parity + bitCount(delta)) & 1) == 0) {
            // m = total magnetization of neighbors; in range [-2 dim, +2 dim]
            int m = 2*(bitsPick4(acc, delta) - p.dim);
            // s = spin at this site (+/- 1)
            int s = 2*((cube >> delta) & 1) - 1;
            internal_energy = - m * s;
        }
    }
}

__global__ void isingCuda_internalEnergyKernel
(IsingCudaParams p, int * odata) {

    unsigned int tid = threadIdx.x;
    unsigned int ip = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;

    int acc = 0;
    int ie; // internal energy
    while (ip < p.nblocks) {
        isingCuda_localInternalEnergy (p, ip, ie);
        acc += ie;
        ip += gridSize;
    }

    extern __shared__ int t[];
    t[tid] = acc;
    __syncthreads();

    // do thread reduction
    
    if (REDUCE_THREADS >= 256) { if (tid < 128) { t[tid] += t[128+tid]; } __syncthreads(); }
    if (REDUCE_THREADS >= 128) { if (tid <  64) { t[tid] += t[ 64+tid]; } __syncthreads(); }
    if (tid < 32) {
        t[tid] += t[32+tid];
        t[tid] += t[16+tid];
        t[tid] += t[ 8+tid];
        t[tid] += t[ 4+tid];
        t[tid] += t[ 2+tid];
        t[tid] += t[ 1+tid];
    }
    if (tid == 0) odata[blockIdx.x] = t[tid];
}


// ----------------------------------------------------------------------------
// Cuda update implementation

// JEE should flip spin from table function is above
__device__ inline int shouldFlipSpin(IsingCudaParams p, Rand48 &rng, int s, int m) {
    float dE = 2*s*(m + p.h);
#ifdef DETERMINISTIC
    float r = 0.1;
#else
    float r = (float)rand48_nextInt(rng) / (unsigned int)(1<<31);
#endif
    return __expf (- dE / p.T) > r;
}

__device__ inline void isingCuda_updateSite(IsingCudaParams p, Rand48 &rng, int ip) {

    int parity = 0;
    Bits128 acc = {0, 0, 0, 0};
    int lenp_d = 1;
    Bits128 n1 = bitsExpand(p.blocks[ip]);
    
    for (int d = 0; d < p.dim; d++) {
        int lenp = (d < 5) ? p.len>>2 : p.len;
        int xp = (ip / lenp_d) % lenp;
        parity += (d < 5) ? 0 : xp;
        
        int dx2 = (xp+1+lenp)%lenp - xp;
        int dx0 = (xp-1+lenp)%lenp - xp;
        Bits128 n2 = bitsExpand(p.blocks[ip+dx2*lenp_d]);
        Bits128 n0 = bitsExpand(p.blocks[ip+dx0*lenp_d]);
        
        if (d < 5) {
            int shift = 4 << d; // 4, 8, 16, 32, 64
            acc = bitsAdd(acc, bitsMaskShiftL(bitsAdd(n1,n2), shift));
            acc = bitsAdd(acc, bitsMaskShiftR(bitsAdd(n1,n0), shift));
        }
        else {
            acc = bitsAdd(bitsAdd(n0,n2), acc);
        }
        
        lenp_d *= lenp;
    }
    
    int deltaMax = p.dim < 5 ? (1 << p.dim) : 32;
    int cube = p.blocks[ip];
    for (int delta = 0; delta < deltaMax; delta++) {
        if (((parity + bitCount(delta)) & 1) == p.parityTarget) {
            // m = total magnetization of neighbors; in range [-2 dim, +2 dim]
            int m = 2*(bitsPick4(acc, delta) - p.dim);
            // s = spin at this site (+/- 1)
            int s = 2*((cube >> delta) & 1) - 1;
            //if (shouldFlipSpin(p, rng, s, m)) {
            if (shouldFlipSpin_table(p, rng, s, m)) {
                cube ^= (1 << delta);
            }
        }
    }
    p.blocks[ip] = cube;
}

__global__ void isingCuda_update(IsingCudaParams p, Rand48 rng) {
    rand48_loadState(rng);

    unsigned int ip = blockIdx.x*(blockDim.x) + threadIdx.x;
    unsigned int gridSize = gridDim.x*blockDim.x;
    while (ip < p.nblocks) {
        isingCuda_updateSite(p, rng, ip);
        ip += gridSize;
    }
    
    rand48_storeState(rng);
}


// ----------------------------------------------------------------------------
// Ising class interface
// JEE
// functions for lookup table above

IsingCuda::IsingCuda(int len, int dim, float h, float T) : Ising(len, dim, h, T) {
    assert(len % 2 == 0);
    assert(dim <= 7);
    len = len;
    dim = dim;
    n = (int)powl(len, dim);
    nblocks = n >> min(5,dim);
    
    int nbytes = nblocks*sizeof(unsigned int);
    blocks = (unsigned int *)malloc(nbytes);
    cudaMalloc((void**)&d_blocks, nbytes);
    
    for (int i = 0; i < nblocks; i++) {
        blocks[i] = 0;
    }
    transferHostToDevice();
    
    rng = new Rand48();
    rng->init(GRID_DIM*BLOCK_DIM, 0); // initialize random numbers

    // JEE Initialize lookup table
    buildLookupTable ();
}

IsingCuda::~IsingCuda() {
    free(blocks);
    cudaFree(d_blocks);
    rng->destroy();
    delete rng;

    // Cleanup lookup table
    cudaFree (table);
}


void IsingCuda::completeNeighborSum(int *sum) {
    assert(0==1);
}

void IsingCuda::update(int parityTarget) {
    IsingCudaParams p;
    p.len = len;
    p.dim = dim;
    p.nblocks = nblocks;
    p.blocks = d_blocks;
    p.h = h;
    p.T = T;
    p.parityTarget = parityTarget;

    int sharedBytes = 0;
    
    isingCuda_update <<<GRID_DIM, BLOCK_DIM, sharedBytes>>> (p, *rng);
}

double IsingCuda::magnetization() {
    return 2.0*isingCuda_bitCount(d_blocks, nblocks) - n;
}

double IsingCuda::energy () 
{
    IsingCudaParams p;
    p.len = len;
    p.dim = dim;
    p.nblocks = nblocks;
    p.blocks = d_blocks;
    p.h = h;
    p.T = T;
    p.parityTarget = 0; // Unnecessary

    int h_odata [BLOCK_DIM];
    int *d_odata;
    cudaMalloc((void**) &d_odata, BLOCK_DIM*sizeof(int));

    int sharedBytes = REDUCE_THREADS * sizeof(unsigned int);
    
    isingCuda_internalEnergyKernel 
        <<<GRID_DIM, BLOCK_DIM, sharedBytes>>> (p, d_odata);

    cudaMemcpy 
        (h_odata, d_odata, BLOCK_DIM*sizeof (int), cudaMemcpyDeviceToHost);
    cudaFree (d_odata); // dont need these any more
    double ie = 0;
    for (int i=0; i<BLOCK_DIM; ++i)
        ie += (double) h_odata[i];

    double m = magnetization ();

    return ie - m * h;
}

void IsingCuda::transferHostToDevice() {
    cudaMemcpy(d_blocks, blocks, nblocks*sizeof(unsigned int), cudaMemcpyHostToDevice);
}

void IsingCuda::transferDeviceToHost() {
    cudaMemcpy(blocks, d_blocks, nblocks*sizeof(unsigned int), cudaMemcpyDeviceToHost);
}


// given index 'i' into the full lattice, return compressed index 'ip'
// and bit index 'delta'.
void IsingCuda::index(int i, int *ip, int *delta) {
    int len_d = 1;
    int lenp_d = 1;
    *ip = 0;
    *delta = 0;
    
    for (int d = 0; d < dim; d++) {
        int x = (i / len_d) % len;
        int xp = (d < 5) ? x/2 : x;
        int del = (d < 5) ? x%2 : 0;
        
        *delta += (del << d) ;
        *ip += xp*lenp_d;
        
        int lenp = (d < 5) ? len/2 : len;
        len_d *= len;
        lenp_d *= lenp;
    }
    
    assert(*ip < nblocks);
    assert(*delta < 32);
}

void IsingCuda::set(int i, int s) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    int mask = ~(1 << delta);
    blocks[ip] = (blocks[ip] & mask) | (s << delta);
}

int IsingCuda::get(int i) {
    int ip, delta;
    index(i, &ip, &delta);
    assert(ip < nblocks);
    return (blocks[ip]>>delta) & 1;
}

// ----------------------------------------------------------------------------
// Cuda utility methods
//

void initCuda(int argc, char *argv[]) {
    int deviceCount;
    cudaGetDeviceCount(&deviceCount);
    if (deviceCount == 0) {
        fprintf(stderr, "No devices supporting CUDA.\n");
        exit(EXIT_FAILURE);
    }
    int dev = deviceCount - 1;
    if (argc > 1) {
        sscanf(argv[1], "%d", &dev);
    }
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    if (deviceProp.major < 1) {
        fprintf(stderr, "Device %d does not support CUDA.\n", dev);
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "Using device %d: %s\n", dev, deviceProp.name);
    cudaSetDevice(dev);
}
