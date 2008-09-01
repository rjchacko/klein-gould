
struct Rand48 {
    // strided iteration constants (48-bit, distributed on 2x 24-bit)
    uint2 A, C;
    // CUDA array -- random numbers for all threads
    uint2 *state;
    // random number for a single thread (used by CUDA device functions only)
    uint2 state0;
    
    // magic constants for rand48
    static const unsigned long long a = 0x5DEECE66DLL, c = 0xB;
    
    void init(int nThreads, int seed) {
        uint2* seeds = new uint2[ nThreads ];
        
        cudaMalloc((void**) &state, sizeof(uint2)*nThreads);
        
        // calculate strided iteration constants
        unsigned long long A, C;
        A = 1LL; C = 0LL;
        for (unsigned int i = 0; i < nThreads; ++i) {
            C += A*c;
            A *= a;
        }
        this->A.x = A & 0xFFFFFFLL;
        this->A.y = (A >> 24) & 0xFFFFFFLL;
        this->C.x = C & 0xFFFFFFLL;
        this->C.y = (C >> 24) & 0xFFFFFFLL;
        
        // prepare first nThreads random numbers from seed
        unsigned long long x = (((unsigned long long)seed) << 16) | 0x330E;
        for (unsigned int i = 0; i < nThreads; ++i) {
            x = a*x + c;
            seeds[i].x = x & 0xFFFFFFLL;
            seeds[i].y = (x >> 24) & 0xFFFFFFLL;
        }
        
        cudaMemcpy(state, seeds, sizeof(uint2)*nThreads, cudaMemcpyHostToDevice);
        
        delete[] seeds;
    }
    
    void destroy() {
        cudaFree((void*) state);
    }
};

__device__ inline void rand48_loadState(Rand48 &r) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    r.state0 = r.state[i];
}

__device__ inline void rand48_storeState(Rand48 &r) {
    int i = threadIdx.x + blockIdx.x*blockDim.x;
    r.state[i] = r.state0;
}

__device__ inline void rand48_iterate(Rand48 &r) {
    // state0 is 2x 24bit to handle overflows optimally, i.e.
    // in one operation.
    
    // the multiplication commands however give the low and hi 32 bit,
    // which have to be converted as follows:
    // 48bit in bytes = ABCD EF (space marks 32bit boundary)
    // R0             = ABC
    // R1             =    D EF
    
    unsigned int R0, R1;
    
    // low 24-bit multiplication
    const unsigned int lo00 = __umul24(r.state0.x, r.A.x);
    const unsigned int hi00 = __umulhi(r.state0.x, r.A.x);
    
    // 24bit distribution of 32bit multiplication results
    R0 = (lo00 & 0xFFFFFF);
    R1 = (lo00 >> 24) | (hi00 << 8);
    
    R0 += r.C.x; R1 += r.C.y;
    
    // transfer overflows
    R1 += (R0 >> 24);
    R0 &= 0xFFFFFF;
    
    // cross-terms, low/hi 24-bit multiplication
    R1 += __umul24(r.state0.y, r.A.x);
    R1 += __umul24(r.state0.x, r.A.y);
    
    R1 &= 0xFFFFFF;
    
    r.state0 = make_uint2(R0, R1);
}

__device__ inline int rand48_nextInt(Rand48 &r) {
    // get upper 31 (!) bits of the 2x 24bits
    int res = ( r.state0.x >> 17 ) | ( r.state0.y << 7 );
    rand48_iterate(r);
    return res;
}

// returns a float in the range (0, 1]
__device__ inline int rand48_nextFloat(Rand48 &r) {
    return ((float)rand48_nextInt(r)+1.0f) / (1<<31);
}
