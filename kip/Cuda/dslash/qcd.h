
#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L4 64 // "time" dimension

#define L1h (L1/2) // half of the full "x" dimension, useful for even/odd lattices
#define L (L1h*L2*L3*L4) // total number of even/odd lattice points

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (L/BLOCK_DIM) // there are L threads in total

//---------- dslash_reference.cpp

extern "C" void constructGaugeField(float **resEven, float **resOdd);
extern "C" void constructSpinorField(float *res);
extern "C" void dslashReference(float *res, float **gaugeEven, float **gaugeOdd, float *spinorField, int oddBit);
extern "C" void printSpinorField(float *spinor);

//---------- dslash_cuda.cu

extern "C" void sendGaugeField(float **gaugeEven, float **gaugeOdd);
extern "C" void sendSpinorFieldEven(float *spinorEven);
extern "C" void sendSpinorFieldOdd(float *spinorOdd);

extern "C" void retrieveSpinorFieldEven(float *res);
extern "C" void retrieveSpinorFieldOdd(float *res);

extern "C" void initializeCuda(int argc, char** argv);
extern "C" void releaseCuda();

extern "C" void dslashCuda(int oddBit);
