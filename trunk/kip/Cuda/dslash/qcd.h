
#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L4 64 // "time" dimension

#define L1h (L1/2) // half of the full "x" dimension, useful for even/odd lattice indexing
#define L (L1h*L2*L3*L4) // poorly named; total number of even/odd lattice points


// ---------- dslash_reference.cpp ----------

extern "C" void constructUnitGaugeField(float **resEven, float **resOdd);
extern "C" void constructGaugeField(float **resEven, float **resOdd);
extern "C" void constructPointSpinorField(float *resEven, float *resOdd, int i0, int s0, int c0);
extern "C" void constructSpinorField(float *res);

extern "C" void dslashReference(float *res, float **gaugeEven, float **gaugeOdd, float *spinorField, int oddBit);

extern "C" void printSpinorElement(float *spinorEven, float *spinorOdd, int X);
extern "C" void printSpinorHalfField(float *spinor);


// ---------- dslash_cuda.cu ----------

extern "C" void sendGaugeField(float **gaugeEven, float **gaugeOdd);
extern "C" void sendSpinorFieldEven(float *spinorEven);
extern "C" void sendSpinorFieldOdd(float *spinorOdd);

extern "C" void retrieveSpinorFieldEven(float *res);
extern "C" void retrieveSpinorFieldOdd(float *res);

extern "C" void initializeCuda(int argc, char** argv);
extern "C" void releaseCuda();

extern "C" void dslashCuda(int oddBit);
