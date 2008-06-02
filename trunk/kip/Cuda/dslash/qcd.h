
#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L4 64 // "time" dimension
#define L (L1*L2*L3*L4) // total number of lattice points

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (L/BLOCK_DIM) // there are L threads in total

#define SPINOR_SIZE (24) // spinors have 4*3*2 floats
#define PACKED_GAUGE_SIZE (4*20) // gauge matrices rounded up to fit float4 elements

#define SPINOR_BYTES (SPINOR_SIZE*sizeof(float))
#define PACKED_GAUGE_BYTES (PACKED_GAUGE_SIZE*sizeof(float))

//---------- qcd_gold.cpp

extern "C" void constructGaugeField(float **res);
extern "C" void constructSpinorField(float *res);
extern "C" void computeGold(float *res, float **gauge, float *spinor);

extern "C" void printSpinorField(float *spinor);
extern "C" void packGaugeField(float4 *res, float **gauge);
extern "C" void packSpinorField(float4 *res, float *spinor);
extern "C" void unpackSpinorField(float *res, float4 *spinorPacked);
