
#define L0 16 // "time" dimension
#define L1 16 // "x" dimension
#define L2 16 // "y" dimension
#define L3 16 // "z" dimension
#define L (L0*L1*L2*L3) // total number of lattice points

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (L/BLOCK_DIM) // there are L threads in total

#define SPINOR_SIZE (24) // spinors have 4*3*2 floats
#define GAUGE_SIZE (4*18) // each gauge matrix has 3*3*2 floats and there are four x,y,z,t directions.
#define PACKED_GAUGE_SIZE (4*20) // gauge matrices rounded up to fit float4 elements

#define SPINOR_BYTES (SPINOR_SIZE*sizeof(float))
#define GAUGE_BYTES (GAUGE_SIZE*sizeof(float))
#define PACKED_GAUGE_BYTES (PACKED_GAUGE_SIZE*sizeof(float))

//---------- qcd_gold.cpp

extern "C" void constructGaugeField(float *res);
extern "C" void constructSpinorField(float *res);
extern "C" void packGaugeField(float4 *res, float *gauge);
extern "C" void packSpinorField(float4 *res, float *spinor);
extern "C" void unpackSpinorField(float *res, float4 *spinorPacked);
extern "C" void printSpinorField(float *spinor);
extern "C" void testSpinorField(float *spinor);

extern "C" void computeGold(float* res, float* gauge, float *spinor);

