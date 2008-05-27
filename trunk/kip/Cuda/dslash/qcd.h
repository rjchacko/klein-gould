
#define DIM0 8 // "time" dimension
#define DIM1 16 // "x" dimension
#define DIM2 16 // "y" dimension
#define DIM3 16 // "z" dimension
#define L (DIM0*DIM1*DIM2*DIM3) // total number of lattice points

#define BLOCK_DIM (64) // threads per block
#define GRID_DIM (L/BLOCK_DIM) // there are L threads in total

#define SPINOR_SIZE (24) // spinors have 4*3*2 floats
#define GAUGE_SIZE (4*20) // each gauge matrix has 3*3*2 floats, rounding up gives 20, and there are four x,y,z,t directions.
#define SPINOR_BYTES (SPINOR_SIZE*sizeof(float))
#define GAUGE_BYTES (GAUGE_SIZE*sizeof(float))

// get a lattice index
#define IDX(z,y,x,t) (z*(DIM0*DIM1*DIM2) + y*(DIM0*DIM1) + x*DIM0 + t)
