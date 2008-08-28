

// ----------------------------------------------------------------------------
// NVIsing
//

typedef unsigned int bits32;

typedef struct {
    // a lattice containing (n = len^dim) spins
    int len, dim, n;
    unsigned int *spins; // n spins stored using n/32 unsigned ints
    float h, T; // external field and temperature
} NVIsing;


NVIsing nvAllocate(int len, int dim, int n);
void    nvFree(NVIsing ising);
void    nvCopySpins(NVIsing ising, unsigned int *spins);
void    nvUpdate(NVIsing ising, int parityTarget);
