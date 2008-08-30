#include <math.h>
#include "ising.h"

NVIsing nvAllocate(int len, int dim) {
    NVIsing ret;
    ret.len = len;
    ret.dim = dim;
    ret.n = (int)powl(len, dim);
    cudaMalloc((void**)&(ret.spins), (ret.n/32)*sizeof(unsigned int));
    return ret;
}

void nvFree(NVIsing ising) {
    cudaFree(ising.spins);
}

void nvLoadSpins(NVIsing ising, unsigned int *spins) {
    cudaMemcpy(ising.spins, spins, (ising.n/32)*sizeof(unsigned int), cudaMemcpyHostToDevice);
}

void nvRetrieveSpins(NVIsing ising, unsigned int *spins) {
    cudaMemcpy(spins, ising.spins, (ising.n/32)*sizeof(unsigned int), cudaMemcpyDeviceToHost);
}





void nvUpdate(NVIsing ising, int parityTarget) {
    
}
