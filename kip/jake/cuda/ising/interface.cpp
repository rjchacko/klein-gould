#include "interface.h"
#include "ising.h"

#include <stdio.h>

void ising_cuda_init_device() {
    char *argv[] = {""};
    initCuda(1, argv);
}

Ising *ising_cuda_create(int len, int dim, float h, float T) {
    return new IsingCuda(len, dim, h, T);
}

void ising_cuda_destroy(Ising *ising) {
    delete ising;
}

unsigned int ising_cuda_magnetization(Ising *ising) {
    // to do
    return 0;
}

void ising_cuda_get_spins(Ising *ising, int *spins) {
    ising->getSpins(spins);
}

void ising_cuda_set_spins(Ising *ising, int *spins) {
    ising->setSpins(spins);
}

void ising_cuda_step(Ising *ising, int parity) {
    ising->update(parity);
}

void ising_cuda_quench(Ising *ising, float h, float T) {
    ising->h = h;
    ising->T = T;
}
