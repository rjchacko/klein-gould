#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "ising.h"


Ising1::Ising1(int len, int dim, float h, float T) : Ising(len, dim, h, T) {
    len = len;
    dim = dim;
    n = (int)powl(len, dim);
    spins = (int *)malloc(n*sizeof(int));
    
    for (int i = 0; i < n; i++) {
        spins[i] = 0;
    }
}

Ising1::~Ising1() {
    free(spins);
}

void Ising1::set(int i, int s) {
    spins[i] = s;
}

int Ising1::get(int i) {
    return spins[i];
}

int Ising1::neighborSum(int i) {
    int acc = 0;
    int len_d = 1;
    for (int d = 0; d < dim; d++) {
        int x = (i / len_d) % len;
        for (int dir = -1; dir <= 1; dir += 2) {
            int xp = (x + dir + len) % len;
            int dx = xp - x;
            acc += get(i + dx * len_d);
        }
        len_d *= len;
    }
    return acc;
}

void Ising1::completeNeighborSum(int *sum) {
    for (int i = 0; i < n; i++) {
        sum[i] = neighborSum(i);
    }
}

void Ising1::update(int parityTarget) {
    for (int i = 0; i < n; i++) {
        if (indexParity(i) == parityTarget) {
            int m = 2 * (neighborSum(i) - dim);
            int s = 2 * get(i) - 1;
            if (shouldFlipSpin(s, m)) {
                set(i, (1-s)/2);
            }
        }
    }
}

/*

void printNeighbors1(Ising1 self) {
    for (int y = 0; y < self.len; y++) {
        for (int x = 0; x < self.len; x++) {
            printf("%d ", self.sum[y*self.len+x]);
        }
        printf("\n");
    }
    printf("\n");
}


*/
