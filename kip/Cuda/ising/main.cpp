#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "ising.h"



void testSum(Ising &ising1, Ising &ising2) {
    int n = ising1.n;
    int *sum1 = (int *)malloc(n*sizeof(int));
    int *sum2 = (int *)malloc(n*sizeof(int));
    ising1.completeNeighborSum(sum1);
    ising2.completeNeighborSum(sum2);
    
    printf("Testing summation\n");
    printf("Lattices\n");
    ising1.print();
    ising2.print();
    printf("Sums\n");
    ising1.printLattice(sum1);
    ising2.printLattice(sum2);
    
    for (int i = 0; i < n; i++) {
        if (sum1[i] != sum2[i]) {
            printf("Mismatch i=%d, sum %d != %d\n\n", i, sum1[i], sum2[i]);
            return;
        }
    }
    printf("All site sums match!\n\n");
}


void testUpdate(Ising &ising1, Ising &ising2) {
    ising1.update(0);
    ising2.update(0);
    
    ising1.update(1);
    ising2.update(1);
    
    printf("Updated lattices\n");
    ising1.print();
    ising2.print();
    
    ising1.compare(&ising2);
}

int main (int argc, char * const argv[]) {
    initPentacubeParity();
    
    int len = 6;
    int dim = 7;
    float h = 0;
    float T = 0.2;
    
    Ising1 ising1 = Ising1(len, dim, h, T);
    Ising2 ising2 = Ising2(len, dim, h, T);
    
    ising1.randomizeSpins();
    ising2.randomizeSpins();
    
    testSum(ising1, ising2);
    testUpdate(ising1, ising2);
    
    return 0;
}
