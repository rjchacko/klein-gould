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
    printf("Neighbor sums are equal\n\n");
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


void test1() {
    int len = 6;
    int dim = 7;
    float h = 0;
    float T = 2;
    
    Ising1 ising1 = Ising1(len, dim, h, T);
    Ising2 ising2 = Ising2(len, dim, h, T);
    
    ising1.randomizeSpins();
    ising2.randomizeSpins();
    
    testSum(ising1, ising2);
    testUpdate(ising1, ising2);
}


void test2() {
    int len = 16;
    int dim = 6;
    float h = 0;
    float T = 0.5;
    
    Ising2 ising1 = Ising2(len, dim, h, T);
    IsingCuda ising2 = IsingCuda(len, dim, h, T);

    printf("Randomizing %d spins\n", ising1.n);
    ising1.randomizeSpins();
    ising2.randomizeSpins();

    int iters = 2;
    
    printf("Updating 1\n");
    srand48(0);
    for (int i = 0; i < iters; i++) {
        ising1.update(0);
        ising1.update(1);
    }
    
    printf("Updating 2\n");
    for (int i = 0; i < iters; i++) {
        ising2.update(0);
        ising2.update(1);
    }

    printf("\nLattices\n");
    ising1.print();
    ising2.print();
    
    ising1.compare(&ising2);
    
    printf("Magnetizations: %g %g\n", ising1.magnetization(), ising2.magnetization());
}


float meanMagnetization(Ising &ising) {
    float m = 0;
    float iters = 100000;
    for (int i = 0; i < iters; i++) {
        ising.update(0);
        ising.update(1);
        m += ising.magnetization();
    }
    return m / iters;
}

void test3() {
    int len = 10;
    int dim = 2;
    float h = -0.5;
    float T = 2.0;
    
    IsingCuda ic = IsingCuda(len, dim, h, T);
    Ising2 ih = Ising2(len, dim, h, T);

    printf("Cuda: <m> = %.3f\n", meanMagnetization(ic));
    printf("Host: <m> = %.3f\n", meanMagnetization(ih));
}


int main (int argc, char *argv[]) {
    initCuda(argc, argv);
    
    //test1();
    //test2();
    test3();
    
    return 0;
}
