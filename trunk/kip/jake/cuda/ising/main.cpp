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
    double m = 0;
    double iters = 10000;
    double mt;
    for (int i = 0; i < iters; i++) {
        ising.update(0);
        ising.update(1);
        mt = ising.magnetization ();
        if (mt > ising.n || mt < -ising.n)
            printf ("Error: magnetization measurement out of range.");
        m += mt;
    }
    return (float) (m / iters);
}

void test3() {
    int len = 10;
    int dim = 2;
    float h = -0.4;
    float T = 2.0;
    
    IsingCuda ic = IsingCuda(len, dim, h, T);
    //Ising2 ih = Ising2(len, dim, h, T);

    printf("Cuda: <m> = %.3f\n", meanMagnetization(ic));
    printf("%d\n", ic.n);
    //printf("Host: <m> = %.3f\n", meanMagnetization(ih));
}

void test4 ()
{
    const char basename[] = "data/trial";
    const char ext[] = "cud";
    char out_name[100];
    int ntrials = 50;

    int len = 128;
    int dim = 2;
    float h = -0.51;
    float T = 2.269*4/9;

    FILE * out;
    double m;

    IsingCuda ic = IsingCuda (len, dim, h, T);

    for (int i=0; i<ntrials; ++i)
    {
        // Initialize data output
        sprintf (out_name, "%s_%.2f_%d.%s", basename, -h, i, ext);
        out = fopen (out_name, "w");
       
        // Initialize the run
        int step = 0;
        ic.allSpinsUp ();
        do
        {
            m = ic.magnetization () / ic.n;
            // Data output to be compatible with my other script
            fprintf (out, "%f\t%f\t%f\t%f\t%d\n", h, T, m, m*m, step++);
            ic.update (0);
            ic.update (1);
        }
        while (m > 0);
        fclose (out);
    }
}

void test5 ()
{
    int len = 10;
    int dim = 7;
    float h = -0.40;
    float T = 2.0;

    int iters = 10000;
    
    IsingCuda ic = IsingCuda(len, dim, h, T);

    for (int i=0; i<iters; ++i)
    {
        ic.update (0);
        ic.update (1);
    }

}

int main (int argc, char *argv[]) {
    initCuda(argc, argv);
    
    //test1();
    //test2();
    //test3();
    //test4 ();
    test5 ();
    
    return 0;
}
