#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

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
    double h = 0;
    double T = 2;
    
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
    double h = 0;
    double T = 0.5;
    
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


double meanMagnetization(Ising &ising) {
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
    return (double) (m / iters);
}

void test3() {
    int len = 10;
    int dim = 2;
    double h = -0.4;
    double T = 2.0;
    
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
    double h = -0.45;
    double T = 2.269*4/9;

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

struct timeval startTime;

void stopwatchStart() {
    gettimeofday(&startTime, NULL);
}

double stopwatchStop() {
    struct timeval endTime;
    gettimeofday( &endTime, 0);

    long ds = endTime.tv_sec - startTime.tv_sec;
    long dus = endTime.tv_usec - startTime.tv_usec;
    return ds + 0.000001*dus;
}


void test5 () // Timing test
{
    int len = 8;
    int dim = 7;
    double h = -0.40;
    double T = 12.0;

    int iters = 1000;
    
    IsingCuda ic = IsingCuda(len, dim, h, T);

    stopwatchStart (); 
    for (int i=0; i<iters; ++i)
    {
        ic.update (0);
        ic.update (1);
    }
    double spinsPerSecond = (double) ic.n*iters/stopwatchStop ();
    printf ("%.2e spin evaluations per second\n", spinsPerSecond);

}

void test6 () // Find m and probe for Tc
{
    int iters = 1000;

    int len = 10;
    int dim = 7;
    double h = 0.00;
    double Tmin = 5.0;
    double Tmax = 15.0;
    double dT = 0.1;

    for (double T=Tmin; T<Tmax; T+=dT)
    {
        IsingCuda ic = IsingCuda (len, dim, h, T);
        ic.allSpinsUp ();
        for (int i=0; i<iters; ++i)
        {
            ic.update ();
        }
        printf ("T = %.2e, m after %d steps = %.3e\n", 
            ic.T, iters, ic.magnetization () / ic.n);
    }
    
}

void test7 () // testing energy calculation
{
    int len = 4;
    int dim = 7;
    double h = 0.00;
    double T = 8.50;

    int iters = 1000;

    IsingCuda ic = IsingCuda(len, dim, h, T);
    ic.allSpinsDown ();
    printf ("M=%f, E=%f\n", ic.magnetization (), ic.energy ());
    for (int i=0; i<iters; ++i)
    {
        ic.update (0);
        ic.update (1);
    }
    printf ("M=%f, E=%f\n", ic.magnetization (), ic.energy ());

}

double mean (double * arr, int length)
{
    double tot = 0;
    for (int i=0; i<length; ++i)
        tot += arr[i];
    return tot/length;
}

double abs_var (double * arr, int length)
{
    double mu = 0;
    double mu2 = 0;
    for (int i=0; i<length; ++i)
    {
        mu += fabs (arr[i]);
        mu2 += arr[i]*arr[i];
    }
    mu  /= length;
    mu2 /= length;
    return (mu2 - mu*mu);
}

void find_tc ()
{
    // Stauffer Tc ~ 12.869019123
    //int l = 6;
    //int d = 7;
    //double h = 0.00;
    //double Tmin = 11.0;
    //double Tmax = 14.0;
    //double dT = 0.01;
    int l = 256;
    int d = 2;
    double h = 0.00;
    double Tmin = 2.2;
    double Tmax = 2.4;
    double dT = 0.002;

    int relaxTime = 2000;
    int iters = 20000;
    double * M = new double [iters];
    double * E = new double [iters];

    FILE * out;
    out = fopen ("data/tc.dat", "w");
    if (out == NULL)
        printf ("Error opening file.\n");

    for (double T=Tmin; T<=Tmax; T+=dT)
    {
        IsingCuda * ic = new IsingCuda (l, d, h, T);
        //Ising1 * ic = new Ising1 (l, d, h, T);
        //Ising2 * ic = new Ising2 (l, d, h, T);
        ic->allSpinsUp ();
        for (int i=0; i<relaxTime; ++i)
            ic->update ();
        for (int i=0; i<iters; ++i)
        {
            ic->update ();
            M[i] = ic->magnetization () / ic->n;
            E[i] = ic->energy () / ic->n;
        }
        if (!(((int) ((T-Tmin)/dT))%10))
            printf ("finished trial %d/%d\n", (int) ((T-Tmin)/dT), (int)((Tmax-Tmin)/dT));
        double chi_n = abs_var (M, iters)/T;
        double c_n = abs_var (E, iters)/T/T;
        fprintf (out, "%f\t%f\t%f\t%f\t%f\n", T, chi_n, fabs(mean (M, iters)), c_n, mean (E, iters));
    }
    printf ("\n");

    delete M; delete E; fclose (out);
}

int main (int argc, char *argv[]) {
    initCuda(argc, argv);
   
    //test1();
    //test2();
    //test3();
    //test4 ();
    //test5 ();
    //test6 ();
    //test7 ();
    find_tc ();
    
    return 0;
}
