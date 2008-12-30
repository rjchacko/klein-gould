#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "ising.h"
#include "tests.h"

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
    int l = 8;
    int d = 7;
    double h = 0.00;
    double Tmin = 11.0;
    double Tmax = 14.0;
    double dT = 0.005;
    //int l = 256;
    //int d = 2;
    //double h = 0.00;
    //double Tmin = 2.2;
    //double Tmax = 2.4;
    //double dT = 0.002;

    int relaxTime = 1000;
    int iters = 4000;
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
        fprintf (out, "%e\t%e\t%e\t%e\t%e\n", T, chi_n, fabs(mean (M, iters)), c_n, mean (E, iters));
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
