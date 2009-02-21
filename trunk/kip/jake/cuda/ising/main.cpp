#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h>

#include "ising.h"
#include "nt.h"
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
    //int l = 8;
    //int d = 7;
    //double h = 0.00;
    //double Tmin = 11.0;
    //double Tmax = 14.0;
    //double dT = 0.005;
    //int l = 512;
    //int d = 2;
    //double h = 0.001;
    //double Tmin = 2.2;
    //double Tmax = 2.4;
    //double dT = 0.001;
    int l = 10;
    int d = 5;
    double h = 0.00;
    double Tmin = 5.0;
    double Tmax = 10.0;
    double dT = 0.1;

    int relaxTime = 1000;
    int iters = 1000;
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

void find_tc2 ()
{
    // // Stauffer Tc ~ 12.869019123
    //int l = 4;
    //int d = 7;
    //double h = 0.00;
    //double Tmin = 11.0;
    //double Tmax = 14.0;
    //double dT = 0.01;
    //int l = 64;
    //int d = 2;
    //double h = 0.00;
    //double Tmin = 2.2;
    //double Tmax = 2.4;
    //double dT = 0.005;
    //int l = 10;
    //int d = 5;
    //double h = 0.00;
    //double Tmin = 8.7;
    //double Tmax = 8.85;
    //double dT = 0.001;
    // // Stauffer Tc ~ 10.8348231
    int l = 12;
    int d = 6;
    double h = 0.00;
    double Tmin = 10.81;
    double Tmax = 10.84;
    double dT = 0.002;

    int relaxTime = 750;
    int iters = 400;

    double t_m = 0.0;
    double t_e = 0.0;
    double sum_m = 0.0;
    double sum_e = 0.0;
    double sum_m2 = 0.0;
    double sum_e2 = 0.0;
    double t_chi = 0.0;
    double t_cv = 0.0;
    double chi = 0.0;
    double cv = 0.0;
    double chi2 = 0.0;
    double cv2 = 0.0;
    double sum_chi = 0.0;
    double sum_cv = 0.0;
    double sum_chi2 = 0.0;
    double sum_cv2 = 0.0;
    double sig_cv_err;
    double sig_chi_err;

    double err_lim = 0.1;
    int min_samples = 15;
    int max_samples = 40;
    int nsamples = 0;

    FILE * out = NULL; 
    out = fopen ("data/tc.dat", "w");
    if (out == NULL)
        printf ("Error opening file.\n");

    for (double T=Tmin; T<=Tmax; T+=dT)
    {
        IsingCuda * ic = new IsingCuda (l, d, h, T);
        //Ising1 * ic = new Ising1 (l, d, h, T);
        //Ising2 * ic = new Ising2 (l, d, h, T);
        ic->allSpinsUp ();

        sum_chi = sum_cv = sum_chi2 = sum_cv2 = 0.0;
        nsamples = 0;
        for (int i=0; i<relaxTime; ++i)
            ic->update ();
        do
        {
            sum_m = sum_e = sum_m2 = sum_e2 = 0.0;
            for (int i=0; i<iters; ++i)
            {
                ic->update ();
                t_m = fabs (ic->magnetization ()) / ic->n;
                t_e = ic->energy () / ic->n;
                sum_m += t_m; sum_e += t_e;
                t_m *= t_m; t_e *= t_e;
                sum_m2 += t_m; sum_e2 += t_e;
            }
            sum_m /= iters; sum_e /= iters;
            sum_m2 /= iters; sum_e2 /= iters;

            ++nsamples;
            
            t_chi = (sum_m2 - sum_m*sum_m)/ic->T;
            t_cv = (sum_e2 - sum_e*sum_e)/ic->T/ic->T;
            sum_chi += t_chi; sum_cv += t_cv;
            t_chi *= t_chi; t_cv *= t_cv;
            sum_chi2 += t_chi; sum_cv2 += t_cv;
            
            chi = sum_chi/nsamples; cv = sum_cv/nsamples;
            chi2 = sum_chi2/nsamples; cv2 = sum_cv2/nsamples;

            sig_chi_err = sqrt ((chi2 - chi*chi)/nsamples);
            sig_cv_err = sqrt ((cv2 - cv*cv)/nsamples);
            //printf ("%.4e %.4e\n", chi, sig_chi_err);
            //printf ("%.4e %d\n", sig_chi/chi, nsamples);
            if (nsamples > max_samples) break;
        } while (sig_chi_err/chi > err_lim || nsamples < min_samples);

        if (!(((int) ((T-Tmin)/dT))%10))
            printf ("finished trial %d/%d\n", (int) ((T-Tmin)/dT), (int)((Tmax-Tmin)/dT));
        //double chi_n = abs_var (M, dataLen)/T;
        //double c_n = abs_var (E, dataLen)/T/T;
        fprintf (out, "%e\t%e\t%e\t%e\t%e\n", T, chi, sig_chi_err, cv, sig_cv_err);
    }
    printf ("\n");

    fclose (out);
}

void mgraph ()
{
    const int l = 8;
    const int d = 7;
    const double T = 4*12.8690191/9;
    const double h = 4.0215;

//    double mu_m, std_m;
    
    IsingCuda ic = IsingCuda (l, d, h, T);

//    est_mean_m (ic, mu_m, std_m);

//    int relaxTime = 50;


    char name [100];

    int ntrials = 4;
    int stepMax = 1000;

    for (int i=1; i<ntrials+1; ++i)
    {
        sprintf (name, "data/m%d.dat", i);
        
        FILE * out = fopen (name, "w");
        assert (out);

        ic.allSpinsUp (); // Prepare the lattice
        ic.downH (); // now, m*h < 0
        
        int step = 0;
        double m;
        do
        {   
            ic.update ();
            ++step;
            m = ic.magnetization ()/ic.n;
            fprintf (out, "%d\t%.5e\n", step, m);
        //} while (fabs(m-mu_m) < 5*std_m || step < 15);
        } while (m*ic.h < 0 && step < stepMax+1);
        
        fclose (out);
    }
}

void find_hs ()
{
    // // Stauffer 7d Tc ~ 12.8690191
    // // Stauffer 6d Tc ~ 10.8348231
    // //          5d Tc ~ 8.77847518
    int d = 5;
    int l = 10;
    const double T = 8.77847518*4/9;
    double h0 = 2.4;
    double hf = 3.0;
    double dh = 0.01;

    int relaxTime = 50;
    int iters = 500;

    double t_m, t_e, sum_m, sum_e, sum_m2, sum_e2, t_chi, t_cv, chi, cv,
    chi2, cv2, sum_chi, sum_cv, sum_chi2, sum_cv2, sig_cv_err, sig_chi_err;

    double err_lim = 0.05;
    int min_samples = 15;
    int max_samples = 50;
    int nsamples = 0;
    int nnucleate = 0;

    FILE * out = NULL; 
    out = fopen ("data/hs.dat", "w");
    if (out == NULL)
        printf ("Error opening file.\n");

    for (double h=h0; h<=hf; h+=dh)
    {
        IsingCuda * ic = new IsingCuda (l, d, h, T);
        //Ising1 * ic = new Ising1 (l, d, h, T);
        //Ising2 * ic = new Ising2 (l, d, h, T);

        sum_chi = sum_cv = sum_chi2 = sum_cv2 = 0.0;
        nsamples = 0;
        nnucleate = 0;

        do
        {
            ic->allSpinsUp ();
            ic->upH (); // Put field up
            for (int i=0; i<relaxTime; ++i)
                ic->update ();

            ic->downH (); // Put field down

            sum_m = sum_e = sum_m2 = sum_e2 = 0.0;
            for (int i=0; i<iters; ++i)
            {
                ic->update ();
                t_m = fabs (ic->magnetization ()) / ic->n;
                t_e = ic->energy () / ic->n;
                sum_m += t_m; sum_e += t_e;
                t_m *= t_m; t_e *= t_e;
                sum_m2 += t_m; sum_e2 += t_e;
            }

            if ( ic->h * sum_m > 0 )
            {
                ++nnucleate;
                //continue; // Nucleated -- ignore trial
            }

            sum_m /= iters; sum_e /= iters;
            sum_m2 /= iters; sum_e2 /= iters;

            ++nsamples;
            
            t_chi = (sum_m2 - sum_m*sum_m)/ic->T;
            t_cv = (sum_e2 - sum_e*sum_e)/ic->T/ic->T;
            sum_chi += t_chi; sum_cv += t_cv;
            t_chi *= t_chi; t_cv *= t_cv;
            sum_chi2 += t_chi; sum_cv2 += t_cv;
            
            chi = sum_chi/nsamples; cv = sum_cv/nsamples;
            chi2 = sum_chi2/nsamples; cv2 = sum_cv2/nsamples;

            sig_chi_err = sqrt ((chi2 - chi*chi)/nsamples);
            sig_cv_err = sqrt ((cv2 - cv*cv)/nsamples);
            //printf ("%.4e %.4e\n", chi, sig_chi_err);
            //printf ("%.4e %d\n", sig_chi/chi, nsamples);
            if (nsamples > max_samples) break;
            if (nnucleate > max_samples/2) break;
        } while (sig_chi_err/chi > err_lim || nsamples < min_samples);

        if (!(((int) ((h-h0)/dh))%10))
            printf ("finished trial %d/%d\n", (int) ((h-h0)/dh), (int)((hf-h0)/dh));
        //double chi_n = abs_var (M, dataLen)/T;
        //double c_n = abs_var (E, dataLen)/T/T;
        fprintf (out, "%e\t%e\t%e\t%e\t%e\n", h, chi, sig_chi_err, cv, sig_cv_err);
    }
    printf ("\n");

    fclose (out);
}

void ntControl ()
{    
    // // 7d Tc ~ 12.8690191
    // // 6d Tc ~ 10.8348231
    // // 5d Tc ~ 8.77847518
    
    //const int l = 10;
    //const int d = 5;
    //const double h = 2.558;
    //const double T = 8.77847518*4/9;
    //nt n = nt (6, d, h, T);
    //for (int l=8; l<=24; l+=2)
    //{
    //    n.lChange (l);
    //    n.nucleationTimes ();
    //}

    //const int l = 12;
    const int d = 7;
    const double T = 4*12.8690191/9;
    const double h = 4.015;
    //const double h = 4.023;
    nt n = nt (6, d, h, T);
    for (int l=6; l<=14; l+=2)
    {
        n.lChange (l);
        n.sim ();
    }
}

void chi ()
{
    //int l = 10; int d = 7; double T = 4*12.869/9; 
    //double h = 4.019; double delta_h = 0.0003;
    
    int l=16; int d=5; double T=8.778*4/9; 
    double h=2.53; double delta_h=0.002;

    nt n = nt (l, d, h, T);
    for (int i=0; i<30; ++i)
    {
        n.hChange (h+i*delta_h);
        n.sim (3000);
    }
}

int main (int argc, char *argv[]) {
    initCuda(argc, argv);
   
    //test1 ();
    //test2 ();
    //test3 ();
    //test4 ();
    //test5 ();
    //test6 ();
    //test7 ();
    //find_tc ();
    //find_tc2 ();
    //find_hs ();
    //mgraph ();
    //ntControl ();
    chi ();
    
    return 0;
}
