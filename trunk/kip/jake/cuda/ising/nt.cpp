#include <cmath>
#include <cstdio>
#include <assert.h>
#include "ising.h"
#include "nt.h"

nt::nt (int l, int d, double h, double T)
{
    ic = new IsingCuda (l, d, h, T);

    char mName [100];
    char ntName [100];
    sprintf (mName, "./data/res/%dd_m.dat", d);
    sprintf (ntName, "./data/res/%dd_nt.dat", d);

    mOut = fopen (mName, "wb"); assert (mOut);
    ntOut = fopen (ntName, "w"); assert (ntOut);
    headFlag =  10000000.0f;

    // Initialize m header
    head[0]      = headFlag;
    head[HLEN-1] = headFlag;
    head[1]      = (float) l;
    head[2]      = (float) d;
    head[3]      = (float) h;
    head[4]      = (float) T;
}

nt::~nt ()
{
    delete ic;
    fclose (mOut);
    fclose (ntOut);
}

void nt::lChange (int l)
{
    int d = ic->dim; double h = ic->h; double T = ic->T;
    delete ic;
    ic = new IsingCuda (l, d, h, T);
    head[1] = (float) l;
}

void nt::estMean_m (double & mu_m, double & std_m)
{
    /**
     * Estimate mean magentization. Used as a tool to help find nucleation
     * time.
     */

    int relaxTime = 25;
    int iters = 200; //int iters = 50;
    int nucTestTime = 30; // int nucTestTime = 15;
    double m, m2;

    while (1)
    {
        ic->allSpinsUp (); ic->downH ();
        m = m2 = 0.0;
        for (int i=0; i<relaxTime; ++i) ic->update (); // Settle into metast.
        for (int i=0; i<iters; ++i)
        {
            ic->update ();
            double t = ic->magnetization ()/ic->n;
            m += t; t *= t; m2 += t;
        }
        for (int i=0; i<nucTestTime; ++i) ic->update (); // Make sure no nuc.
        if (ic->magnetization () * ic->h < 0) // Not nucleated
        {
            mu_m = m/iters;
            std_m = sqrt (m2/iters - m*m/iters/iters );
            return;
        }
        else if (iters > 50) --iters; // Try shorter this time
        else
        {
            // handle case where we had nucleation too early
        }
    }
}

void nt::nucleationTimes ()
{
    double mu_m, std_m;

    int relaxTime = 25;
    int afterTime = 25;
    int ntrials = 1600;
    int histSize = 1000;

    const int nbins = 25; //histSize/25;
    int * bins = new int [nbins];

    double m;

    double meanNucleation = 0;
    double meanNucleationError = 0;
    int step;

    int decays = 0;

    estMean_m (mu_m, std_m); // Get estimates of mean and std

    #ifdef MOUT
    fwrite (head, sizeof(float), HLEN, mOut);
    #endif

    for (int i=0; i<ntrials; ++i)
    {
        ic->allSpinsUp (); ic->downH (); // Prepare the lattice
        for (int j=0; j<relaxTime; ++j) // Enter metastable state
            ic->update ();
        
        step = 0;
        m = ic->magnetization () / ic->n;
        while (fabs(m - mu_m) < 5*std_m)
        {   
            ic->update ();
            m = ic->magnetization () / ic->n;
            ++step;
            #ifdef MOUT
            fwrite (&m, sizeof(float), 1, mOut);
            #endif
        }
        
        for (int j=0; j<afterTime; ++j) // Make sure metastable decayed
            ic->update ();
        
        // ignore premature decays and ones which were false positives
        if (0 < step && ic->magnetization ()*ic->h > 0) 
        {
            ++decays;
            meanNucleation += (double) step;
            if (step < histSize+1) // Does not exceed histogram
                ++ bins [(step-1)*nbins/histSize]; 
        }
        
        printf ("\r                                                \r");
        printf ("(L=%d, d=%d) Finished %d/%d", ic->len, ic->dim, i, ntrials);
        fflush (stdout);
    }
    printf ("\n");

    meanNucleation /= decays;
    meanNucleationError = meanNucleation / sqrt(ntrials);
    
    fprintf (ntOut, "%d\t%.5e\t%.5e\n", 
        ic->len, meanNucleation, meanNucleationError);
    fflush (ntOut);

    delete [] bins;
}
