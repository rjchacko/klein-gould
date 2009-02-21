#include <cmath>
#include <cstdio>
#include <assert.h>
#include "ising.h"
#include "nt.h"

nt::nt (int l, int d, double h, double T, char * base, int ovr)
{
    relaxTime = 40;

    ic = new IsingCuda (l, d, h, T);

    char mName [100];
    char ntName [100];
    sprintf (mName, "./data/res/%s%dd_m.dat", base, d);
    sprintf (ntName, "./data/res/%s%dd_nt.dat", base, d);

    if (ovr)
    {
        mOut = fopen (mName, "wb"); assert (mOut);
        ntOut = fopen (ntName, "w"); assert (ntOut);
    }
    else
    {
        ntOut = fopen (ntName, "a"); assert (ntOut);
        mOut = fopen (mName, "ab"); assert (mOut);
    }

    // Initialize m header
    head[0]      = HEADFLAG;
    head[1]      = (float) l;
    head[2]      = (float) d;
    head[3]      = (float) h;
    head[4]      = (float) T;
    // Initialize m footer
    foot[0]      = FOOTFLAG;
}

nt::nt (int l, int d, double h, double T)
{
    nt (l, d, h, T, "", 1);
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

void nt::hChange (double h)
{
    int d = ic->dim; int l = ic->len; double T = ic->T;
    delete ic;
    ic = new IsingCuda (l, d, h, T);
    head[3] = (float) h;
}

void nt::estMean_m (double & mu_m, double & std_m)
{
    /**
     * Estimate mean magentization. Used as a tool to help find nucleation
     * time.
     */

    int iters = 250; //int iters = 50;
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
        else if (iters > 100) --iters; // Try shorter this time
        else
        {
            // handle case where we had nucleation too early
        }
    }
}

void nt::sim ()
{
    int ntrials = 1600;
    sim (ntrials);
}

void nt::sim (int ntrials)
{
    double mu_m, std_m;

    int afterTime = 15;
    //int ntrials = 1600;
    int histSize = 1000;

    const int nbins = 25; //histSize/25;
    int * bins = new int [nbins];

    float m;

    double meanNucleation = 0;
    double meanNucleationError = 0;
    double r, rError;
    int step;

    int decays = 0;

    estMean_m (mu_m, std_m); // Get estimates of mean and std

    for (int i=0; i<ntrials; ++i)
    {
        #ifdef MOUT
        fwrite (head, sizeof(float), HLEN, mOut);
        #endif

        ic->allSpinsUp (); ic->downH (); // Prepare the lattice
        for (int j=0; j<relaxTime; ++j) // Enter metastable state
            ic->update ();
        
        step = 0;
        m = ic->magnetization () / ic->n;
        while (fabs(m - mu_m) < 5*std_m)
        {   
            ic->update ();
            m = (float) ic->magnetization () / ic->n;
            ++step;
            // All metastable data is good, so we do not filter any as is
            // done with nucleation times.
            #ifdef MOUT 
            fwrite (&m, sizeof(float), 1, mOut); 
            #endif
        }
        
        // ignore premature decays and ones which were false positives
        if (0 < step) 
        {
            for (int j=0; j<afterTime; ++j) // Make sure metastable decayed
                ic->update ();

            if (ic->magnetization ()*ic->h > 0)
            {
                ++decays;
                meanNucleation += (double) step;
                if (step < histSize+1) // Does not exceed histogram
                    ++ bins [(step-1)*nbins/histSize]; 
            }
        }
        
        printf ("\r                                                     \r");
        printf ("(L=%d, d=%d) Trial %d/%d :: %d decays", 
            ic->len, ic->dim, i+1, ntrials, decays);
        fflush (stdout);
        
        #ifdef MOUT
        fwrite (foot, sizeof(float), FLEN, mOut);
        #endif
    }
    printf ("\n");

    meanNucleation /= decays;
    meanNucleationError = meanNucleation / sqrt(decays);
    r = 1/meanNucleation;
    rError = 1/(meanNucleation-meanNucleationError) - 1/meanNucleation;
    
    fprintf (ntOut, "%d\t%.5e\t%.5e\t%.5e\t%.5e\t%.5e\n", 
        ic->len, fabs (ic->h), meanNucleation, meanNucleationError, r, rError);
    fflush (ntOut);

    delete [] bins;
}
