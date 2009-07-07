#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "generate.h"

generate::generate (IsingCuda * ic, int theSeed)
{
    this->ic = ic;
    this->theSeed = theSeed;
}

int generate::initialGuess ()
{
    int t = 0;
    ic->rngSeed (theSeed);
    ic->allSpinsUp (); ic->upH ();
    for (int i=0; i<upRelax; ++i) ic->update ();
    ic->flipH ();
    while (ic->magnetization ()/ic->n > thresh) { ++t; ic->update (); }
    return t;
}

void generate::genDroplet ()
{
    genDropletLog (initialGuess (), dropTimeSearch, 2);
}

void generate::genDropletLog (int t, int search, int twoPow)
{
    int numNucleate = 0;
    const double limit = 0.1;

    // for hypothesis test
    double pbar;
    double sigma = 0.5;
    double z_alpha = 1.96;

    // Prepare ising class
    ic->rngSeed (theSeed);
    ic->allSpinsUp (); ic->upH ();
    for (int i=0; i<upRelax; ++i) ic->update ();
    ic->flipH ();

    for (int i=0; i<t; ++i) ic->update (); // regenerate the state
    int * s = new int [ic->n]; // yeah... not so efficient, but whatever
    ic->getSpins (s); // save the intermediate state

    ic->rngSeed (time (NULL)); // randomize the seed
    for (int i=0; i<numTests; ++i)
    {
        ic->setSpins (s); // reset the intermediate state
        for (int j=0; j<nucCheck; ++j)
        {
            ic->update ();
            if (ic->magnetization () / ic->n < thresh) break;
        }
        if (ic->magnetization () / ic->n < thresh) ++numNucleate;
    }

    pbar = (double) numNucleate/numTests;
    if ((double) search/twoPow < limit && pbar > 0.5)
    {
        ic->setSpins (s);
        return;
    }

    if (pbar > 0.5 + z_alpha*sigma/sqrt (numTests))
        genDropletLog (t-search/twoPow-1, search, twoPow*2);
    else if (pbar < 0.5 - z_alpha*sigma/sqrt (numTests))
        genDropletLog (t+search/twoPow+1, search, twoPow*2);
    else // restore droplet
        ic->setSpins (s);
    delete s;
}
