#ifndef GENERATE_H
#define GENERATE_H

#include "ising.h"

class generate
{
public:
    static const int numTests = 100;
    static const int upRelax = 15;

    // BELOW IS SUBJECT TO CHANGE DEPENDING ON r OF SYSTEM
    static const int nucCheck = 25; // number of steps to allow m<thresh
    static const double thresh = 0.75;

    // The amound of time possibilities for the log-generate function
    static const int dropTimeSearch = 64;

    IsingCuda * ic;
    int theSeed;

    generate (IsingCuda *, int);
    void genDroplet ();

private:
    int initialGuess ();
    void genDropletLog (int, int, int);
};

#endif // GENERATE_H
