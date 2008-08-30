#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "ising.h"



// ------------------------------------------------------------------------------------------
// Main
//

int main (int argc, char * const argv[]) {
    initPentacubeParity();
    
    int len = 6;
    int dim = 6;
    float h = 0;
    float T = 0.2;
    
    Ising1 ising1 = Ising1(len, dim, h, T);
    ising1.randomizeSpins();
    ising1.update(0);
    ising1.update(1);
    ising1.print();

    Ising2 ising2 = Ising2(len, dim, h, T);
    ising2.randomizeSpins();
    ising2.update(0);
    ising2.update(1);
    ising2.print();
    
    ising1.compare(&ising2);
    return 0;
}
