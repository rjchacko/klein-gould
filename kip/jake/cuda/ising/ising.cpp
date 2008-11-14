#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "ising.h"


Ising::Ising(int len, int dim, float h, float T)
    : len(len), dim(dim), n((int)powl(len, dim)), h(h), T(T) {}


double Ising::magnetization() {
    int ret = 0;
    for (int i = 0; i < n; i++)
        ret += get(i);
    return 2.0*ret-n;
}

void Ising::randomizeSpins() {
    // srand(0);
    for (int i = 0; i < n; i++) {
        set(i, (i%5+i%4)%2);
        //set(i, rand()%2);
        //set(i, 0);
    }
    transferHostToDevice();
}

void Ising::allSpinsUp ()
{
    for (int i=0; i<n; ++i)
        set (i, 1);
    transferHostToDevice ();
}

void Ising::allSpinsDown ()
{
    for (int i=0; i<n; ++i)
        set (i, 0);
    transferHostToDevice ();
}

void Ising::setSpins(int *spins) {
    for (int i = 0; i < n; i++) {
        set(i, spins[i]);
    }
    transferHostToDevice();
}

void Ising::getSpins(int *spins) {
    transferDeviceToHost();
    for (int i = 0; i < n; i++) {
        spins[i] = get(i);
    }
}

void Ising::printLattice(int *a) {
    for (int y = 0; y < ((dim==1)?1:len); y++) {
        for (int x = 0; x < len; x++) {
            printf("%d ", a[y*len+x]);
        }
        printf("\n");
    }
    printf("\n");
}

void Ising::print() {
    transferDeviceToHost();
    for (int y = 0; y < ((dim==1)?1:len); y++) {
        for (int x = 0; x < len; x++) {
            printf("%d ", get(y*len+x));
        }
        printf("\n");
    }
    printf("\n");
}

void Ising::compare(Ising *that) {
    transferDeviceToHost();
    that->transferDeviceToHost();
    for (int i = 0; i < n; i++) {
        if (get(i) != that->get(i)) {
            printf("Lattices are not equal at %d\n", i);
            return;
        }
    }
    printf("Lattices are equal\n");
}

int Ising::shouldFlipSpin(int s, int m) {
    float dE = 2*s*(m + h);
    if (dE < 0)
        return 1;
    else {
#ifdef DETERMINISTIC
        float r = 0.1;
#else
        float r = (float)lrand48() / (unsigned int)(1<<31);
#endif
        return exp(- dE / T) > r;
    }
}

int Ising::indexParity(int i) {
    int acc = 0;
    int len_d = 1;
    for (int d = 0; d < dim; d++) {
        int x = (i / len_d) % len;
        acc += x;
        len_d *= len;
    }
    return acc % 2;
}
