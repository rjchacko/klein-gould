#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

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
        //set(i, (i%5+i%4)%2);
        set(i, rand()%2);
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

void Ising::flipH ()
{
    h = -h;
}

void Ising::downH ()
{
    h = (h < 0) ? h : -h;
}

void Ising::upH ()
{
    h = (h > 0) ? h : -h;
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

// JEE

void Ising::saveIsing (char * filename)
{
    FILE * fp = fopen (filename, "wb");
    assert (fp);
    int * s = new int [n];
    char * sc = new char [n];
    getSpins (s);
    for (int i=0; i<n; ++i) sc[i] = s[i];
    double params [] = {(double) len, (double) dim, h, T};
    
    // Reverse the parameter bytes so java program can read them
    char * f, * ftemp;
    double * fpoint;
    ftemp = new char [8];
    for (int i=0; i<4; ++i)
    {
        f = (char *) &params[i];
        for (int j=0; j<8; ++j)
            ftemp[j] = f[8-j-1];
        fpoint = (double *) ftemp;
        params[i] = *fpoint;
    }

    fwrite (params, sizeof(double), 4, fp);
    fwrite (sc, sizeof(char), n, fp);
    fclose (fp);
    delete ftemp; 
    delete s; delete sc;
}
