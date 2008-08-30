#include <math.h>

// ----------------------------------------------------------------------------
// Ising (Abstract)
//

class Ising {
public:
    int len, dim, n;
    float h, T;
    
    Ising(int len, int dim, float h, float T)
            : len(len), dim(dim), n((int)powl(len, dim)), h(h), T(T) {
    }
    
    virtual ~Ising() {
    }
    
    void randomizeSpins() {
        for (int i = 0; i < n; i++) {
            set(i, (i%5+i%4)%2);
        }
    }
    
    virtual void set(int i, int s) = 0;
    virtual int get(int i) = 0;
    virtual void completeNeighborSum(int *sum) = 0;
    virtual void update(int parityTarget) = 0;
    
    void setSpins(int *spins) {
        for (int i = 0; i < n; i++) {
            set(i, spins[i]);
        }
    }
    
    void getSpins(int *spins) {
        for (int i = 0; i < n; i++) {
            spins[i] = get(i);
        }
    }
    
    void print() {
        for (int y = 0; y < len; y++) {
            for (int x = 0; x < len; x++) {
                printf("%d ", get(y*len+x));
            }
            printf("\n");
        }
        printf("\n");
    }
    
    void compare(Ising *that) {
        for (int i = 0; i < n; i++) {
            if (get(i) != that->get(i)) {
                printf("Lattices are not equal at %d\n", i);
                return;
            }
        }
        printf("Lattices are equal\n");
    }

protected:
    int shouldFlipSpin(int s, int m) {
        float dE = 2*s*(m + h);
        if (dE < 0)
            return 1;
        else {
            float r = 0.5; // rand() / (float)RAND_MAX;
            return exp(- dE / T) > r;
        }
    }
    
    int indexParity(int i) {
        int acc = 0;
        int len_d = 1;
        for (int d = 0; d < dim; d++) {
            int x = (i / len_d) % len;
            acc += x;
            len_d *= len;
        }
        return acc % 2;
    }
};


// ----------------------------------------------------------------------------
// Ising1 -- slow implementation
//

class Ising1 : virtual public Ising {
public:
    int *spins;
    
    Ising1(int len, int dim, float h, float T);
    ~Ising1();
    
    void set(int i, int s);
    int get(int i);
    void completeNeighborSum(int *sum);
    void update(int parityTarget);
    
private:
    int neighborSum(int i);
};


// ----------------------------------------------------------------------------
// Ising2 -- bits packed implementation
//

void initPentacubeParity();

class Ising2 : virtual public Ising {
public:
    unsigned int *blocks; // n/32 blocks represents n spins
    
    Ising2(int len, int dim, float h, float T);
    ~Ising2();
    
    void set(int i, int s);
    int get(int i);
    void completeNeighborSum(int *sum);
    void update(int parityTarget);
    
private:
    void index(int i, int *ip, int *delta);
    int reverseIndex(int ip, int delta);
};


// ----------------------------------------------------------------------------
// NVIsing
//

typedef unsigned int bits32;

typedef struct {
    // a lattice containing (n = len^dim) spins
    int len, dim, n;
    unsigned int *spins; // n spins stored using n/32 unsigned ints
    float h, T; // external field and temperature
} NVIsing;


NVIsing nvAllocate(int len, int dim, int n);
void    nvFree(NVIsing ising);
void    nvLoadSpins(NVIsing ising, unsigned int *spins);
void    nvRetrieveSpins(NVIsing ising, unsigned int *spins);
void    nvUpdate(NVIsing ising, int parityTarget);
