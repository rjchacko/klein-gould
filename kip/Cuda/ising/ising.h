
// ----------------------------------------------------------------------------
// Ising -- base class
//

class Ising {
public:
    int len, dim, n;
    float h, T;
    
    Ising(int len, int dim, float h, float T);
    virtual ~Ising() {}
    
    virtual void completeNeighborSum(int *sum) = 0;
    virtual void update(int parityTarget) = 0;
    
    void randomizeSpins();
    void setSpins(int *spins);
    void getSpins(int *spins);
    void printLattice(int *a);
    void print();
    void compare(Ising *that);

protected:
    virtual void transferHostToDevice() {}
    virtual void transferDeviceToHost() {}
    
    virtual void set(int i, int s) = 0;
    virtual int get(int i) = 0;

    int shouldFlipSpin(int s, int m);
    int indexParity(int i);
};


// ----------------------------------------------------------------------------
// Ising1 -- slow implementation
//

class Ising1 : virtual public Ising {
public:
    int *spins;
    
    Ising1(int len, int dim, float h, float T);
    ~Ising1();
    
    void completeNeighborSum(int *sum);
    void update(int parityTarget);

protected:
    void transferHostToDevice() {}
    void transferDeviceToHost() {}
    void set(int i, int s);
    int get(int i);

private:
    int neighborSum(int i);
};


// ----------------------------------------------------------------------------
// Ising2 -- bits packed implementation
//

class Ising2 : virtual public Ising {
public:
    int nblocks; // nblocks = n/32 when dim >= 5
    unsigned int *blocks;
    
    Ising2(int len, int dim, float h, float T);
    ~Ising2();
    
    void completeNeighborSum(int *sum);
    void update(int parityTarget);

protected:
    void set(int i, int s);
    int get(int i);

private:
    void index(int i, int *ip, int *delta);
    int reverseIndex(int ip, int delta);
};


// ----------------------------------------------------------------------------
// IsingCuda -- bits packed implementation on Cuda
//

class IsingCuda : virtual public Ising {
public:
    IsingCuda(int len, int dim, float h, float T);
    ~IsingCuda();
    
    void completeNeighborSum(int *sum);
    void update(int parityTarget);

protected:
    int nblocks; // nblocks = n/32 when dim >= 5
    unsigned int *blocks;
    unsigned int *d_blocks;
    
    void transferHostToDevice();
    void transferDeviceToHost();
    void set(int i, int s);
    int get(int i);

private:
    void index(int i, int *ip, int *delta);
};

