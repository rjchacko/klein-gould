#ifndef NT_H
#define NT_H

// Leave this in if you want m output. WARNING: files could be large!
#define MOUT

#define HEADFLAG  10000000.0f
#define FOOTFLAG -10000000.0f
#define HLEN 10 // 40 byte header
#define FLEN 1  //  4 byte footer
/**
 * HEADER PROTOCOL
 * 0 - flag
 * 1 - l
 * 2 - d
 * 3 - h
 * 4 - T
 * end - flag
 *
 * (the rest are auxiliary)
 */

class nt
{
public:
    nt (int, int, double, double);
    ~nt ();
    void sim ();
    void sim (int);
    void lChange (int);
    void hChange (double);
protected:
private:
    IsingCuda * ic; // Need ptr, otherwise constructor is called
    FILE * mOut, * ntOut;
    float foot[FLEN];
    float head[HLEN]; // 100B header

    void estMean_m (double &, double &);
};

#endif // NT_H
