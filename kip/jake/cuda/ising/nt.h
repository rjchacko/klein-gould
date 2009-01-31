#ifndef NT_H
#define NT_H

//#define MOUT
#define HLEN 25
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
    void nucleationTimes ();
    void lChange (int);
protected:
private:
    IsingCuda * ic; // Need ptr, otherwise constructor is called
    FILE * mOut, * ntOut;
    float headFlag; // For specifying beg/end data blocks in mOut
    float head[HLEN]; // 100B header

    void estMean_m (double &, double &);
};

#endif // NT_H
