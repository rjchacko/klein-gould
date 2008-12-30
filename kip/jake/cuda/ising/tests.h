#ifndef TESTS_H
#define TESTS_H

#include "ising.h"

void testSum (Ising &, Ising &);
void testUpdate (Ising &, Ising &);
void test1 (void);
void test2 (void);
double meanMagnetization (Ising &);
void test4 (void);
void test5 (void);
void test6 (void);
void test7 (void);

void stopwatchStart (void);
double stopwatchStop (void);

#endif // TESTS_H
