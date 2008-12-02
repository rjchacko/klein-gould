#include <stdio.h>
#include "ising.h"


void accumulateStatistics(Ising *ising, double *m_acc, double *msqr_acc) {
    double m = ising->magnetization() / ising->n;
    *m_acc += m;
    *msqr_acc += m*m;
}

// find_tc() assumes initCuda() already called
void find_tc() {
    int l = 128;
    int d = 2;
    double h = 0.001;
    double Tmin = 2.2;
    double Tmax = 2.4;
    double dT = 0.005;
    
    Ising1 *ic1 = new Ising1(l, d, h, Tmin);
    Ising2 *ic2 = new Ising2(l, d, h, Tmin);
    IsingCuda *ic3 = new IsingCuda(l, d, h, Tmin);

    int relaxTime = 1000;
    int iters = 100000;
    
    printf("# Comparing Ising1, Ising2, IsingCuda.\n");
    printf("# d=%d, L=%d, iters=%d, h=%f\n", d, l, iters, h);
    printf("# format: T <m_1> <m_2> <m_3> <chi^2_1> <chi^2_2> <chi^2_3>\n");
    
    for (double T=Tmin; T<Tmax; T+=dT) {
        ic1->T = T;
        ic2->T = T;
        ic3->T = T;
        ic1->setSpinsUp();
        ic2->setSpinsUp();
        ic3->setSpinsUp();
        
        for (int i=0; i<relaxTime; ++i) {
            ic1->update (0); ic1->update (1);
            ic2->update (0); ic2->update (1);
            ic3->update (0); ic3->update (1);
        }
        
        double m_1=0, msqr_1=0;
        double m_2=0, msqr_2=0;
        double m_3=0, msqr_3=0;
        
        for (int i=0; i<iters; ++i) {
            ic1->update (0); ic1->update (1);
            ic2->update (0); ic2->update (1);
            ic3->update (0); ic3->update (1);
            
            accumulateStatistics(ic1, &m_1, &msqr_1);
            accumulateStatistics(ic2, &m_2, &msqr_2);
            accumulateStatistics(ic3, &m_3, &msqr_3);
        }
        m_1 /= iters;
        m_2 /= iters;
        m_3 /= iters;
        double chi2_1 = (msqr_1 / iters) - m_1*m_1;
        double chi2_2 = (msqr_2 / iters) - m_2*m_2;
        double chi2_3 = (msqr_3 / iters) - m_3*m_3;
        
//        printf ("%f %f %f %f %f\n", T, m_1, m_2, chi2_1, chi2_2);
        printf ("%f %f %f %f %f %f %f\n", T, m_1, m_2, m_3, chi2_1, chi2_2, chi2_3);
    }
    printf ("\n");
}
