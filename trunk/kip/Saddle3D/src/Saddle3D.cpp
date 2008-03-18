
#include <complex>
#include <iostream>
#include <fstream>
#include <fftw3.h>
#include <math.h>
#include <time.h>


using namespace std;

const double DENSITY = 1;
const double KR_SP = 5.76345919689454979140645;
const double T_SP = 0.08617089416190739793014991;

const int Lp = 64;
const int Lp3 = Lp*Lp*Lp;
const double L = 16000; 
const double R = 1300;
const double dx = L/Lp;
const double dt = 0.5;
const double T = 0.095;  
const double packingFraction = 0.05;

double phi[Lp3];
double phi_bar[Lp3];
double del_phi[Lp3];
fftw_complex* phi_fft;
fftw_plan fftForward, fftBackward;
double rms_dF_dphi;
double freeEnergyDensity;
bool rescaleClipped;
double t = 0;


enum SeedType {BCC, TRIANGLE};

double sqr(double x) {
	return x*x;
}

void initializeFFT() {
    phi_fft = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Lp3);
    fftForward = fftw_plan_dft_3d(Lp, Lp, Lp, phi_fft, phi_fft, FFTW_FORWARD, FFTW_MEASURE);
    fftBackward = fftw_plan_dft_3d(Lp, Lp, Lp, phi_fft, phi_fft, FFTW_BACKWARD, FFTW_MEASURE);
}

void destroyFFT() {
    fftw_destroy_plan(fftForward);
    fftw_destroy_plan(fftBackward);
    fftw_free(phi_fft);
}

void initializeField(SeedType type) {
	for (int i = 0; i < Lp3; i++) {
		double x = dx*(i%Lp - Lp/2);
		double y = dx*((i%(Lp*Lp))/Lp - Lp/2);
		double z = dx*((i/(Lp*Lp)) - Lp/2);
		double field = 0;
		double k = KR_SP/R;
		if (type == BCC) {
			// BCC (reciprocal lattice is FCC)
			field += cos(k * ( x + z) / sqrt(2));
			field += cos(k * (-x + z) / sqrt(2));
			field += cos(k * ( y + z) / sqrt(2));
			field += cos(k * (-y + z) / sqrt(2));
		}
		
		double r = sqrt(x*x+y*y+z*z);
		double mag = 0.2 / (1+sqr(r/R));
		phi[i] = DENSITY*(1+mag*field);
	}
}

double arrayMax(double a[]) {
	double ret = a[0];
	for (int i = 1; i < Lp3; i++)
		ret = max(ret, a[i]);
	return ret;
}

double arrayMin(double a[]) {
	double ret = a[0];
	for (int i = 1; i < Lp3; i++)
		ret = min(ret, a[i]);
	return ret;
}

double mean(double a[]) {
	double acc = 0;
	for (int i = 0; i < Lp3; i++)
		acc += a[i];
	return acc / Lp3;
}

double meanSquared(double a[]) {
	double acc = 0;
	for (int i = 0; i < Lp3; i++)
		acc += a[i]*a[i];
	return acc / Lp3;
}

double variance(double a[]) {
	double m = mean(a);
	return meanSquared(a) - m*m;
}

double phiVariance() {
	return variance(phi);
}

void scaleField(double scale) {
	// phi will not be scaled above PHI_UB or below PHI_LB
	double PHI_UB = 20*DENSITY;
	double PHI_LB = 0.01*DENSITY;
	double s1 = (PHI_UB-DENSITY)/(arrayMax(phi)-DENSITY+1e-10);
	double s2 = (PHI_LB-DENSITY)/(arrayMin(phi)-DENSITY-1e-10);
	rescaleClipped = scale > min(s1,s2);
	if (rescaleClipped)
		scale = min(s1,s2);
	for (int i = 0; i < Lp*Lp*Lp; i++) {
		phi[i] = (phi[i]-DENSITY)*scale + DENSITY;
	}
}

double potential(double kR) {
	return (kR == 0) ? 1 : (3/(kR*kR))*(sin(kR)/kR - cos(kR));
}

void convolveWithPotential() {
	for (int i = 0; i < Lp3; i++) {
		phi_fft[i][0] = phi[i];
		phi_fft[i][1] = 0;
	}

	fftw_execute(fftForward);
    
	for (int x3 = -Lp/2; x3 < Lp/2; x3++) {
		for (int x2 = -Lp/2; x2 < Lp/2; x2++) {
			for (int x1 = -Lp/2; x1 < Lp/2; x1++) {
				int i = Lp*Lp*((x3+Lp)%Lp) + Lp*((x2+Lp)%Lp) + (x1+Lp)%Lp;
				double x = sqrt(x1*x1+x2*x2+x3*x3);
				double k = 2*M_PI*x/L;
				double J = potential(k*R);
				phi_fft[i][0] *= J;
				phi_fft[i][1] *= J;
			}
		}
	}
    
    fftw_execute(fftBackward);
 
    for (int i = 0; i < Lp3; i++)
		phi[i] = phi_fft[i][0] / Lp3;

//	fft.convolve(phi, phi_bar, new Function3D() {
//	public double eval(double k1, double k2, double k3) {
//		return potential(hypot(k1*Rx,k2*Ry,k3*Rz));
//	}
//	});	
}

double entropy(double phi) {
	if (packingFraction > 0) {
		double a = 1/packingFraction;
		return phi*log(phi) + (a-phi)*log(a-phi);
	}
	else {
		return phi*log(phi);
	}
}

double dentropy_dphi(double phi) {
	if (packingFraction > 0) {
		double a = 1/packingFraction;
		return log(phi) - log(a-phi);
	}
	else {
		return log(phi);
	}
}

double freeEnergyBackground() {
	return 0.5*sqr(DENSITY) + T*entropy(DENSITY);
}

void simulate() {
	convolveWithPotential();

	for (int i = 0; i < Lp3; i++) {
		del_phi[i] = - dt*(phi_bar[i]+T*dentropy_dphi(phi[i]));
	}
	double mu = mean(del_phi)-(DENSITY-mean(phi));
	for (int i = 0; i < Lp3; i++) {
		// clip del_phi to ensure phi(t+dt) > phi(t)/2
		del_phi[i] = max(del_phi[i]-mu, -phi[i]/2.);
	}

	rms_dF_dphi = 0;
	freeEnergyDensity = 0;
	for (int i = 0; i < Lp3; i++) {
		rms_dF_dphi += sqr(del_phi[i] / dt);
		freeEnergyDensity += 0.5*phi[i]*phi_bar[i]+T*entropy(phi[i]);
		phi[i] += del_phi[i];
	}
	rms_dF_dphi = sqrt(rms_dF_dphi/Lp3);
	freeEnergyDensity /= Lp3;
	freeEnergyDensity -= freeEnergyBackground();
	t += dt;
}


void writeArray(char *fname, double *a, int w, int h, int d) {
	 ofstream file;
	 file.open (fname, ios::out | ios::binary);
	 file.write((char*)&w, sizeof(int));
	 file.write((char*)&h, sizeof(int));
	 file.write((char*)&d, sizeof(int));
	 file.write((char*)a, (w*h*d)*sizeof(double));
	 file.close();
}

double time() {
	return ((double)clock())/CLOCKS_PER_SEC;
}

int main() {
	cout << "initialize\n";
	initializeFFT();
	initializeField(BCC);
	convolveWithPotential();
	
	for (int i = 0; i < 10; i++) {
		double t1 = time();
		simulate();
		double t2 = time();
		cout << (t2-t1) << endl;
	}
//	writeArray("/Users/kbarros/Desktop/test", phi, Lp, Lp, Lp);
	return 0;
}
