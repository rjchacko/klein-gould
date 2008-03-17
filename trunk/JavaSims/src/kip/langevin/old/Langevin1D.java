package kip.langevin.old;

import org.opensourcephysics.controls.SimControl;
import static java.lang.Math.*;
import kip.util.Random;

public class Langevin1D implements Cloneable {
    public int N;
    public double[] ψ, φ;
    public double L, dx, R, h, ε, t, dt, λ, alpha;
    public double ψcutoff;
    public double M = 1;
    public int randomSeed;
    
    Random random = new Random();
    double[] ψnew;
    
    
    public Object clone() {
        try {
            Langevin1D c = (Langevin1D)super.clone();
            c.ψ = (double[])ψ.clone();
            c.φ = (double[])φ.clone();
            c.ψnew = new double[N];
            c.random = (Random)random.clone();
            return c;
        } catch (Exception e) {
            return null;
        }
    }
    
    
    public void initialize(SimControl control) {
        randomSeed = control.getInt("Random seed");
        random.setSeed(randomSeed);        
        ψcutoff = control.getDouble("Crude cutoff");
        L = control.getDouble("Length/R");
        dx = control.getDouble("dx/R");
        getParameters(control);
        
        t = 0;
        N = (int) (L / dx);
        
        if (ψ == null || ψ.length != N) {
            ψ = new double[N];
            ψnew = new double[N];
            φ = new double[N];
        }
        else {
            for (int i = 0; i < N; i++)
                ψ[i] = φ[i] = 0;
        }
    }
    
    
    public void getParameters(SimControl control) {
        R = control.getDouble("R");
        h = control.getDouble("h");
        ε = control.getDouble("\u03b5");
        λ = control.getDouble("\u03bb");
        alpha = control.getDouble("\u03b1");
        dt = control.getDouble("dt");
    }


    //                                            ⌠t' -α(t-t')
    // ∂ψ/∂t = -M(-R²∇²ψ + 2εψ + 4ψ³ - h) + η + λ⎮ e         ψ(t') dt'
    //                                            ⌡t=0
    //
    // where the final term is identified with the field φ(t), and
	//
    //	      ⌠ t+dt
    //        |      η dt = sqrt(dt / dx R) guassian_noise()
    //        ⌡ t
    //
    public void step() {
        for (int i = 0; i < N; i++) {
            int ip = (i+1) % N;
            int im = (i-1+N) % N;
            double ψ3 = ψ[i] * ψ[i] * ψ[i];
            double laplace = (ψ[ip]-2*ψ[i]+ψ[im])/(dx*dx);
            double dψ_dt = -M * (-laplace + 2*ε*ψ[i] + 4*ψ3 - h);
            double dt_eta = sqrt(dt/(dx*R)) * noise();
            φ[i] += -alpha*φ[i]*dt + λ*ψ[i]*dt;
            ψnew[i] = ψ[i] + dt*dψ_dt + dt_eta + dt*φ[i];
        }
        for (int i = 0; i < N; i++) {
            ψ[i] = ψnew[i];
        }
        
        t += dt;
    }
    
    
    public double averageψ() {
        double acc = 0;
        for (int i = 0; i < N; i++)
            acc += ψ[i];
        return acc / N;
    }

    
    // if the system were mean field, then the grad2 would
    // disappear from free energy, and psi becomes just a number.
    // evaluate such a free energy for given psi
    private double evalMeanField(double psi) {
        return (2*ε-λ/alpha)*psi + 4*psi*psi*psi - h;
    }
    
    
    // assuming the field ψ is spatially constant, return
    // its value in metastable equilibrium
    public double metastableψ() {
        double hi = -sqrt(-(2*ε-λ/alpha)/12);
        double lo = -10;
        
        assert (0 < evalMeanField(hi)) : "Below the spinodal!";
        assert (0 > evalMeanField(lo));
        
        while (hi - lo > 1e-12) {
            double mid = (lo + hi) / 2;
            if (evalMeanField(mid) < 0)
                lo = mid;
            else
                hi = mid;
        }
        return (lo + hi) / 2;
    }
    
    
    // evaluate noise in η
    private double noise() {
        return random.nextGaussian();
        // return (random.nextDouble()*2-1)*sqrt(3);
    }
    
    
    // choose a new random number seed
    public void perturb() {
        random = new Random();
    }
    

    // has system nucleated?
    public boolean nucleated() {
        return ψ[maxIndex(ψ)] > ψcutoff;
    }
    
    
    // has the system nucleated within 'range' of position 'x'
    public boolean nucleatedInRegion(double x, double range) {
        int ilo = x2i(x-range/2);
        int ihi = x2i(x+range/2);
        for (int i = ilo; i != ihi; i = (i+1)%N) {
            if (ψ[i] > ψcutoff)
                return true;
        }
        return false;
    }
    
    
    // the 'x' position of droplet center
    public double dropletLocation() {
        int i = maxIndex(ψ);
        return (ψ[i] > ψcutoff) ? i*dx : Double.NaN;    
    }
    
    
    // given a guess x for the droplet center position, return
    // delta = (x' - x), the correction, where x' is a more careful
    // measurement.
    //
    // unfortunately, droplet profiles can be pretty weird, so
    // it's hard to say exactly where the center is.  this method
    // attempts to take an accurate measure by averaging all
    // points where the psi field crosses the nucleating threshhold.
    public double dropletDelta(double x, double range) {
        double acc = 0, cnt = 0;
        int ilo = x2i(x-range/2);
        int ihi = x2i(x+range/2);
        
        for (int i = ilo; i != ihi; i = (i+1)%N) {
            if (ψ[i] > ψcutoff) {
                acc += deltaX(i*dx, x);
                cnt++;
            }
        }
        return acc/cnt; // NaN if no droplet
    }
    
    public double deltaX(double x1, double x2) {
        double d = x1 - x2;
        if (d < -L/2) return d + L;
        if (d > L/2) return d - L;
        return d;
    }
    
    public double distance(double x1, double x2) {
        return abs(deltaX(x1, x2));
    }
    
    public int x2i(double x) {
        int i = (int)round(x / dx);
        if (i < 0) return i+N;
        if (i >= N) return i-N;
        return i;
    }
    
    public double meanψ() {
        double acc = 0;
        for (int i = 0; i < N; i++)
            acc += ψ[i];
        return acc / N;
    }
    
    
    public double sigma2ψ() {
        double meanψ = meanψ();
        double acc = 0;
        for (int i = 0; i < N; i++) {
            double dψ = ψ[i]-meanψ;
            acc += dψ*dψ;
        }
        return acc / N;
    }
    
    public void dump(String path, double[] _ψ) {
        String filename = path + "/_h" + h + "_r" + randomSeed + "_t" + (int)t + "_dt" + dt + "_R" + R;
        if (_ψ == null)
            _ψ = ψ;
        scikit.util.Dump.dumpColumns(filename, _ψ, 1);
    }
    
    private static int maxIndex(double[] a) {
        double a_max = a[0];
        int i_max = 0;
        for (int i = 1; i < a.length; i++) {
            if (a[i] > a_max) {
                a_max = a[i];
                i_max = i;
            }
        }
        return i_max;
    }
    
}

