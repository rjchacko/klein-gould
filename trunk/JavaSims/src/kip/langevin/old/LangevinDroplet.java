package kip.langevin.old;

import org.opensourcephysics.controls.SimControl;
import org.opensourcephysics.frames.*;
import static java.lang.Math.*;


public class LangevinDroplet {
    PlotFrame dropletPlot = new PlotFrame("x", "droplet fields", "Critical Droplet");

    SimControl control;
    
    double dx;
    int N;
    
    // cumulative droplet fields
    double[] ψdroplet, φdroplet;
    int dropletCount;
    
    // distance at which two nucleating droplets are considered to be from the same place
    // (scaled by R)
    double dropletRange = 10;
    
    // time duration to wait for nucleation
    double t_wait;
    
    
    public void initialize(SimControl control) {
        dropletPlot.setMarkerColor(0, java.awt.Color.BLACK);
        dropletPlot.setMarkerColor(1, java.awt.Color.BLUE);
        dropletPlot.getDataset(2).setLineColor(java.awt.Color.RED);
        dropletPlot.getDataset(3).setLineColor(java.awt.Color.GREEN);
        dropletPlot.setMarkerShape(2, org.opensourcephysics.display.Dataset.PIXEL);
        dropletPlot.setMarkerShape(3, org.opensourcephysics.display.Dataset.PIXEL);
        dropletPlot.setConnected(2, true);
        dropletPlot.setConnected(3, true);
        // invisible datasets 4 and 5 store error bars
        dropletPlot.setMarkerShape(4, org.opensourcephysics.display.Dataset.NO_MARKER);
        dropletPlot.setMarkerShape(5, org.opensourcephysics.display.Dataset.NO_MARKER);
        
        this.control = control;
        double overshoot = control.getDouble("Intervention overshoot");
        double L = control.getDouble("Length/R");
        dx = control.getDouble("dx/R");
        N = (int) (L / dx);
        t_wait = 2*overshoot;
        ψdroplet = new double[N];
        φdroplet = new double[N];
        dropletCount = 0;
    }


    void plotDroplet() {
        dropletPlot.clearData();
        int di = N < 400 ? 1 : (N / 400);
        for (int i = 0; i < N; i += di) {
            if (dropletCount > 0) {
                dropletPlot.append(2, i*dx, getψ(i));
                dropletPlot.append(3, i*dx, getφ(i));
            }
        }
        dropletPlot.setMessage("count = " + dropletCount);
        dropletPlot.repaint();
    }
    
    
    void plotTestDroplet(Langevin1D sim, double dropLoc) {
        plotDroplet();
        int j = sim.x2i(dropLoc);
        int di = N < 400 ? 1 : (N / 400);
        for (int i = 0; i < N; i += di) {
            dropletPlot.append(0, i*dx, sim.ψ[(j+i+sim.N/2)%sim.N]);
            dropletPlot.append(1, i*dx, sim.φ[(j+i+sim.N/2)%sim.N]);
        }
        dropletPlot.setMessage("t = " + sim.t);
        dropletPlot.repaint();        
    }
    
    
    public void findDroplet(Langevin1D oldsim, double t_nuc, double dropLoc) {
        // set up time bounds for intervention search
        double timerange = 1;
        Langevin1D sim;
        double[] a;
        double t_lo  = max(t_nuc - 0.8*(t_wait/2), 0);
        double t_hi  = t_nuc + 0.8*(t_wait/2);
        control.println("["+(int)t_lo + ", " + (int)t_hi + "]");
        //        visualIntervention(t_lo, t_hi, t_wait, control.getDouble("Intervention dt"), dropLoc);
        
        // binary search for nucleation time
        while (true) {
            if (t_hi - t_lo < timerange/4) {
                t_lo -= timerange;
                t_hi += timerange;
            }
            t_nuc = (t_lo + t_hi) / 2;
            sim = averagedSimulation(oldsim, t_nuc, timerange);
            plotTestDroplet(sim, dropLoc);
            a = smartNucleationProbabilityAndLocation(sim, dropLoc);
            
            if (a[0] < 0.45)
                t_lo = t_nuc;
            else if (a[0] > 0.55)
                t_hi = t_nuc;
            else
                break;
        }

        // accumulate result into nucleating droplet
        System.out.println("Seed: " + sim.randomSeed + " Time: " + sim.t);
        control.println("Old/new loc " + dropLoc + " / " + a[1]);
        
        double[] ψ = new double[N];
        double[] φ = new double[N];
        int j = sim.x2i(a[1]);
        for (int i = 0; i < sim.N; i++) {
            ψ[i] = sim.ψ[(j+i+sim.N/2)%sim.N];
            φ[i] = sim.φ[(j+i+sim.N/2)%sim.N];
            ψdroplet[i] += ψ[i];
            φdroplet[i] += φ[i];
        }
        dropletCount++;
        control.println();
        
        String path = control.getString("Data path");
        if (path.length() > 0)
            sim.dump(path, ψ);
    }
    
    
    Langevin1D averagedSimulation(Langevin1D oldsim, double t, double timerange) {
        assert oldsim.t <= t-timerange/2;        

        Langevin1D sim = (Langevin1D)oldsim.clone();
        double[] ψavg = new double[sim.N];
        double[] φavg = new double[sim.N];
        int cnt = 0;
        
        while (sim.t < t-timerange/2)
            sim.step();
        while (sim.t < t+timerange/2) {
            for (int i = 0; i < sim.N; i++) {
                ψavg[i] += sim.ψ[i];
                φavg[i] += sim.φ[i];
            }
            cnt++;
            sim.step();
        }
        for (int i = 0; i < sim.N; i++) {
            sim.ψ[i] = ψavg[i] / cnt;
            sim.φ[i] = φavg[i] / cnt;
        }        
        sim.t = t;
        return sim;
    }
    
    
    // returns {probability, new drop location}
    double[] nucleationProbabilityAndLocation(Langevin1D sim, double dropLoc, int trials) {
        double count = 0;
        double dropDelta = 0;
        for (int i = 0; i < trials; i++) {
            Langevin1D sim2 = (Langevin1D)sim.clone();
            sim2.perturb();
            while (sim2.t - sim.t < t_wait) {
                sim2.step();
                if (sim2.nucleatedInRegion(dropLoc, dropletRange)) {
                    count++;
                    dropDelta += sim2.dropletDelta(dropLoc, dropletRange);
                    break;
                }
            }
        }
        double delta = count == 0 ? 0 : dropDelta / count;
        if (abs(delta / sim.L) > 0.1) {
            System.out.println("yow! delta = " + delta + " seed = " + sim.randomSeed);
        }
        return new double[] {count/trials, dropLoc + delta};
    }
    
    
    double[] smartNucleationProbabilityAndLocation(Langevin1D sim, double dropLoc) {
        control.println("t = " + sim.t);

//        double[] cutoff = new double[] {0.41, 0.31};
//        int[] trials = new int[]         {10,  20};
        double[] cutoff = new double[] {0.41, 0.31, 0.21};
        int[] trials = new int[]         {10,  20,    40};
        
        double a[] = nucleationProbabilityAndLocation(sim, dropLoc, trials[0]);
        double p = a[0];
        control.println("p = (     " + trials[0] + " * " + p + ") -> " + p);
        
        for (int i = 0; i < cutoff.length; i++) {
            if (abs(p-0.5) > cutoff[i])
                break;
            a = nucleationProbabilityAndLocation(sim, dropLoc, trials[i]);
            p = (p + a[0]) / 2;
            control.println("         ( + " + trials[i] + " * " + a[0] + ") -> " + p);
        }
        a[0] = p;
        return a;
    }

/*    
    PlotFrame interventionPlot = new PlotFrame("t", "prob. nuc.", "Intervention");
    void visualIntervention(double t_lo, double t_hi, double t_wait, double dt, double dropLoc) {
        interventionPlot.setPreferredMinMax(t_lo, t_hi, 0, 1);
        int N = (int)round((t_hi - t_lo) / dt);
        double[] p = new double[N];
        
        for (int count = 1; animationThread == Thread.currentThread(); count++) {
            control.println("["+(int)t_lo + ", " + (int)t_hi + "], t = " + sim.t);        
            Langevin1D sim = newSimulationAtTime(t_lo);
            for (int i = 0; i < N; i++) {
                while (sim.t < t_lo + i*dt)
                    sim.step();
                p[i] += nucleationProbabilityAndLocation(sim, t_wait, dropLoc, 10)[0];
            }
            interventionPlot.clearData();
            for (int i = 0; i < N; i++)
                interventionPlot.append(0, t_lo + i*dt, p[i]/count);
            interventionPlot.setMessage("averages : " + count*10);
            interventionPlot.repaint();
            control.clearMessages();
            Thread.yield();
        }
    }
*/
    
    public double getψ(int i) {
        return ψdroplet[i]/dropletCount;
    }
    
    public double getφ(int i) {
        return φdroplet[i]/dropletCount;
    }
}