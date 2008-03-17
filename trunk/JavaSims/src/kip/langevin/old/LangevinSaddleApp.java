package kip.langevin.old;

import org.opensourcephysics.controls.*;
import org.opensourcephysics.frames.*;
import static java.lang.Math.*;

public class LangevinSaddleApp extends AbstractSimulation {
   PlotFrame fieldPlot = new PlotFrame("x", "fields", "Langevin");
   Langevin1D sim = new Langevin1D();
   
   int stepsPerDisplay;
   double metaψ;
   double correction = -1.5;
   
   void plotFields() {
      fieldPlot.clearData();
      int di = sim.N < 400 ? 1 : (sim.N / 400);
      for (int i = 0; i < sim.N; i += di) {
         fieldPlot.append(0, i*sim.dx, sim.ψ[i]);
         fieldPlot.append(1, i*sim.dx, sim.φ[i]);
      }
      // fieldPlot.setMessage("t = " + (int)sim.t);      
   }
   
   
   public void initialize() {
      sim.initialize(control);
      metaψ = sim.metastableψ();
      
      for (int i = 0; i < sim.N; i++)
         sim.ψ[i] = sim.φ[i] = metaψ;
      sim.ψ[sim.N/2] += 0.1;
      
      fieldPlot.setPreferredMinMax(0, sim.L, -1, 1);      
	  
	  control.println("Metastable field: " + metaψ);
   }
   
   
   public void doStep() {
      for (int i = 0; i < stepsPerDisplay; i++) {
         double avg1 = sim.averageψ();
         sim.step();
         double avg2 = sim.averageψ();
                  
         if (correction != 0) {
            double scale = (correction*(avg2-avg1) + avg2 - metaψ) / (avg2 - metaψ);
            
            for (int j = 0; j < sim.N; j++) {
               sim.ψ[j] = scale * (sim.ψ[j] - metaψ) + metaψ;
               sim.ψ[j] = min(max(sim.ψ[j], metaψ), -metaψ);
            }
            
            sim.ψ[sim.N/2] = max(sim.ψ[sim.N/2], metaψ+0.005);
            sim.ψ[0] = sim.ψ[sim.N-1] = metaψ;

            fieldPlot.setMessage("scale "+scale);
         }
      }
      
      plotFields();
   }
   
   
   public void startRunning() {
      sim.getParameters(control);
      stepsPerDisplay = control.getInt("Steps per display");
   }
   
   
   public void reset() {
      control.setAdjustableValue("Steps per display", 2000);
	  control.setAdjustableValue("\u03bb", 0); // λ
	  control.setAdjustableValue("h", 0.223);
      control.setValue("Length", 50);
      control.setValue("dx", 0.5);
      control.setValue("Random seed", 0);
      control.setValue("Crude cutoff", 0);
      control.setAdjustableValue("dt", 0.01);
	  control.setAdjustableValue("R\u00b2", 0.1); // R2
	  control.setAdjustableValue("\u03b5", "-5/9"); // ε
      control.setAdjustableValue("\u03b1", 1); // α
      control.setAdjustableValue("\u0393", 0.0); // Γ
   }
   
   public void freeDynamics() {
      correction = (correction == 0) ? -1.5 : 0;
   }
   
   public static void main(String[] args) {
      SimulationControl c = SimulationControl.createApp(new LangevinSaddleApp());
      c.addButton("freeDynamics", "Free");
   }
}
