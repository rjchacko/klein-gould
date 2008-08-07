package chris.foo.ofc.apps.old;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.imageio.ImageIO;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.foo.ofc.old.SimpleDamage2D;


public class SimpleAnimationOFCApp extends Simulation{

	Grid grid1 = new Grid ("Stress Lattice");
	Grid grid2 = new Grid ("Failed Sites");
	Plot plot = new Plot("Histogram of Number of Showers");
	SimpleDamage2D model;
	String outdir;
	long counter;
	DecimalFormat fmt = new DecimalFormat("000000");

	public static void main(String[] args) {
		new Control(new SimpleAnimationOFCApp(), "OFC Model");
	}
	
	public void load(Control c) {
		params.add("Random Seed",0);
		params.add("Output directory", new DirectoryValue("/Users/cserino/CurrentSemester/Research/"));
		params.addm("Auto Scale", new ChoiceValue("Yes", "No"));
		params.add("Lattice Size",1<<9);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Stress Distribution", new ChoiceValue("Flat","Hammer Blow"));
		params.add("Noise", new ChoiceValue("On","Off"));
		params.add("Noise Width", new DoubleValue(0.05,0,0.1).withSlider());
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond"));
		params.addm("Interaction Radius (R)",(int)(50));		
		params.addm("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		c.frameTogether("OFC Model with Damage (Between Avalanches)", grid1, grid2);
		c.frame(plot);
		
	}
	
	
	
	public void animate() {

		int i;
		counter++;
		String tmp;
		
		tmp = fmt.format(counter);
		
		plot.setAutoScale(false);
		if (params.sget("Auto Scale").equals("Yes")) {
			plot.setAutoScale(true);
		}
		

		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);
		
		ColorGradient smooth = new ColorGradient();
		for (i=0 ; i<model.N ; i++){
			smooth.getColor(model.stress[i],-2,model.Sc);
		}
		
		grid1.setColors(smooth);
		grid2.setColors(palette);
		grid1.registerData(model.L,model.L,model.stress);
		grid2.registerData(model.L, model.L, model.alive);	
		
		params.set("Number of Resets",model.time);
		params.set("Number of Showers",model.showernumber);

		plot.registerBars("histNshowers", model.histNS, Color.RED);
		
		// Animation Stuff
		String imgName = outdir+File.separator+"OFC-"+tmp+".png";
		try {
			ImageIO.write(grid2.getImage(), "png", new File(imgName));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + imgName);
		}
	}

	public void clear() {
		
		plot.clear();
		grid1.clear();
		grid2.clear();
	
	}

	public void run() {
		counter=0;
		
		outdir = params.sget("Output directory");
		
		model = new SimpleDamage2D(params);
		
	
		model.Initialize(params.sget("Stress Distribution"));
		
		while(!(model.crack)){
			model.Avalanche();
			// calculate R_gyr here
			Job.animate();
		}
		Job.animate();
		
	}

	
	
	
}
