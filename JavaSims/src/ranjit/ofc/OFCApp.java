package ranjit.ofc;


import java.awt.Color;
import java.io.IOException;

import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;

public class OFCApp extends Simulation{
	
    Lattice OFCModel;
    Site Initiator;
    int avalancheSize;
    Grid grid1 = new Grid ("Stress Lattice");
    Histogram avalancheSizes=new Histogram(1);
    Plot avSizes = new Plot("Avalanche Sizes");
    
	public static void main(String[] args) throws IOException {
		new Control(new OFCApp(),"OFC model");	    
	}

	public void load(Control c){
		params.add("L",256);
		params.add("R",20);
		params.add("alpha", 0.99);
		params.add("number of bins",20000);
		params.add("boundary conditions", new ChoiceValue("open", "closed"));
		params.add("number of steps", 10000);
		
		c.frame(grid1,avSizes);
	}
	@Override
	public void animate() {
		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);
		
		ColorGradient smooth = new ColorGradient();
		for (int i=0 ; i<OFCModel.sites.length ; i++){
			smooth.getColor(OFCModel.sites[i].stress,-2,OFCModel.maxStress);
		}
		
		grid1.setColors(smooth);
		double x[]=new double[OFCModel.sites.length];
		for(int i=0;i<x.length;i++){
			x[i]=OFCModel.sites[i].stress;
		}
		grid1.registerData(OFCModel.size,OFCModel.size,x);
		avSizes.registerBars("Avalanche Sizes", avalancheSizes, Color.RED);
		
	}

	@Override
	public void clear() {
		grid1.clear();
		avalancheSizes.clear();
	}

	@Override
	public void run() {
		
		int size = params.iget("L");
	    int range = params.iget("R");
	    int binnumber = params.iget("number of bins");	
	    double ALPHA= params.fget("alpha");
	    double ninR=(double) (2*range+1)*(2*range+1)-1;
	    double alpha= ALPHA*(1.0/ninR);
	    boolean boundary=false;	    
	    if(params.sget("boundary conditions")=="open") boundary= true;
	    
		OFCModel = new Lattice(size, range, binnumber, boundary, alpha); 
	    
	    OFCModel.initializeLattice();
	    Job.animate();
	    int maxSteps=params.iget("number of steps");
	    
		for(int i=0;i<maxSteps;i++){	    
	    	OFCModel.iteration=i;
	    	Initiator=OFCModel.findInitiator();
	    	OFCModel.maxStress=Initiator.stress;
	    	OFCModel.minStress=OFCModel.maxStress-1;
	    	OFCModel.maxBin=(int) (Math.floor(binnumber*OFCModel.maxStress)%binnumber);
	    	
	   	 	avalancheSize=OFCModel.Avalanche(Initiator);
	   	 	avalancheSizes.accum(avalancheSize);
	   	 	Job.animate();
		}
	   
		
//	    for(int i=0;i<size*size;i++){
//	    	System.out.println(OFCModel.sites[i].bin + "," + OFCModel.sites[i].failCounter);
//	    }
	}
    
}
