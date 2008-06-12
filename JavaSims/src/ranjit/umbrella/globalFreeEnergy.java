package ranjit.umbrella;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.IOException;

import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class globalFreeEnergy extends Simulation {
	Accumulator temp=new Accumulator();
	Accumulator freeEnergy1=new Accumulator();
	Accumulator freeEnergy2=new Accumulator();
	Accumulator freeEnergy3=new Accumulator();
	Accumulator freeEnergy4=new Accumulator();
	Accumulator fslope=new Accumulator();
	Plot freeEnergyPlot=new Plot("Free Energy");
	int firstWindow,lastWindow;
	String filename;
	public globalFreeEnergy() {
	}

	public void load(Control c){
		params.add("directory1", "/Users/rjchacko/Desktop/L256h0.4w32k/Free Energy");
		params.add("directory2", "/Users/rjchacko/Desktop/L256h0.45w32k/Free Energy");
		params.add("directory3", "/Users/rjchacko/Desktop/L256h0.5w32k/Free Energy");
		params.add("directory4", "/Users/rjchacko/Desktop/L256h0.55w32k/Free Energy");
		params.add("first window",1);
		params.add("last window",180);
		params.add("current window");
		c.frame(freeEnergyPlot);
	}
	@Override
	public void animate() {
		freeEnergyPlot.registerPoints("free energy1", freeEnergy1, Color.RED);
		freeEnergyPlot.registerPoints("free energy2", freeEnergy2, Color.BLUE);
		freeEnergyPlot.registerPoints("free energy3", freeEnergy3, Color.GREEN);
		freeEnergyPlot.registerPoints("free energy4", freeEnergy4, Color.YELLOW);
	}

	@Override
	public void clear() {
		freeEnergy1.clear();
		freeEnergyPlot.clear();
	}

	@Override
	public void run() {
		filename=params.sget("directory1");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			double data[]=temp.copyData();
			for(int j=0;j<data.length/2-1;j++){	
				double slope=data[2*j+3]-data[2*j+1];
				fslope.accum((data[2*j+2]+data[2*j])/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration1=0;
		double data1[]=fslope.copyData();
	    for(int i=0;i<data1.length/2;i++){
	    	integration1+=data1[2*i+1];
	    	freeEnergy1.accum(data1[2*i],integration1);
	    }
	    Job.animate();
	    filename=params.sget("directory2");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			double data[]=temp.copyData();
			for(int j=0;j<data.length/2-1;j++){	
				double slope=data[2*j+3]-data[2*j+1];
				fslope.accum((data[2*j+2]+data[2*j])/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration2=0;
		double data2[]=fslope.copyData();
	    for(int i=0;i<data2.length/2;i++){
	    	integration2+=data2[2*i+1];
	    	freeEnergy2.accum(data2[2*i],integration2);
	    }
	    Job.animate();
	    
	    filename=params.sget("directory3");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			double data[]=temp.copyData();
			for(int j=0;j<data.length/2-1;j++){	
				double slope=data[2*j+3]-data[2*j+1];
				fslope.accum((data[2*j+2]+data[2*j])/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration3=0;
		double data3[]=fslope.copyData();
	    for(int i=0;i<data3.length/2;i++){
	    	integration3+=data3[2*i+1];
	    	freeEnergy3.accum(data3[2*i],integration3);
	    }
	    Job.animate();
	    
	    filename=params.sget("directory4");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			double data[]=temp.copyData();
			for(int j=0;j<data.length/2-1;j++){	
				double slope=data[2*j+3]-data[2*j+1];
				fslope.accum((data[2*j+2]+data[2*j])/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration4=0;
		double data4[]=fslope.copyData();
	    for(int i=0;i<data4.length/2;i++){
	    	integration4+=data4[2*i+1];
	    	freeEnergy4.accum(data4[2*i],integration4);
	    }
	    Job.animate();
	}
	
	public void readHistograms(int i) throws FileNotFoundException, IOException {
		temp=new Accumulator();
		java.io.BufferedReader reader;
		String s;
		double numberOfOccurences;
		double binNumber;
		reader = new java.io.BufferedReader(new java.io.FileReader(filename+i+".txt"));
		while((s = reader.readLine())!=null) {
		      s = s.trim();
		      if(s.equals("")||(s.charAt(0)=='#')|| (s.charAt(0)=='=') ) { // ignore empty lines and lines beginning with #
		        continue;
		      }
		      try {
		        java.util.StringTokenizer st = new java.util.StringTokenizer(s, " ");
		        binNumber = Double.parseDouble(st.nextToken());
		        numberOfOccurences = Double.parseDouble(st.nextToken());
		        temp.accum(binNumber, numberOfOccurences);
		      } catch(java.util.NoSuchElementException nsee) {
		        nsee.printStackTrace();
		      } catch(NumberFormatException nfe) {
		        nfe.printStackTrace();
		      }
		}
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Control(new globalFreeEnergy(),"Free Energy");
	}

	
}
