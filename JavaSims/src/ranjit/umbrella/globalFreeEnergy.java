package ranjit.umbrella;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.IOException;

import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

// TODO: check all instances of data.copyData() - kip
public class globalFreeEnergy extends Simulation {
	Accumulator temp=new Accumulator();
	Accumulator freeEnergy1=new Accumulator();
	Accumulator freeEnergy2=new Accumulator();
	Accumulator freeEnergy3=new Accumulator();
	Accumulator freeEnergy4=new Accumulator();
	Accumulator freeEnergy5=new Accumulator();
	Accumulator fslope=new Accumulator();
	Plot freeEnergyPlot=new Plot("Free Energy");
	int firstWindow,lastWindow;
	String filename;
	public globalFreeEnergy() {
	}

	public void load(Control c){
//		params.add("directory1", "/Users/rjchacko/Desktop/L32h0.4w512/Free Energy");
//		params.add("directory2", "/Users/rjchacko/Desktop/L32h0.45w512/Free Energy");
//		params.add("directory3", "/Users/rjchacko/Desktop/L32h0.5w512/Free Energy");
//		params.add("directory4", "/Users/rjchacko/Desktop/L32h0.55w512/Free Energy");
//		params.add("directory5", "/Users/rjchacko/Desktop/L32h0.6w512/Free Energy");
		params.add("directory1", "/Users/rjchacko/Desktop/eq initial conditions/L256h0.4w32kr1/Free Energy");
		params.add("directory2", "/Users/rjchacko/Desktop/eq initial conditions/L256h0.4w32kr1/Free Energy");
		params.add("directory3", "/Users/rjchacko/Desktop/eq initial conditions/L256h0.4w32kr1/Free Energy");
		params.add("directory4", "/Users/rjchacko/Desktop/eq initial conditions/L256h0.4w32kr1/Free Energy");
		params.add("directory5", "/Users/rjchacko/Desktop/eq initial conditions/L256h0.4w32kr1/Free Energy");
		params.add("first window",0);
		params.add("last window",100);
		params.add("current window");
		c.frame(freeEnergyPlot);
	}
	@Override
	public void animate() {
		freeEnergyPlot.registerPoints("free energy1", freeEnergy1, Color.RED);
		freeEnergyPlot.registerPoints("free energy2", freeEnergy2, Color.BLUE);
		freeEnergyPlot.registerPoints("free energy3", freeEnergy3, Color.GREEN);
		freeEnergyPlot.registerPoints("free energy4", freeEnergy4, Color.RED);
		freeEnergyPlot.registerPoints("free energy5", freeEnergy5, Color.BLUE);
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
		
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=temp.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration1=0;
		DatasetBuffer data1=fslope.copyData();
	    for(int i=0;i<data1.size();i++){
	    	integration1+=data1.y(i);
	    	freeEnergy1.accum(data1.x(i),integration1);
	    }
	    Job.animate();
	    fslope.clear();
	    
	    filename=params.sget("directory2");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=temp.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration2=0;
		DatasetBuffer data2=fslope.copyData();
	    for(int i=0;i<data2.size();i++){
	    	integration2+=data2.y(i);
	    	freeEnergy2.accum(data2.x(i),integration2);
	    }
	    Job.animate();
	    
	    fslope.clear();
	    filename=params.sget("directory3");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=temp.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration3=0;
		DatasetBuffer data3=fslope.copyData();
	    for(int i=0;i<data3.size();i++){
	    	integration3+=data3.y(i);
	    	freeEnergy3.accum(data3.x(i),integration3);
	    }
	    Job.animate();
	    
	    fslope.clear();
	    filename=params.sget("directory4");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=temp.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration4=0;
		DatasetBuffer data4=fslope.copyData();
	    for(int i=0;i<data4.size();i++){
	    	integration4+=data4.y(i);
	    	freeEnergy4.accum(data4.x(i),integration4);
	    }
	    Job.animate();
	    
	    fslope.clear();
	    filename=params.sget("directory5");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=temp.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
			
			params.set("current window",i);
		}
		double integration5=0;
		DatasetBuffer data5=fslope.copyData();
	    for(int i=0;i<data5.size();i++){
	    	integration5+=data5.y(i);
	    	freeEnergy5.accum(data5.x(i),integration5);
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
