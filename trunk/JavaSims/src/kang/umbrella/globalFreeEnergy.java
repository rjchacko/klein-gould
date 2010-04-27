package kang.umbrella;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.IOException;

import scikit.dataset.Accumulator;
import scikit.dataset.DatasetBuffer;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.numerics.fft.managed.ComplexDoubleFFT_Mixed;

// TODO: check all instances of data.copyData() - kip
public class globalFreeEnergy extends Simulation {
	Accumulator freeEnergy=new Accumulator();
	Accumulator fslope=new Accumulator();
	Accumulator fslope2=new Accumulator();
	Accumulator fftIAccum=new Accumulator();
	Accumulator fftRAccum=new Accumulator();
	Accumulator fftMagAccum=new Accumulator();
	Plot freeEnergyPlot=new Plot("Free Energy");
	Plot fslopePlot=new Plot("First Derivative of Free Energy");
	Plot fslope2Plot=new Plot("Second Derivative of Free Energy");
	Plot fftIPlot=new Plot("Fourier Transform I");
	Plot fftRPlot=new Plot("Fourier Transform R");
	int firstWindow,lastWindow;
	String filename;
	
	double fftdata[]=null;
	ComplexDoubleFFT_Mixed fft=null;
	
	public globalFreeEnergy() {
	}

	public void load(Control c){
		params.add("directory", "/Users/rjchacko/Desktop/NN2DIsingData/data64/Magnetization");
		params.add("first window",0);
		params.add("last window",1023);
		c.frame(freeEnergyPlot, fslopePlot, fslope2Plot, fftIPlot, fftRPlot);
	}
	@Override
	public void animate() {
		freeEnergyPlot.registerPoints("free energy", freeEnergy, Color.RED);
		fslopePlot.registerPoints("free energy slope", fslope, Color.RED);
		fslope2Plot.registerPoints("inverse susceptibility", fslope2, Color.RED);
		fftIPlot.registerPoints("fftI", fftIAccum, Color.RED);
		fftIPlot.registerPoints("fftR", fftRAccum, Color.BLUE);
//		fftPlot.registerPoints("fftMag", fftMagAccum, Color.BLUE);
		
	}

	@Override
	public void clear() {
		freeEnergy.clear();
		freeEnergyPlot.clear();
		fftIAccum.clear();
		fftIPlot.clear();
		fftRPlot.clear();
	}

	@Override
	public void run() {
		filename=params.sget("directory");
		firstWindow=params.iget("first window");
		lastWindow=params.iget("last window");
		Accumulator tempAccum=null;
		for(int i=firstWindow;i<=lastWindow;i++){
			try {
				tempAccum=readHistograms(i);
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			DatasetBuffer data=tempAccum.copyData();
			for(int j=0;j<data.size()-1;j++){	
				double slope=data.y(j+1)-data.y(j);
				fslope.accum((data.x(j+1)+data.x(j))/2.,slope);
			}
		}
		
		double integration=0;
		DatasetBuffer data=fslope.copyData();
		fftdata=new double[2*data.size()];
		fft=new ComplexDoubleFFT_Mixed(fftdata.length/2);
	    for(int i=0;i<data.size();i++){
	    	integration+=data.y(i);
	    	freeEnergy.accum(data.x(i),integration);
	    	fftdata[2*i]=integration;
	    }
	    
	    DatasetBuffer data2=fslope.copyData();
		for(int j=0;j<data2.size()-1;j++){	
			double slope=data2.y(j+1)-data2.y(j);
			fslope2.accum((data2.x(j+1)+data2.x(j))/2.,slope);
		}
	    
	    Job.animate();
	    
	    fft.transform(fftdata);
	 
	    for(int i=0;i<fftdata.length/2;i++){
	    	fftIAccum.accum(i, fftdata[2*i]);
	    	fftRAccum.accum(i,fftdata[2*i+1]);
//	    	fftMagAccum.accum(i, fftdata[2*i+1]*fftdata[2*i+1]+fftdata[2*i]*fftdata[2*i]);
	    }
	    Job.animate();
	}
	
	public Accumulator readHistograms(int i) throws FileNotFoundException, IOException {
		Accumulator temp=new Accumulator();
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
		return temp;
	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Control(new globalFreeEnergy(),"Free Energy");
	}

	
}
