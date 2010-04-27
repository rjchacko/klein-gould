package kang.umbrella;

import java.awt.Color;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.HashMap;

import scikit.dataset.Accumulator;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;

public class multiHistogram extends Simulation {
	String filename;
	double T;
	int windowWidth,totalWindows, numberOfMCS, L;
	Accumulator x=new Accumulator(1);
	Plot xPlot=new Plot("Adjusted Probability");
	Accumulator totHist=new Accumulator(1);
	Plot totHistPlot=new Plot("Total Histogram");
	public multiHistogram() {
		// TODO Auto-generated constructor stub
	}
	
	public void load(Control c){
		params.add("T",1.0);
		params.add("L",32);
		params.add("number of windows", 256);
		params.add("MCS per window", 1000);
		params.add("window width",7);
		params.add("prefix","/Users/rjchacko/Desktop/data/Free Energy");
		params.add("iterations");
		c.frame(totHistPlot,xPlot);
	}
	
	
	public void animate() {
		xPlot.registerPoints("Adjusted P", x, Color.RED);
		totHistPlot.registerPoints("Total Histogram", totHist, Color.RED);
	}

	
	public void clear() {
		x.clear();
		xPlot.clear();
	}

	
	public void run() {
		HashMap<Double, Double> totalHistogram= new HashMap<Double, Double>();
		T=params.fget("T");
		totalWindows=params.iget("number of windows");
		windowWidth=params.iget("window width");
		numberOfMCS=params.iget("MCS per window");
		L=params.iget("L");
		filename=params.sget("prefix");
		int numberOfWindows=256;
		
//		int samples=numberOfMCS*L*L;
		try {
			makeTotalHistogram(numberOfWindows, totalHistogram);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		for(Double key: totalHistogram.keySet()){
			double value = totalHistogram.get(key);
			totHist.accum(key, value);
		}
		Job.animate();
		double[] z=new double[numberOfWindows];
		for(int i=0;i<z.length;i++){
			z[i]=1;
		}
		for(int k=0;k<50;k++){
			params.set("iterations", k);
			for(int j=0;j<1000;j++){
				z=updateZArray(z.clone(),totalHistogram);		
			}
					
			HashMap<Double, Double> p=adjustP(totalHistogram,z);

			for (Double key : p.keySet()) {
				double value = p.get(key);
				x.accum(key, T*Math.log(value));
			}
			Job.animate();
		}
	}

	public double biasFunction(int i, double phi){
		double x=0;
		double phi0=L*L-(i+1)*2*L*L/totalWindows;
		if(Math.abs(phi-phi0)<=windowWidth) x=1;
		else x=0;
		return x;
	}
	
	public void makeTotalHistogram(int numberOfWindows,HashMap<Double, Double> th) throws FileNotFoundException, IOException {
		for(int i=0;i<=numberOfWindows;i++){   					
			th=addToTotalHistogram(i, th);	 
		}

	}
	public HashMap<Double,Double> addToTotalHistogram(int i, HashMap<Double, Double> th) throws FileNotFoundException, IOException {
		String path=filename+i+".txt";
		java.io.BufferedReader reader=new java.io.BufferedReader(new java.io.FileReader(path));
		
		double numberOfOccurences;
		double binNumber;
		
		String s;
		while((s = reader.readLine())!=null) {
		      s = s.trim();
		      if(s.equals("")||(s.charAt(0)=='#')|| (s.charAt(0)=='=') ) { // ignore empty lines and lines beginning with #
		        continue;
		      }
		      
		      try {
		        String data[]=s.split(" ");
		        binNumber = Double.parseDouble(data[0]);
		        numberOfOccurences = Double.parseDouble(data[1]);
		        Double tempY= th.get(binNumber);
		        th.put(binNumber, (tempY==null)?numberOfOccurences:numberOfOccurences+tempY);
		      } catch(java.util.NoSuchElementException nsee) {
		        nsee.printStackTrace();
		      } catch(NumberFormatException nfe) {
		        nfe.printStackTrace();
		      }
		}
	
		return th;
	}
	
	public HashMap<Double, Double> adjustP(HashMap<Double, Double> th, double[] zz ){
		HashMap<Double, Double> p= new HashMap<Double, Double>();
		
		for(Double x: th.keySet()){	
			double denominator = calcDenom(zz, x);
			Double y=th.get(x);
			double p0=y/denominator;
			p.put(x, p0);
		}
		return p;
	}

	private double calcDenom(double[] zz, Double x) {
		double denominator=0;
		for(int j=0; j<zz.length; j++){
			denominator+=biasFunction(j,x)*zz[0]/zz[j];
		}
		return denominator;
	}
	
	public double computeZ(int i, double[] zz,HashMap<Double, Double> th){
		
		double integral=0;
		for(Double phi: th.keySet()){
			Double y=th.get(phi);
			double numerator=biasFunction(i,phi)*y;
			double denominator=0;
			for(int k=0;k<zz.length;k++){
				denominator+=biasFunction(k,phi)/zz[k];
			}
			integral+=numerator/denominator;
		}		
		return integral;
	}
	
	
	public double[]  updateZArray(double[] oldZ, HashMap<Double, Double> th){
		double[] newZ=new double[oldZ.length];
		newZ[0]=oldZ[0];
		for(int i=1;i<oldZ.length;i++){
			newZ[i]=computeZ(i,oldZ,th);
		}
		return newZ;
	}
	public static void main(String[] args) {
		new Control(new multiHistogram(),"Multiple Histogram Method");
	}

}
