package chris.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

public class FitUtil{
	
	private double x, xs, y, xy, Xdata[], Ydata[], sigma, m, b, dm, db, chis;
	private int Ndata;
	private static int dlength = 100000;

	
	public FitUtil(){
		
		resetFitter();
		return;
	}
	
	public void resetFitter(){
		
		m    = 0;
		b    = 0;
		dm   = 0;
		db   = 0;
		chis = 0;
		return;
	}
	
	
	public void fit(String fin, int skip, int xid, int yid, double err){

		int c1, c2, counter;
		boolean order;
		Xdata   = new double[dlength];
		Ydata   = new double[dlength];
		counter = 0;
		x       = 0;
		xs      = 0;
		y       = 0;
		xy      = 0;
		Ndata   = 0;
		sigma   = err;
		
		if(xid < yid){
			order = true;
			c1    = xid;
			c2    = yid;
		}
		else{
			order = false;
			c1    = yid;
			c2    = xid;
		}
	
		try {
			
			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));
			
			String rin;
			int pd;
			
			for (int jj = 0 ; jj < skip ; jj++){
				rin = bir.readLine();
			}
			
			while ( (rin = bir.readLine()) != null ){
				pd = rin.indexOf('\t');
				for(int jj = 1 ; jj < c1 ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(order){
					if(pd == -1){
						Xdata[counter] = Double.parseDouble(rin);
					}
					else{
						Xdata[counter] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				else{
					if(pd == -1){
						Ydata[counter] = Double.parseDouble(rin);
					}
					else{
						Ydata[counter] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				for(int jj = c1 ; jj < c2 ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(order){
					if(pd == -1){
						Ydata[counter++] = Double.parseDouble(rin);
					}
					else{
						Ydata[counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				else{
					if(pd == -1){
						Xdata[counter++] = Double.parseDouble(rin);
					}
					else{
						Xdata[counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				if(counter%dlength == 0){
					calcFitVals(counter);
					counter = 0;
				}
				Ndata++;
			}
			calcFitVals(counter);
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		
		getFitVals();
		chis = sigma;
		if(sigma > 0) calcChis(fin, skip, xid, yid);

		return;
	}

	
	public void fit(double[] XX, double[] YY, double err){
		
		int counter;
		
		counter = 0;
		x       = 0;
		xs      = 0;
		y       = 0;
		xy      = 0;
		Ndata   = 0;
		sigma   = err;
		
		if(XX.length != YY.length){
			// ERROR!!!
			return;
		}
		Ndata = XX.length;
		if(Ndata > dlength){
			Xdata = new double[dlength];
			Ydata = new double[dlength];
			for(int jj = 0 ; jj < Ndata ; jj++){
				Xdata[counter] = XX[jj];
				Ydata[counter++] = YY[jj];
				if(counter%dlength == 0){
					calcFitVals(dlength);
					counter = 0;
				}
			}
			calcFitVals(counter);
		}
		else{
			Xdata = XX;
			Ydata = YY;
			calcFitVals(Ndata);
		}
		
		getFitVals();
		chis = sigma;
		if(sigma > 0) calcChis(XX, YY);
		
		return;
	}
	
			
	private void calcFitVals(int counter){
		
		int ub = counter%dlength;
		
		if(ub == 0) ub = dlength;
		
		for(int jj = 0 ; jj < ub ; jj++){
			xs += Xdata[jj]*Xdata[jj];
			x  += Xdata[jj];
			y  += Ydata[jj];
			xy += Xdata[jj]*Ydata[jj];
		}
		return;
	}
	
	private void getFitVals(){
		
		b  = (xs*y-x*xy)/(Ndata*xs-x*x);  	// Melissinos pg 450 (NB the equatiosn for a* and b* are mixed 
		m  = (Ndata*xy-x*y)/(Ndata*xs-x*x); // up but the equations for the errors in a* and b* are correct)		
		db = xs/(Ndata*xs-x*x);		
		db = sigma*Math.sqrt(db);
		dm = Ndata/(Ndata*xs-x*x);
		dm = sigma*Math.sqrt(dm);
		return;
	}
	
	private void calcChis(String fin, int skip, int xid, int yid){
		
		double psum = 0;
		int c1, c2, counter;
		boolean order;
		counter = 0;
		
		if(xid < yid){
			order = true;
			c1    = xid;
			c2    = yid;
		}
		else{
			order = false;
			c1    = yid;
			c2    = xid;
		}
	
		try {
			
			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));
			
			String rin;
			int pd;
			
			for (int jj = 0 ; jj < skip ; jj++){
				rin = bir.readLine();
			}
			
			while ( (rin = bir.readLine()) != null ){
				pd = rin.indexOf('\t');
				for(int jj = 1 ; jj < c1 ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(order){
					if(pd == -1){
						Xdata[counter] = Double.parseDouble(rin);
					}
					else{
						Xdata[counter] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				else{
					if(pd == -1){
						Ydata[counter++] = Double.parseDouble(rin);
					}
					else{
						Ydata[counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				for(int jj = c1 ; jj < c2 ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(order){
					if(pd == -1){
						Ydata[counter++] = Double.parseDouble(rin);
					}
					else{
						Ydata[counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				else{
					if(pd == -1){
						Xdata[counter++] = Double.parseDouble(rin);
					}
					else{
						Xdata[counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
				if(counter%dlength == 0){
					psum = doPsum(counter,psum);
					counter = 0;
				}
			}
			psum = doPsum(counter,psum);
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		
		chis = psum/(sigma*sigma);
		return;
	}
	
	private void calcChis(double[] XX, double[] YY){
		
		double psum = 0;
		
		if(Ndata > dlength){
			for(int jj = 0 ; jj < Ndata ; jj++){
				Xdata[jj%dlength] = XX[jj];
				Ydata[jj%dlength] = YY[jj];
				if(((jj+1)%dlength) == 0) psum = doPsum(dlength, psum); 
			}
			psum = doPsum(Ndata, psum);
		}
		else{
			Xdata = XX;
			Ydata = YY;
			psum = doPsum(Ndata, psum);
		}
		
		chis = psum/(sigma*sigma);
		return;
	}
	
	private double doPsum(int counter, double psum){
		
		int ub = counter%dlength;
		
		if(ub == 0) ub = dlength;
		
		for(int jj = 0 ; jj < ub ; jj++){
			psum += MathUtil.sqr(Ydata[jj] - (m*Xdata[jj]+b));
		}

		return psum;
	}
	
	public double getSlope(){
		
		return m;
	}
	
	public double getIntercept(){
		
		return b;
	}
	
	public double getSlopeErr(){
		
		return dm;
	}
	
	public double getInterceptErr(){
		
		return db;
	}
	
	public double getChis(){
		
		return chis;
	}
	
	
}
