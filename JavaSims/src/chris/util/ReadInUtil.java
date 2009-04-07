package chris.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

public class ReadInUtil {

	private String fin;
	
	public ReadInUtil(String filein){
		
		fin = filein;
		
		return;
	}
	
	public double[][] getData(int[] cns, int skip){
		
		int counter = 0;
		int Ncols   = cns.length;
		double[][] values = new double[Ncols][1000000];
		
		for(int ii = 0 ; ii < Ncols ; ii++){
		
			int cn = cns[ii];
			counter = 0;
			
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
					for(int jj = 1 ; jj < cn ; jj++){
						rin = rin.substring(pd + 1);
						pd = rin.indexOf('\t');
					}
					if(pd == -1){
						values[ii][counter++] = Double.parseDouble(rin);
					}
					else{
						values[ii][counter++] = Double.parseDouble(rin.substring(0,pd));	
					}
				}
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
		
		double[][] ret = new double[Ncols][counter];
		
		for (int jj = 0 ; jj < counter ; jj++){
			for(int kk = 0 ; kk < Ncols ; kk++){
				ret[kk][jj] = values[kk][jj];
			}
		}
		
		return ret;

	}
		
	public void getData(int[] cns, int skip, double[] retX, double[] retY){
		
		int overflow = retX.length;
		
		int counter = 0;
		
		int xn = cns[0];
		int yn = cns[1];
		
		try {
			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));

			String rin;
			int pd;

			for (int jj = 0 ; jj < skip ; jj++){
				rin = bir.readLine();
				overflow--;
			}

			while ( (rin = bir.readLine()) != null ){

				if(overflow-- <= 0){
					break;
				}
				
				pd = rin.indexOf('\t');
				for(int jj = 1 ; jj < xn ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				retX[counter] = Double.parseDouble(rin.substring(0,pd));
				
				for(int jj = (xn-1) ; jj < yn ; jj++){
					rin = rin.substring(pd + 1);
					pd = rin.indexOf('\t');
				}
				if(pd == -1){
					retY[counter++] = Double.parseDouble(rin);
				}
				else{
					retY[counter++] = Double.parseDouble(rin.substring(0,pd));	
				}
			}
		} 
		catch (FileNotFoundException e) {
			e.printStackTrace();
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
	
		return;
	}
	
}
