package chris.util;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;

import scikit.jobs.params.Parameters;

public class ReadInUtil {

	private String fin;

	public ReadInUtil(String filein){

		fin = filein;

		return;
	}
	public double[] getData(int skip){

		int counter = 0;
		double[] values = new double[1000000];
		String rin;
		
		try {

			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));

			for (int jj = 0 ; jj < skip ; jj++){
				rin = bir.readLine();
			}
			while ( (rin = bir.readLine()) != null ){
				values[counter++] = Double.parseDouble(rin);
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}


		double[] ret = new double[counter];

		for (int jj = 0 ; jj < counter ; jj++){
				ret[jj] = values[jj];
		}
		return ret;
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

	public int countTo(int skip, String str){
		
		int counter = 0;

		try {

			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));

			String rin;
			for (int jj = 0 ; jj < skip ; jj++)
				rin = bir.readLine();
			
			while ( (rin = bir.readLine()) != null){
				if(rin.equals(str))
					break;
				counter++;
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}


		return counter;
	}
	
	public double[][] getDataBeforeString(int[] cns, int skip, String str){
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

				while ( (rin = bir.readLine()) != null){
					pd = rin.indexOf('\t');
					for(int jj = 1 ; jj < cn ; jj++){
						rin = rin.substring(pd + 1);
						pd = rin.indexOf('\t');
					}
					if(pd == -1){
						if(rin.equals(str))
							break;
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
	// ReadIn(directory, keyword, phase data, parameter data)

	public void getDataandParams(int[] cns, int skip, double[][] data, String str, Parameters p){
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

				while ( (rin = bir.readLine()) != null){
					pd = rin.indexOf('\t');
					for(int jj = 1 ; jj < cn ; jj++){
						rin = rin.substring(pd + 1);
						pd = rin.indexOf('\t');
					}
					if(pd == -1){
						if(rin.equals(str)){
							getParams(bir, p);
							break;
						}
							
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

		for (int jj = 0 ; jj < counter ; jj++){
			for(int kk = 0 ; kk < Ncols ; kk++){
				data[kk][jj] = values[kk][jj];
			}
		}
		return;
	}
	
	// ONLY FOR LENNARD-JONES APPS
	private void getParams(BufferedReader bir, Parameters p){
		
		String rin;
		int pd;
		
		try{
			if(p.containsKey("L")){
				rin = bir.readLine(); // seed, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("L",Double.parseDouble(rin.substring(pd + 2)));
				rin = bir.readLine(); 
				pd  = rin.indexOf('=');
				p.set("Boundary Conditions",rin.substring(pd + 2));
				rin = bir.readLine(); // IC, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("ODE Solver",rin.substring(pd + 2));
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("N",Integer.parseInt(rin.substring(pd + 2)));
				rin = bir.readLine(); // M, throw it away
				rin = bir.readLine(); // R, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("dt",Double.parseDouble(rin.substring(pd + 2)));
				// throw away everything else
			}
			else{
				rin = bir.readLine(); // seed, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("Lx",Double.parseDouble(rin.substring(pd + 2)));
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("Ly",Double.parseDouble(rin.substring(pd + 2)));
				rin = bir.readLine(); 
				pd  = rin.indexOf('=');
				p.set("Boundary Conditions",rin.substring(pd + 2));
				rin = bir.readLine(); // IC, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("ODE Solver",rin.substring(pd + 2));
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("N",Integer.parseInt(rin.substring(pd + 2)));
				rin = bir.readLine(); // M, throw it away
				rin = bir.readLine(); // R, throw it away
				rin = bir.readLine();
				pd  = rin.indexOf('=');
				p.set("dt",Double.parseDouble(rin.substring(pd + 2)));
				// throw away everything else
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}

		return;
	}
	
	public Parameters getOFCparams(){

		Parameters p   = new Parameters();
		Parameters isD = new Parameters();
		Parameters isI = new Parameters();
		
		isD.add("Failure Stress (?_f)");
		isD.add("?_f width");
		isD.add("Residual Stress (?_r)");
		isD.add("?_r width");
		isD.add("Dissipation (?)");
		isD.add("? width");
		isI.add("Random Seed");
		isI.add("Interaction Radius (R)");
		isI.add("Lattice Size");
		isI.add("Equil Time");
		isI.add("Sim Time");
		
		try {

			FileInputStream fis = new FileInputStream(fin);
			BufferedInputStream bis = new BufferedInputStream(fis);
			BufferedReader bir = new BufferedReader(new InputStreamReader(bis));

			String rin;
			int pd;

			while ( (rin = bir.readLine()) != null ){
				pd = rin.indexOf('=');
				if(pd < 0)
					break;
				if(isD.containsKey(rin.substring(0,pd-1))){
					if(rin.substring(0,pd-1).indexOf("Failure") >= 0){
						p.add("Failure Stress (\u03C3_f)", Double.parseDouble(rin.substring(pd+1)));
					}
					else if(rin.substring(0,pd-1).indexOf("?_f") >= 0){
						p.add("\u03C3_f width", Double.parseDouble(rin.substring(pd+1)));
					}			
					else if(rin.substring(0,pd-1).indexOf("Residual") >= 0){
						p.add("Residual Stress (\u03C3_r)", Double.parseDouble(rin.substring(pd+1)));
					}
					else if(rin.substring(0,pd-1).indexOf("?_r") >= 0){
						p.add("\u03C3_r width", Double.parseDouble(rin.substring(pd+1)));
					}					
					else if(rin.substring(0,pd-1).indexOf("Dissipation") >= 0){
						p.add("Dissipation (\u03B1)",Double.parseDouble(rin.substring(pd+1)));
					}
					else if(rin.substring(0,pd-1).indexOf("?") >= 0){
						p.add("\u03B1 width", Double.parseDouble(rin.substring(pd+1)));
					}
					else{
						
					}
				}
				else if(isI.containsKey(rin.substring(0,pd-1))){
					p.add(rin.substring(0,pd-1),Integer.parseInt(rin.substring(pd+2)));
				}
				else{
					p.add(rin.substring(0,pd-1),rin.substring(pd+2));
				}			
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		return p;
	}
	
}
