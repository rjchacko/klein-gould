package chris.TFB.Apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import chris.TFB.TFB;
import chris.util.PrintUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;

public class newTFBapp extends Simulation{

	int mct;
	double data[][];
	TFB model;
	
	public static void main(String[] args) {
		new Control(new newTFBapp(), "TFB Model");
	}
	
	public void animate() {
	
		return;
	}

	public void clear() {
		return;
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop"));
		params.add("File Name","tfb");
		params.add("seed",(int) 0);
		params.add("Sim Time", (int) 1e4); //1e6
		params.add("N", (int)(256*256)); // 256*256
		params.add("stress / fibre", 1e-5); 
		params.add("K", 1);
		params.add("D", 1);
		params.add("T", 0.5);
		params.add("phi");
		params.add("mct");
	}

	public void run() {
		
		params.set("mct","Initializing");
		Job.animate();
		
		mct   = params.iget("Sim Time");
		data  = new double[4][mct];
		model = new TFB(params);
		
		for (int jj = 1 ; jj < mct ; jj++){
			
			model.nextBundle();
			model.calcOmega(jj);
			data[0][jj] = model.getPhi();
			data[1][jj] = model.getE();
			data[2][jj] = model.getOmega1inv();
			data[3][jj] = model.getOmega2inv();
			if(jj % 1000 == 0){
				params.set("phi", data[0][jj]);
				params.set("mct",jj);
				Job.animate();
			}
		}
		
		PrintUtil.printlnToFile(params.sget("Data Directory")+File.separator+"Params_"+params.sget("File Name")+".txt",params.toString());
		try{
			File file = new File(params.sget("Data Directory")+File.separator+params.sget("File Name")+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("mct");
			pw.print("\t");
			pw.print("phi");
			pw.print("\t");
			pw.print("E");
			pw.print("\t");
			pw.print("1/Omega_e");
			pw.print("\t");
			pw.println("1/Omega_phi");	
			pw.print(0);
			pw.print("\t");
			pw.print(1);
			pw.print("\t");
			pw.print(params.fget("stress / fibre")*params.fget("stress / fibre")/(2*params.fget("K")*params.iget("N"))-params.fget("D"));
			pw.print("\t");
			pw.print("NaN");
			pw.print("\t");
			pw.println("NaN");	
			for (int jj = 1 ; jj < mct ; jj++){
				pw.print(jj);
				pw.print("\t");
				pw.print(data[0][jj]);
				pw.print("\t");
				pw.print(data[1][jj]);
				pw.print("\t");
				pw.print(data[2][jj]);
				pw.print("\t");
				pw.println(data[3][jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		params.set("mct","Done");
		Job.animate();
		Job.signalStop();
		Job.animate();
		
	}

}
