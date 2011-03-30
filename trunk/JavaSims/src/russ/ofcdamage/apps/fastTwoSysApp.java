package russ.ofcdamage.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import chris.ofcdamage.damage2Dfast;
import chris.util.dummyParamUtil;

public class fastTwoSysApp extends Simulation{
	
	private int simt, eqt;
	private static int dl = 250000;
	private DecimalFormat pfmt = new DecimalFormat("00000");
	private damage2Dfast model[] = new damage2Dfast[2];
	private double st[];
	private double data[] = new double[dl];
	
	public static void main(String[] args) {
		new Control(new fastTwoSysApp(), "Damage Parameters");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Random Seed", (int) 0);
//		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond","All Sites"));
		params.add("Interaction Radius (R)", (int) 10);
		params.add("Lattice Size", (int) 256);
		params.add("Pass Stress to", new ChoiceValue("Live Sites","All Sites"));
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Open"));
		params.add("Equil Time", 500000);
		params.add("Sim Time", 500000);
//		params.add("Number of Lives",(int) 1);
//		params.add("NL width", (int) 0);
		params.add("Failure Stress (\u03C3_f)", 2.);
//		params.add("\u03C3_f width", 0.);
		params.add("Residual Stress (\u03C3_r)", 1.);
		params.add("\u03C3_r width", 0.1);
		params.add("Dissipation (\u03B1)", 0.05);
//		params.add("\u03B1 width", 0.);
		params.add("Mode");
		params.add("\u03D5");
		params.add("Status");
	}
	
	public void run() {

		double dphi = 1/((double)(params.iget("Lattice Size"))*(params.iget("Lattice Size")));
		double phin = 0.95;
		double phis = phin;
		int count   = (int)((1-phin)*(params.iget("Lattice Size"))*params.iget("Lattice Size"));
		simt = params.iget("Sim Time");
		eqt  = params.iget("Equil Time");
	
		
		while (phin < 1){
			
			params.set("Mode", "Freeze");
			params.set("\u03D5",count);

			model[0]  = new damage2Dfast(dummyParamUtil.ofcParams(params));
			model[1]  = new damage2Dfast(dummyParamUtil.ofcParams(params));
			// re-initialize the stress in second system with a new seed
			model[1].reRandomize((params.iget("Random Seed") + 1), model[0].getDorA());
			
			if(phin == phis) model[0].PrintParams(model[0].getOutdir()+File.separator+"Params_"+model[0].getBname()+".log",params,model[0]);				
			st = new double[2*model[0].getN()];

			// equilibrate model			
			params.set("Mode", "Equilibrating");
			//params.set("\u03D5",count);
			Job.animate();
			for (int jj = 0 ; jj < eqt ; jj++){
				for(int kk = 0 ; kk < 2 ; kk++)
					model[kk].evolveEQ(jj,false);
				if(jj%1000 == 0){
					params.set("Status", (jj-eqt));
					Job.animate();
				}
			}
			// simulate model for data
			params.set("Mode", "Simulating");
			params.set("Status", 0);
			params.set("\u03D5",count);
			Job.animate();
			model[0].setBname(params.sget("Data File")+"_1_"+pfmt.format(count),false);
			model[1].setBname(params.sget("Data File")+"_2_"+pfmt.format(count),false);
			writeHeaders();
			
			if(params.sget("Pass Stress to").equals("Live Sites")){
				for (int jj = 0 ; jj < simt ; jj++){
					for(int kk = 0 ; kk < 2 ; kk++)
						model[kk].evolveEQ(jj,true);
					calcMetric(jj);
					if(jj%1000 == 0){
						params.set("Status", (jj));
						Job.animate();
					}
				}
			}
			else{
				for (int jj = 0 ; jj < simt ; jj++){
					for(int kk = 0 ; kk < 2 ; kk++)
						model[kk].evolve(jj,false);
					calcMetric(jj);
					if(jj%1000 == 0){
						params.set("Status", (jj));
						Job.animate();
					}
				}
			}
			if((simt-1)%dl != 0) writeData(simt);
			phin += dphi;
			count--;
			//phin = (double)(Math.round(100*phin))/100;
			params.set("Random Seed", params.iget("Random Seed") + 2);
//			model = null; // clear model for avoid memory overflow
//			model = new damage2Dfast[2]; // re-initialize
		}
		params.set("Mode", "Done");
		Job.animate();
	
		return;
	}
	
	public void animate() {
		
		return;
	}

	public void clear() {

		return;
	}
	
	private void calcMetric(int t){
		
		int N = model[0].getN();
		double omega = 0;
				
		for(int jj = 0 ; jj < N ; jj++){
			for(int kk = 0 ; kk < 2 ; kk++)
				st[jj+kk*N] += model[kk].getStress(jj);
			omega += (st[jj] - st[jj+N])*(st[jj] - st[jj+N]);
		}
		omega /= (double)(t)*t;
		if(t%dl == 0 && t > 0){
			writeData(t);
		}
		data[t%dl] = 1/omega;
		
		return;
	}
	
	private void writeData(int mct){
		
		int ub;
		int offset = (int)((mct-1)/dl);
		offset = offset*dl;

		ub = mct%dl;
		if(ub==0) ub = dl;
		
		try{
			File file = new File(model[0].getOutdir()+File.separator+model[0].getBname()+"_II_"+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			for (int jj = 0 ; jj < ub ; jj++){
				pw.print(jj+offset);
				pw.print("\t");
				pw.println(data[jj]);
			}			
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}

	private void writeHeaders(){
		try{
			File file = new File(model[0].getOutdir()+File.separator+model[0].getBname()+".txt");
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			pw.print("time");
			pw.print("\t");
			pw.println("1/Omega");
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		return;
	}
	
}
