package chris.queue.Apps;

import java.io.File;

import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.router2D;
import chris.util.PrintUtil;

public class router2Dapp extends Simulation{

	Grid grid = new Grid("Buffers");
	int L, N;
	router2D model;
	private final boolean draw = false;
	String pth, fout;
	
	public static void main(String[] args) {
		new Control(new router2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("L",32);
		params.add("l",5);
		params.add("\u03BB",new DoubleValue(0.05,0,1));
		params.add("seed",0);
		params.add("messages");
		params.set("messages",0);
		params.add("t");
		params.set("t",0);

		c.frame(grid);
	}

	public void run() {
		
		Histogram foo = new Histogram(1000.);

		L     = params.iget("L");
		N     = L*L;
		model = new router2D(params); 
		model.writeDataHeader();

		setupcolors(params.fget("\u03BB"),params.iget("l"));
		
		pth   = params.sget("Data Directory");
		fout  = params.sget("Data File");
		
		
		int tss  = 50000;
		int tmax = 500000;
		int t    = 0;
		
		while(t++ < tss){
			model.step(1,false,foo);
			if(t%1e3 ==0){
				params.set("t",t-tss);
				Job.animate();
			}
		}
		t = 0;
		int[] ts     = new int[tmax];
		Histogram hN = new Histogram(1.);
		Histogram ht = new Histogram(1.);
		while(t++ < tmax){
			model.step(1,false,ht);
			hN.accum(model.getNmsg());
			ts[t-1] = model.getNmsg();
			if(t%1e4 ==0){
				params.set("t",t);
				Job.animate();
			}
		}
		printData(hN, ht, ts, model.getNmsgHist());
		
		params.set("t","Done");
		Job.signalStop();
		Job.animate();
		
	}
	
	private void printData(Histogram h1, Histogram h2, int[] v, Histogram h3){
		String f1 = pth + File.separator + fout+"_nhist_.txt";
		String f2 = pth + File.separator + fout+"_thist_.txt";
		String f3 = pth + File.separator + fout+"_tsers_.txt";
		String f4 = pth + File.separator + fout+"_nindv_.txt";
		
		PrintUtil.printHistToFile(f1, h1);
		PrintUtil.printHistToFile(f2, h2);	
		PrintUtil.printVectorToFile(f3,v);
		PrintUtil.printHistToFile(f4, h3);	

		return;
	}
	
	private void setupcolors(double lambda, int l){

		int nmbar = (int)((lambda*lambda*(l-1)/(1-l*lambda) + l*lambda));
		if(nmbar == 0){
			nmbar = 5;
		}
		else{
			nmbar = 2*(nmbar+1)*N;
		}
		setupcolors(nmbar);
		return;
	}

	private void setupcolors(int mx){
		
		ColorGradient cg = new ColorGradient();
		ColorPalette cp  = new ColorPalette();

		for (int jj = 0 ; jj <= mx ; jj++){
			cp.setColor(mx-jj,cg.getColor(jj, 0, mx));
		}
		grid.setColors(cp);
		return;
	}
	
	public void animate() {
		
		//params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		if(draw) 
			grid.registerData(L,L,model.getDensity());
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}	
}
