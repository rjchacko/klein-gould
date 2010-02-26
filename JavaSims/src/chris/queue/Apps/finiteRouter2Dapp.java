package chris.queue.Apps;

import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.queue.finiteRouter2D;
import chris.util.movieUtil;

public class finiteRouter2Dapp extends Simulation{

	Grid grid  = new Grid("Buffers");
	Grid movie = new Grid("buffers");
	private int L, N, Mx;
	private String pth;
	private finiteRouter2D model;
	private boolean draw, rec;

	public static void main(String[] args) {
		new Control(new finiteRouter2Dapp(), "2D Router Network");
	}
	
	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Data File", "default");
		params.add("Draw", new ChoiceValue("Animate","Record","Off"));
		params.addm("M",10);
		params.add("L",32);
		params.add("l",5);
		params.addm("\u03BB",new DoubleValue(0.05,0,1).withSlider());
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
		Mx    = N*params.iget("M");
		//msg   = new int[10000];
		model = new finiteRouter2D(params); 
		pth   = params.sget("Data Directory");
        draw  = !(params.sget("Draw").equals("Off")); 
		rec   = false;
		
        if(draw)
        	setupcolors(params.iget("M"));
		
		while(model.step(1,false,foo) < Mx){
			//msg[model.getT()] = model.getNmsg();
			Job.animate();
		}
		Job.animate();

		//PrintUtil.printVectorToFile("/Users/cserino/Desktop/NmsgT.txt",msg,(int)(model.getT()));
		
		params.set("t","Done");
		Job.signalStop();
		Job.animate();
		
	}
	
	private void setupcolors(int mx){
		
		ColorGradient cg = new ColorGradient();
		ColorPalette cp  = new ColorPalette();

		for (int jj = 0 ; jj <= mx ; jj++){
			cp.setColor(mx-jj,cg.getColor(jj, 0, mx));
		}
		//public static void saveLegend(int min, int max, ColorPalette cp, String pth){

		grid.setColors(cp);
		movie.setColors(cp);
		
		rec = params.sget("Draw").equals("Record"); 
		return;
	}
	
	public void animate() {

		model.setM(params.iget("M"));
		model.resetLambda(params.fget("\u03BB"));
		params.set("t", model.getT());
		params.set("messages",model.getNmsg());
		if(draw){ 
			grid.registerData(L,L,model.getDensity());
			if(rec) 
				movieUtil.saveImage(model.getDensity(), L, L, 20, movie, pth, model.getT());
		}
		return;
	}

	public void clear() {
		
		grid.clear();
		return;
	}	
}
