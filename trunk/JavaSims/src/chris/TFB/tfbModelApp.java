package chris.TFB;

import java.awt.Color;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

import java.io.IOException;
import javax.imageio.ImageIO;

public class tfbModelApp extends Simulation{
	
	private Grid grid = new Grid("Fibre State");
	private ColorPalette palette1;
	
	private double stress, Stot, Sf[], Pheal, kappa, beta;
	private int L, N, Nalive, state[], t, Nh, Nf;
	private Random rand;
	private boolean draw;
	private String dir;
	
	private DecimalFormat fmt = new DecimalFormat("000000");
	
	
	public static void main(String[] args) {
		new Control(new tfbModelApp(), "TFB Model");
	}
	
	
	public void load(Control c) {
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("Random Seed",0);
		params.add("Lattice Size",1<<8);
		params.add("Stress / site",1.0);		
		params.add("Failure Stress",5.0);
		params.add("\u03C3_f width",(double) 0);
		params.add("Kappa", (double) 1);
		params.addm("Temperature", new DoubleValue(1, 0, 5000).withSlider());
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("MC Time Step");
		params.add("Stress");
		params.add("Num Alive");
		
		c.frame(grid);
	}
	
	public void run() {
		
		Init();
		
		while(evolve()){
			t++;
			takeData();
			Job.animate();
		}
		
	}
	
	public void animate() {
		
		params.set("MC Time Step", (int) t);
		params.set("Stress",stress);
		params.set("Num Alive", Nalive);
		
		if(draw){
			grid.registerData(L,L,state);
			if(params.sget("Record").equals("On")) takePicture();
		}
	}

	
	public void clear() {
		grid.clear();
	}
	
	public void Init(){
		
		double Sf0, SfW;
		
		stress = params.fget("Stress / site");
		L      = params.iget("Lattice Size");
		N      = L*L;
		Stot   = N*stress;
		beta   = 1/(params.fget("Temperature"));
		kappa  = params.fget("Kappa");
		rand   = new Random(params.iget("Random Seed"));
		Sf     = new double[N];
		state  = new int[N];
		Sf0    = params.fget("Failure Stress");
		SfW    = params.fget("\u03C3_f width");
		Nalive = N;
		draw   = (params.sget("Animation").equals("On"));
		dir    = params.sget("Data Directory");
		
		if(SfW == 0){
			for(int jj = 0 ; jj < N ; jj++){
				Sf[jj] = Sf0;
				state[jj] = 1;
			}
		}
		else{
			for(int jj = 0 ; jj < N ; jj++){
				Sf[jj] = Sf0 + SfW*(0.5 - rand.nextDouble());
				state[jj] = 1;
			}
		}
		
		palette1  = new ColorPalette();
		palette1.setColor(0,Color.BLACK);
		palette1.setColor(1,Color.WHITE);
		grid.setColors(palette1);
	}
	
	public boolean evolve(){
		
		beta = 1/(params.fget("Temperature"));
		
		Pheal = Math.exp(-beta*stress*stress/kappa);
		
		int site = rand.nextInt(N);
		state[site] = nextState(site);
	
		return (Nalive > 0);
	}
	
	public int nextState(int site){
		
		if(state[site] == 1){
			//try to kill it
			if(rand.nextDouble() < Math.exp(-beta*(Sf[site]-stress)*(Sf[site] - stress)/kappa)){
				// kill site
				Nalive--;
				stress = Stot/Nalive;
				Nf++;
				return 0;
			}
			else{
				// do nothing?
				return 1;
			}
			
		}
		else{
			//try to heal it
			if(rand.nextDouble() < Pheal){
				// heal site
				Nalive++;
				stress = Stot/Nalive;
				Nh++;
				return 1;
			}
			else{
				// do nothing?
				return 0;
			}
			
			
		}
		
	}
	
	public void takeData(){
		
		return;
	}
	
	public void takePicture(){
			


		String SaveAs = dir + File.separator + grid.getTitle()+fmt.format(t)+".png";

		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}

	

		return;
	}

}
