package chris.TFB;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

import scikit.graphics.dim2.Grid;
import scikit.jobs.params.Parameters;

public class TFBmodel {

	@SuppressWarnings("unused")
	private double stress, L, N, Stot, beta, kappa, Sf[], Sf0, SfW, t, metric, sCum[];
	private int state[], Nalive;
	private String dir;
	@SuppressWarnings("unused")
	private boolean draw;
	private Random rand;
	
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	public TFBmodel(Parameters params){
		
		tfbConstructor(params);
		
		return;
	}
	
	public void tfbConstructor(Parameters params){
		
		stress = params.fget("Stress / site");
		L      = params.iget("Lattice Size");
		N      = L*L;
		Stot   = N*stress;
		beta   = 1/(params.fget("Temperature"));
		kappa  = params.fget("Kappa");
		rand   = new Random(params.iget("Random Seed"));
		Sf     = new double[(int) N];
		state  = new int[(int) N];
		sCum   = new double[(int) N];
		Sf0    = params.fget("Failure Stress");
		SfW    = params.fget("\u03C3_f width");
		Nalive = (int) N;
		draw   = (params.sget("Animation").equals("On"));
		dir    = params.sget("Data Directory");
		t      = 0;
		
		initArrays();
	}
	
	private void initArrays(){
		
		for (int jj = 0 ; jj < N ; jj++){
			Sf[jj]    = Sf0 + SfW*(0.5 - rand.nextDouble());
			state[jj] = 2;
			sCum[jj]  = 0;
		}
		
		return;
	}
	
	public boolean nextLattice(){
		
		t++;
		int deltaN = 0;
		
		for(int jj = 0 ; jj < N ; jj++){
			if(nextState(jj)){
				state[jj] = 1-state[jj];	// state = 0 --> 1 // state = 1 --> 0 
				deltaN += 2*state[jj] - 1;  // state = 0 --> -1 // state = 1 --> +1
			}
		}
		
		Nalive += deltaN;
		
		return (Nalive > 0);
	}
	
	private boolean nextState(int st){
		
		if(state[st] == 1){
			// try to fail it	
			
			/*
			 *	If site passes failure requirement, return TRUE   
			 */
		}
		
		else{
			// try to heal it 
			/*
			 *	If site passes healing requirement, return TRUE   
			 */
		}
		
		
		return false;
	}
	public void takeData(){
		
		return;
	}
	
	public void takePicture(Grid grid){

		String SaveAs = dir + File.separator + grid.getTitle()+fmt.format(t)+".png";
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}

		return;
	}
	
	@SuppressWarnings("unused")
	private void ergMetric(){
		
		
		return;
	}
	
}
