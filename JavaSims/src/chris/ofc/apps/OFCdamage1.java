package chris.ofc.apps;

import java.awt.Color;
import java.util.Random;

import kip.util.LatticeNeighbors;
import scikit.dataset.Histogram;
import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DoubleValue;

public class OFCdamage1 extends Simulation {
	
	public int time,L,R,N,imax, showernumber;
	public int alive[], dead[];
	public double alpha,Sc,dS,release,stressMax;
	public double stress[];
	public Boolean lanch;
	Grid grid1 = new Grid ("Stress Lattice");
	Grid grid2 = new Grid ("Failed Sites");
	Plot plot = new Plot("Histogram of Number of Showers");
	Histogram histNS;
	LatticeNeighbors neighbors;

	public static void main(String[] args) {
		new Control(new OFCdamage1(), "OFC Model");
	}
	
	public void load(Control c) {
		params.add("Random Seed",0);
		params.addm("Auto Scale", new ChoiceValue("Yes", "No"));
		params.add("Lattice Size",1<<9);
		params.add("Boundary Condtions", new ChoiceValue("PERIODIC","BORDERED"));
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.addm("Interaction Radius (R)",(int)(50));		
		params.addm("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("Number of Resets");
		params.add("Number of Showers");
		
		
		c.frameTogether("OFC Model with Damage (Between Avalanches)", grid1, grid2);
		c.frame(plot);
		
	}
	
	
	
	
	public void animate() {

		int i;
		
		plot.setAutoScale(false);
		if (params.sget("Auto Scale").equals("Yes")) {
			plot.setAutoScale(true);
		}
		

		ColorPalette palette = new ColorPalette();
		palette.setColor(0,Color.BLACK);
		palette.setColor(1,Color.WHITE);
		
		ColorGradient smooth = new ColorGradient();
		for (i=0 ; i<N ; i++){
			smooth.getColor(stress[i],-2,Sc);
		}
		
		grid1.setColors(smooth);
		grid2.setColors(palette);
		grid1.registerData(L,L,stress);
		grid2.registerData(L, L, alive);	
		
		params.set("Number of Resets",time);
		params.set("Number of Showers",showernumber);

		plot.registerBars("histNshowers", histNS, Color.RED);
		
	}

	public void clear() {
		grid1.clear();
		grid2.clear();
		plot.clear();
	}

	public void run() {
		
		histNS = new Histogram(1.);
		int Nalive;
		
		time=0;
		showernumber=0;
		
		
		L=    params.iget("Lattice Size");
		R=    params.iget("Interaction Radius (R)");
		alpha=params.fget("Dissipation (\u03B1)");
		Sc=   params.fget("Critical Stress (\u03C3_c)");
		N=L*L;
		
		stress = new double[N];
		alive  = new int[N];
		dead   = new int[N];

		if (params.sget("Boundary Condtions").equals("BORDERED")){
			neighbors = new LatticeNeighbors(L, L, 0, R, LatticeNeighbors.Type.BORDERED);	
		}
		else{
			neighbors = new LatticeNeighbors(L, L, 0, R, LatticeNeighbors.Type.PERIODIC);
		}
		Random rand = new Random(params.iget("Random Seed"));
		

		// Initialize Sites
	
		imax=0;	
		for (int i = 0; i<N; i++){
			stress[i]=Sc*rand.nextDouble();		// How random??
			if(stress[i]>stress[imax]) imax=i;
			alive[i]=1;
		}
		alive[imax]=0;
		
//		for (int i = 0; i<N; i++){
//			stress[i]=0.1*Sc*rand.nextDouble();	// How random??
//			alive[i]=1;
//		}
//		imax=(int)(N/2+L/2);
//		stress[imax]=Sc;
//		alive[imax]=0;
		
		
		// Bring most stressed site to failure
		
		stressMax = stress[imax];
		for (int i = 0; i<N; i++){
			stress[i]+=Sc-stressMax;
		}
		alive[imax]=0;		
			
		
		while(true){
			
			//Job.animate();
			
			// Redistribute stress from failed site
			
			int[] nbs = neighbors.get(imax);
	
			Nalive = 0;
			for (int i = 0; i<nbs.length; i++){
				Nalive+=alive[nbs[i]];
			}
			if(Nalive>0){
				release=(1-alpha)*stress[imax]/Nalive;
				for (int i = 0; i<nbs.length; i++){
					stress[nbs[i]]+=release*alive[nbs[i]];
				}
			}
			stress[imax]=-2;
			
			// Look for avalanche
			
			int search=0;
			for (int i = 0; i<N ; i++){
				if(stress[i]>Sc){
					dead[search]=i;
					alive[i]=0;
					search++;
				}
			}
			
			while (search>0){
				showernumber++;
				// Redistribute avalanche stress
				
				for (int i = 0; i<search; i++){
					
					nbs = neighbors.get(dead[i]);
					Nalive=0;
					for (int j = 0; j<nbs.length; j++){
						Nalive+=alive[nbs[j]];
					}			
					if(Nalive>0){									// what do we do with this stress if "if" fails??
						release=(1-alpha)*stress[dead[i]]/Nalive;
						for (int j = 0; j<nbs.length; j++){
							stress[nbs[j]]+=release*alive[nbs[j]];
						}
					}
					stress[dead[i]]=-2;
				}
				
				// Look for subsequent avalanche
				
				search=0;
				for (int i = 0; i<N ; i++){
					if(stress[i]>Sc){
						dead[search]=i;
						alive[i]=0;
						search++;
					}
				}
				histNS.accum((double)(time));
				Job.animate();
			}
		
			// Reset stresses to create failure
			
			for (int i = 0; i<N; i++){
				if(stress[i]>=stress[imax]) imax=i;
			}
			
			stressMax = stress[imax];
			
			if(alive[imax]==0){
				System.out.println("All Sites Failed!");
				break;
			}

			for (int i = 0; i<N; i++){
				stress[i]+=Sc-stressMax;
			}

			alive[imax]=0;
			
			showernumber=0;
			time++;
			Job.animate();
			
		}
		
	}
			
}

	
	
	
