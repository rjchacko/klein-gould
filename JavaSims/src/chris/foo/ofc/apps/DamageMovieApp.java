package chris.foo.ofc.apps;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

import scikit.graphics.ColorGradient;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

public class DamageMovieApp extends Simulation{

	private double stress[], sbar[], sr[], sr0, dsr, sf, alpha, dt, t;
	private int R, N, eqtime, simtime, search[], qN, lives[], maxL, L, ln, mct, Nalive;
	private boolean record, alive[];
	private String fout, imode, bname, BCs, PicDir;
	private Random rand;
	private LatticeNeighbors nbs;
	private ColorGradient cg = new ColorGradient();
	private DecimalFormat fmtS = new DecimalFormat("000");
	private DecimalFormat fmtL = new DecimalFormat("0000000");

	
	Grid gridA = new Grid ("Alive");
	Grid gridR = new Grid ("Density");

	
	public static void main(String[] args) {
		new Control(new DamageMovieApp(), "OFC Parameters");
	}
	
	public void animate() {
		
		if(record){
			
			double[] rho = getRho();
			for (int jj=0 ; jj < N ; jj++){
				rho[jj] = rho[jj]/qN;
				cg.getColor(rho[jj],0.,1.);
			}
			gridR.setColors(cg);  // Do I only need to do this once??? C-H-E-C-K ! ! !
			// do this by setting the scale properly maybe??
			gridA.registerData(L, L, bool2bin(alive));
			gridR.registerData(L, L, rho);
			TakePicture(gridA);
			TakePicture(gridR);
		}
		
		return;
	}

	public void clear() {
		
		if(record){
			gridA.clear();
			gridR.clear();
		}
		return;
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("File Name","FailSafe");
		params.add("Initialization Mode", new ChoiceValue("Random","Site Tracking"));
		params.add("Random Seed",0);
		params.add("Boundary Conditions", new ChoiceValue("Peroidic", "Open"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Equilibrate Time",(int)200000);
		params.add("Simulation Time",(int)200000);
		params.add("Number of Lives", (int) 1);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("Residual Stress (\u03C3_r)",1.0);
		params.add("\u03C3_r width",(double) 0.5);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("Record", new ChoiceValue("On","Off"));
		params.add("Status");
		
		c.frameTogether("Stress",gridA,gridR);
	}

	public void run() {
		
		while(true){
		
			Initialize();

			record = false;
			// equilibrate
			for (int jj = 0 ; jj < eqtime ; jj++){
				Evolve();
				if(jj%200 == 0) params.set("Status",jj - eqtime);
				Job.animate();
			}

			// simulate as EQ model
			for (int jj = 0 ; jj < simtime ; jj++){
				Evolve();
				t += dt;
				mct++;
				//saveData(jj);
				if(jj%200 == 0) params.set("Status", jj);
				Job.animate();
			}

			// simulate with damage		
			if(maxL > 0){	
				record = (params.sget("Record")=="On");
				while (Nalive > 0){
					EvolveD();
					t += dt;
					mct++;
					params.set("Status", mct);
					Job.animate();
				}
			}

			// end
			params.set("Status","Done");
			Job.signalStop();
			Job.animate();

		}
				
	}
	
	private void Initialize() {
		
		int tmp, tx;
		
		fout    = params.sget("Data Directory");
		bname   = params.sget("File Name");
		imode   = params.sget("Initialization Mode");
		tmp     = params.iget("Random Seed");
		BCs     = params.sget("Boundary Conditions");
		R       = params.iget("Interaction Radius (R)");
		L       = params.iget("Lattice Size");
		eqtime  = params.iget("Equilibrate Time");
		simtime = params.iget("Simulation Time");
		maxL    = params.iget("Number of Lives");
		sf      = params.fget("Failure Stress (\u03C3_f)");
		sr0	    = params.fget("Residual Stress (\u03C3_r)");
		dsr     = params.fget("\u03C3_r width");
		alpha   = params.fget("Dissipation (\u03B1)");
		record  = (params.sget("Record")=="On");
		
		PrintUtil.printlnToFile(fout + File.separator + "Params_" + bname + ".txt", params.toString());
		
		mct  = 0;
		t    = 0;
		dt   = 0;
		N    = L*L;
		fout += File.separator + bname + ".txt";
		rand = new Random(tmp);
		
		sbar   = new double[N];
		stress = new double[N];
		sr     = new double[N];
		search = new int[N];
		if(maxL>0){
			Nalive = N;
			lives = new int[N];
			alive = new boolean[N];
		}		
		
		
		for (int jj = 0 ; jj < N ; jj++){
			sbar[jj] = 0;
			if(imode == "Random"){
				stress[jj] = sr0 + (sf - sr0)*rand.nextDouble();
				sr[jj]     = sr0 + dsr*2*(rand.nextDouble() - 0.5);
			}
			else{
				stress[jj] = sr0;
				sr[jj]     = sr0;
				stress[(int)(N/2) + (int)(L/2)] = 1.1*sr0;
			}		
			search[jj] = 0;
			if(maxL>0){
				lives[jj] = maxL;
				alive[jj] = true;
			}
		}
		
		
		if(BCs == "Peroidic"){
			nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
		}
		else if(BCs == "Open"){
			nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
		}
		tx = (int)(L/2);
		qN  = (nbs.get(L*tx + tx)).length;

		if(record){
			PicDir = params.sget("Data Directory") + "/GridPics"; 
			DirUtil.MkDir(PicDir);
			ColorPalette palette1  = new ColorPalette();
			palette1.setColor(1,Color.WHITE);
			palette1.setColor(0,Color.BLACK);
			gridA.setColors(palette1);
		}
	
		return;
	}
	
	private void Evolve(){
		
		int jmax = 0;
		int counter = 0;
		double smax;
		int nlist[], nlength;
		double release;
		
		// force the first failure
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj] > stress[jmax]) jmax = jj;
		}
		smax = stress[jmax];
		dt   = sf-smax;
		for (int jj = 0 ; jj < N ; jj++){
			stress[jj] += (sf - smax);
		}
		
		search[counter++] = jmax;
		
		while(counter > 0){
			
			// distribute stress
			for (int jj = 0 ; jj < counter ; jj++){
				nlist   = nbs.get(search[jj]);
				nlength = nlist.length;
				release = (1-alpha)*(stress[search[jj]] - sr[search[jj]]) / qN;
				for (int kk = 0 ; kk < nlength ; kk++){
					stress[nlist[kk]] += release; 
				}
			}
			
			// reset the failed sites
			for (int jj = 0 ; jj < counter ; jj++){
				sr[search[jj]]     = nextSr();
				stress[search[jj]] = sr[search[jj]];
			}

			// check and see if the lattice is stable
			counter = 0;
			for (int jj = 0 ; jj < N ; jj++){
				if(stress[jj] > sf) search[counter++] = jj;
			}
			
		}
		
		return;
	}
	
	private void EvolveD(){
		
		int jmax = 0;
		int counter = 0;
		double smax;
		int nlist[], nlength, countLive;
		double release;
		ln = 1;
		
		// force the first failure
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj] > stress[jmax]) jmax = jj;
		}
		smax = stress[jmax];
		dt   = sf-smax;
		for (int jj = 0 ; jj < N ; jj++){
			stress[jj] += (sf - smax);
		}
		
		lives[jmax]--;
		if(lives[jmax]==0){
			alive[jmax] = false;
			Nalive--;
		}

		
		search[counter++] = jmax;
		
		while(counter > 0){
			
			// distribute stress
			for (int jj = 0 ; jj < counter ; jj++){
				nlist   = nbs.get(search[jj]);
				nlength = nlist.length;
				countLive = 0;
				for (int kk = 0 ; kk < nlength ; kk++){
					countLive += bool2bin(alive[nlist[kk]]);
				}
				if(countLive==0) continue;
				release = (1-alpha)*(stress[search[jj]] - sr[search[jj]]) / countLive;
				for (int kk = 0 ; kk < nlength ; kk++){
					stress[nlist[kk]] += bool2bin(alive[nlist[kk]])*release; 
				}
			}
			
			// reset the failed sites
			for (int jj = 0 ; jj < counter ; jj++){
				sr[search[jj]]     = nextSr();
				stress[search[jj]] = sr[search[jj]];
			}

			// check and see if the lattice is stable
			counter = 0;
			for (int jj = 0 ; jj < N ; jj++){
				if(stress[jj] > sf){
					search[counter++] = jj;
					lives[jj]--;
					if(lives[jj]==0){
						alive[jj] = false;
						Nalive--;
					}
				}
				
			}
			Job.animate();
			ln++;
		}
		
		return;
	}
	
	private int bool2bin(boolean b){
	
		return b ? 1 : 0;
	}
	
	private int[] bool2bin(boolean b[]){
		
		int lth = b.length;
		int[] ret = new int[lth];
		
		for (int jj = 0 ; jj < lth ; jj++){
			ret[jj] = bool2bin(b[jj]);
		}
		
		return ret;
	}
	
	private double nextSr(){
		
		return (dsr == 0) ? sr0 : sr0 + dsr*2*(rand.nextDouble() - 0.5);
	}
	
	private void TakePicture(Grid grid){
		
		String SaveAs = PicDir + File.separator + grid.getTitle()+fmtL.format(mct)+"_"+fmtS.format(ln)+".png";

		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
	private double[] getRho(){
		
		double[] ret = new double[N];
		int[] dnbs; 
		
		for (int jj = 0 ; jj < N ; jj++){
			if (alive[jj]) continue;
			dnbs = nbs.get(jj);
			for (int kk = 0 ; kk < dnbs.length ; kk++){
				ret[dnbs[kk]]++;
			}
		}

		return ret;
	}

	
}
