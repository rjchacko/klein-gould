package chris.old.ofc.apps;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Random;

import javax.imageio.ImageIO;

import chris.util.DirUtil;
import chris.util.LatticeNeighbors;
import chris.util.PrintUtil;

import scikit.graphics.ColorGradient;
import scikit.graphics.dim2.Grid;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;

public class OFCfailsafe extends Simulation{

	private double stress[], sbar[], sr[], sr0, dsr, sf, alpha, metric[], ss[], ds, dsa[];
	private int R, L, N, eqtime, simtime, search[], qN, GRcount, GR[];
	private String fout, imode, bname, BCs, PicDir; // NNdebug;
	private Random rand;         
	private boolean draw, record;
	private LatticeNeighbors nbs;
	private Grid gridS = new Grid("Stress");
	private ColorGradient cg = new ColorGradient();
	private DecimalFormat fmt = new DecimalFormat("0000000");
	
	private static int dlength = 10;


	public static void main(String[] args) {
		new Control(new OFCfailsafe(), "OFC Parameters");
	}
	
	public void animate() {
		
		if(draw){
			gridS.registerData(L, L, stress);
		}
		
	}

	public void clear() {
		if(draw){
			gridS.clear();
		}
	}

	public void load(Control c) {
		
		params.add("Data Directory",new DirectoryValue("/Users/cserino/Desktop/"));
		params.add("File Name","FailSafe");
		params.add("Initialization Mode", new ChoiceValue("Random","Site Tracking"));
		params.add("Random Seed",0);
		params.add("Boundary Conditions", new ChoiceValue("Peroidic", "Open"));
		params.add("Interaction Radius (R)",(int)(30));
		params.add("Lattice Size",1<<8);
		params.add("Equilibrate Time",(int)400000);
		params.add("Simulation Time",(int)400000);
		params.add("Failure Stress (\u03C3_f)",2.0);
		params.add("Residual Stress (\u03C3_r)",1.);
		params.add("\u03C3_f width",(double) 0);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.01,0,1));
		params.add("Animation", new ChoiceValue("Off","On"));
		params.addm("Record", new ChoiceValue("Off","On"));
		params.add("Status");
		
		c.frame(gridS);
		
	}

	public void run() {
		
		Initialize();
		
		// equilibrate
		draw = false;
		for (int jj = 0 ; jj < eqtime ; jj++){
			Evolve();
			if(jj%200 == 0) params.set("Status",jj - eqtime);
			Job.animate();
		}
		
		// simulate
		draw = (params.sget("Animation") == "On");
		for (int jj = 0 ; jj < simtime ; jj++){
			Evolve();
			saveData(jj);
			if(jj%200 == 0) params.set("Status", jj);
			Job.animate();
			if(record) TakePic(gridS,jj);
		}
		
		// write data
		writeData(simtime);
		params.set("Status","Done");
		Job.animate();
		
	}
	
	private void Initialize() {
		
		int tmp;
		
		fout    = params.sget("Data Directory");
		bname   = params.sget("File Name");
		imode   = params.sget("Initialization Mode");
		tmp     = params.iget("Random Seed");
		BCs     = params.sget("Boundary Conditions");
		R       = params.iget("Interaction Radius (R)");
		L       = params.iget("Lattice Size");
		eqtime  = params.iget("Equilibrate Time");
		simtime = params.iget("Simulation Time");
		sf     = params.fget("Failure Stress (\u03C3_f)");
		sr0	    = params.fget("Residual Stress (\u03C3_r)");
		dsr     = params.fget("\u03C3_f width");
		alpha   = params.fget("Dissipation (\u03B1)");
		draw    = (params.sget("Animation") == "On");
		record  = (params.sget("Record") == "On");
		
		String tmps = fout + File.separator + "Params_" + bname + ".txt";
	//	NNdebug = fout + File.separator + "NNdebug" + ".txt";
		
		PrintUtil.printlnToFile(tmps, params.toString());

		N    = L*L;
		fout += File.separator + bname + ".txt";
		rand = new Random(tmp);
		
		sbar   = new double[N];
		stress = new double[N];
		sr     = new double[N];
		search = new int[N];
		metric = new double[dlength];
		ss     = new double[dlength];
		GR     = new int[dlength];
		dsa    = new double[dlength];
		
		for (int jj = 0 ; jj < N ; jj++){
			sbar[jj] = 0;
			if(imode == "Random"){
				stress[jj] = sr0 + (sf - sr0)*rand.nextDouble();
			    sr[jj]     = sf + dsr*2*(rand.nextDouble() - 0.5);
			}
			else{
				stress[jj] = sf;
				sr[jj]     = sr0;
				stress[(int)(N/2) + (int)(L/2)] = 1.5*sr0;
			}		
			search[jj] = 0;
		}
		
		
		if(BCs == "Peroidic"){
			nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
			qN  = (nbs.get(1)).length;
		}
		else if(BCs == "Open"){
			nbs = new LatticeNeighbors((int) L,(int) L,0,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
			int tx = L/2;
			qN  = (nbs.get(L*tx + tx)).length;
		}
		PrintUtil.printlnToFile(tmps,"Neighbors = ", qN);

			
		if(draw){
			double clr;
			
			gridS.setScale(sr0-dsr,sf);
			for (int jj = 0 ; jj <= N ; jj++){
				clr = sr0-dsr + jj*(sf-sr0+dsr)/N;
				cg.getColor(clr,sr0-dsr,sf);
			}
			gridS.setColors(cg);
	
			if(record){
				PicDir = params.sget("Data Directory") + "/GridPics"; 
				DirUtil.MkDir(PicDir);
			}
		}
			
			
			
		return;
	}
	
	private void Evolve(){
		
		int jmax = 0;
		int counter = 0;
		double dds;
		int nlist[], nlength;
		double release;
		ds = 0;
	//	double smxb = 0;
	//	double smxa = 0;
		
		// force the first failure
		for (int jj = 0 ; jj < N ; jj++){
			if(stress[jj]>stress[jmax]) jmax = jj;
		}
		dds = sf-stress[jmax];
		for (int jj = 0 ; jj < N ; jj++){
			stress[jj] += dds;
		}
		
		search[counter++] = jmax;
		GRcount = 1;
		
		while(counter > 0){
			
			// distribute stress
			for (int jj = 0 ; jj < counter ; jj++){
				nlist   = nbs.get(search[jj]);
				nlength = nlist.length;
				release = (1-alpha)*(stress[search[jj]] - sr[search[jj]]) / qN;
				for (int kk = 0 ; kk < nlength ; kk++){
					//if(stress[nlist[kk]]>smxb) smxb = stress[nlist[kk]];
					stress[nlist[kk]] += release; 
					//if(stress[nlist[kk]]>smxa) smxa = stress[nlist[kk]];
				}
				//if(draw) PrintUtil.printlnToFile(NNdebug,smxb,smxa);
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
					if(stress[jj]-sf > ds) ds = stress[jj]-sf;
				}
			}
			GRcount += counter;
			
		}

		return;
	}
	
	private double nextSr(){
		
		return (dsr == 0) ? sr0 : sr0 + dsr*2*(rand.nextDouble() - 0.5);
	}
	
	private void saveData(double mcStep){
		
		
		if((mcStep)%dlength == 0 && mcStep > 0){
			writeData((int)(mcStep));
		}
		
		double tmpbar, Omega, sumS;
		int index = (int)((mcStep)%dlength);
		
		// calculate metric
		tmpbar = 0;
		sumS   = 0;

		for (int jj = 0 ; jj < N ; jj++){
			sbar[jj] += stress[jj];
			tmpbar   += sbar[jj];
			sumS     += stress[jj];
		}

		tmpbar = tmpbar / N;

		Omega = 0;
		for (int jj = 0 ; jj < N ; jj++){
			Omega += (sbar[jj] - tmpbar)*(sbar[jj] - tmpbar);
		}

		Omega = Omega/((mcStep+1)*(mcStep+1)*N);
		
		// save metric
		metric[index] = Omega;
		ss[index]     = sumS;
		GR[index]     = GRcount;
		dsa[index]    = ds;
		
		return;
	}
	
	
	

	private void writeData(int Nsteps){
		
		int offset = (int)((Nsteps-1)/dlength);
		offset = offset*dlength;
		
		try{
			File file = new File(fout);
			PrintWriter pw = new PrintWriter(new FileWriter(file, true), true);
			if(Nsteps == dlength){
				pw.print("MC Step");
				pw.print("\t");
				pw.print("Metric");
				pw.print("\t");
				pw.print("Total Stress");
				pw.print("\t");
				pw.print("Shower Size");
				pw.print("\t");
				pw.print("Fail Stress");
				pw.println();
			}
			if(Nsteps == simtime){
				for (int jj = 0 ; jj < (Nsteps-offset) ; jj++){		
					pw.print(jj+1+offset);
					pw.print("\t");
					pw.print(metric[jj]);
					pw.print("\t");
					pw.print(ss[jj]);
					pw.print("\t");
					pw.print(GR[jj]);
					pw.print("\t");
					pw.print(dsa[jj]);
					pw.println();
				}
			}
			else{
				for (int jj = 0 ; jj < dlength ; jj++){		
					pw.print(jj+1+offset);
					pw.print("\t");
					pw.print(metric[jj]);
					pw.print("\t");
					pw.print(ss[jj]);
					pw.print("\t");
					pw.print(GR[jj]);
					pw.print("\t");
					pw.print(dsa[jj]);
					pw.println();
				}
			}
			pw.close();
		}
		catch (IOException ex){
			ex.printStackTrace();
		}
		
		return;
	}
	
	private void TakePic(Grid grid, int mcStep){
		
		String SaveAs = PicDir + File.separator + grid.getTitle()+fmt.format(mcStep)+".png";
		
		try {
			ImageIO.write(grid.getImage(), "png", new File(SaveAs));
		} catch (IOException e) {
			System.err.println("Error in Writing File" + SaveAs);
		}
		
		return;
	}
	
}
