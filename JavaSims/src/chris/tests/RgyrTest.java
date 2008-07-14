package chris.tests;



import scikit.dataset.Histogram;
import scikit.graphics.ColorPalette;
import scikit.graphics.dim2.Grid;
import scikit.graphics.dim2.Plot;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.ChoiceValue;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.DoubleValue;
import chris.ofc.old.NfailDamage2D;


public class RgyrTest extends Simulation{

	Grid grid1 = new Grid ("Sites");
	Plot plot1 = new Plot("R_gyr");
	Plot plot2 = new Plot("Alive");
	double chk;
	NfailDamage2D model;
	Histogram histRgyr;
	int loopvar;
	double MI;
	
	
	ColorPalette palette1;
	
	public static void main(String[] args) {
		new Control(new RgyrTest(), "OFC Model");
	}
	
	public void load(Control c) {		

		params.add("Data Directory",new DirectoryValue("/Users/cserino/CurrentSemester/Research/"));
		params.add("Random Seed",0);
		params.addm("Auto Scale", new ChoiceValue("Yes", "No"));
		params.add("Lattice Size",1<<9);
		params.add("Number of Lives",1);
		params.add("Life Style", new ChoiceValue("Constant","Flat","Gaussian"));
		params.add("Nlives Width",0.1);
		params.add("Boundary Condtions", new ChoiceValue("Periodic","Bordered"));
		params.add("Stress Distribution", new ChoiceValue("Flat","Hammer Blow"));
		params.add("Hammer Size",1);	
		params.add("Critical Stress (\u03C3_c)",4.0);
		params.add("\u03C3_c Noise", new ChoiceValue("Off","On"));	
		params.add("\u03C3_c width",Math.sqrt(Math.sqrt(0.4)));
		params.add("Residual Stress (\u03C3_r)",2.0);
		params.add("\u03C3_r Noise", new ChoiceValue("Off","On"));
		params.add("\u03C3_r width",Math.sqrt(Math.sqrt(2)));
		params.add("Interaction Shape", new ChoiceValue("Circle","Square","Diamond"));
		params.add("Interaction Radius (R)",(int)(50));
		params.add("Minimum Interaction Radius (r)",30);
		params.add("Dissipation (\u03B1)",new DoubleValue(0.2,0,1));
		params.add("\u03B1 Noise", new ChoiceValue("On","Off"));
		params.add("\u03B1 Width", 0.05);
		params.add("Number of Resets");
		params.add("Number of Showers");
			
		c.frame(grid1);
	}
	
	public void animate() {
		
		int[] foo = new int[model.N];
		
		for (int i = 0 ; i < model.N ; i++){
			foo[i]=model.alive[i+model.N];
		}
		
		grid1.registerData(model.L,model.L,foo);

		params.set("Number of Resets",MI);

	}

	public void clear() {
		grid1.clear();
		plot1.clear();
	}

	public void run() {
//		
//		int ic;
//	
//		model = new NfailDamage2D(params);
//
//		model.Initialize("Flat");
//
//		if(model.L%2 == 0){
//			ic=(model.N+model.L)/2;
//		}
//		else{
//			ic=model.N/2+1;
//		}
//		
//		palette1 = new ColorPalette();
//		palette1.setColor(0,Color.BLACK);
//		palette1.setColor(1,Color.WHITE);
//		grid1.setColors(palette1);
//		
//		for (int i = 0 ; i < 2*model.N ; i++){
//			model.alive[i]=0;
//		}
	
//
//		int[] cc = model.neighbors.get(ic);
//		
//		for (int i = 0 ; i < cc.length ; i++){
//			model.alive[cc[i]+model.N]=1;
//		}
//		
//		MI=model.radiusGyration(ic);
//		
//		MI=MI*MI;
		
		Job.animate();
		
	}	
		
			
}


//if (BCs.equals("Bordered")){
//	if(shape.equals("Circle")){
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Circle);
//	}
//	else if(shape.equals("Square")){
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Square);
//	}
//	else{
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.BORDERED,LatticeNeighbors.Shape.Diamond);
//	}
//}
//else{
//	if(shape.equals("Circle")){
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Circle);
//	}
//	else if(shape.equals("Square")){
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Square);
//	}
//	else{
//		neighbors = new LatticeNeighbors(L,L,rmin,R,LatticeNeighbors.Type.PERIODIC,LatticeNeighbors.Shape.Diamond);
//	}
//}

