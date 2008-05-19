package ranjit.ofc;


import java.io.IOException;


import scikit.jobs.Control;
import scikit.jobs.Simulation;

public class OFCApp extends Simulation{
	int size = 256;
    int range = 20;
    int binnumber = 20000;
    boolean boundary= true;
    double R=(double) (2*range+1)*(2*range+1)-1;
    double alpha= 0.99*(1.0/R);
    Site Initiator;
    int avalancheSize;
	 
	public static void main(String[] args) throws IOException {
		new Control(new OFCApp(),"OFC model");	    
	}

	public void load(Control c){
		params.add("L",256);
		params.add("R",20);
		params.add("number of bins",20000);
		params.add("boundary conditions", "open");
	}
	@Override
	public void animate() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void clear() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void run() {
		Lattice OFCModel = new Lattice(size, range, binnumber, boundary, alpha); 
	    
	    OFCModel.initializeLattice();
	    
	    System.out.println("Lattice created.");
	    
		for(int i=0;i<1000000;i++){	    
	    	OFCModel.iteration=i;
	    	Initiator=OFCModel.findInitiator();
	    	OFCModel.oldthreshold=OFCModel.maxStress;
	    	OFCModel.maxStress=Initiator.stress;
	    	OFCModel.minStress=OFCModel.maxStress-1;
	    	OFCModel.maxBin=(int) (Math.floor(binnumber*OFCModel.maxStress)%binnumber);
	    
	   	 	avalancheSize=OFCModel.Avalanche(Initiator);
	    	System.out.println( i + "  " + avalancheSize );
		}
		
		
	    for(int i=0;i<size*size;i++){
		System.out.println(OFCModel.site1[i].bin + "," + OFCModel.site1[i].failCounter);
	    }
	}
    
}
