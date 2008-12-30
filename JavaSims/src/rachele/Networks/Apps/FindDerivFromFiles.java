package rachele.Networks.Apps;

//not working.  Use deriv.f90 instead.
// Problem with java interpretation of characters.
//Fortran handles the files much better.

import java.io.File;

import rachele.util.FileUtil;
import scikit.jobs.Control;
import scikit.jobs.Simulation;
import scikit.jobs.params.DirectoryValue;
import scikit.jobs.params.FileValue;

/**
* 
* Used to calculate d log(delta m)/d tau where tau \equiv (h - h_sp)/h_sp.
* 
* Reads in results from two runs: lower run (h = h_sp - \epsilon) and higher
* run (h = h_sp + \epsilon.)
* 
* Calculates the derivative using 2 step difference:
* 
* d log(delta m)/d tau = (log(delta m)( h + \epsilon ) - log(delta m)(h - \epsilon)) / ( 2 * \epsilon )
*  
*/
public class FindDerivFromFiles extends Simulation{
	
	public static void main(String[] args) {
		new Control(new FindDerivFromFiles(), "Find derivative from files");
	}
	

	public void animate() {
	}


	public void clear() {
	}


	public void load(Control c) {
		params.add("Data File high",new FileValue("/home/erdomi/data/spinodal_find/cw/testRuns/test"));
//		params.add("Data File high",new FileValue("/home/erdomi/data/spinodal_find/cw/t49/h0.3167M"));
		params.add("Data File low",new FileValue("/home/erdomi/data/spinodal_find/cw/testRuns/test"));
//		params.add("Data File low",new FileValue("/home/erdomi/data/spinodal_find/cw/t49/h0.3184M"));
		params.add("Output File",new DirectoryValue("/home/erdomi/data/spinodal_find/cw/t49"));
		params.addm("epsilon", 1.0e-4);
		
	}


	public void run() {
		double epsilon = params.fget("epsilon");
		String inFileLo = params.sget("Data File low");
		String inFileHi = params.sget("Data File high");
		double [][] dataLo = FileUtil.readDoubleData(inFileLo);
		System.out.println("done lo");
		double [][] dataHi = FileUtil.readDoubleData(inFileHi);
		System.out.println("done hi");
		int size = Math.min(dataLo[0].length, dataHi[0].length);
		double [] deriv = new double [size];
		for (int i = 0; i < size; i++){
			if(dataLo[0][i] == dataHi[0][i]){
			deriv[i] = (Math.log(dataHi[1][i])-Math.log(dataLo[1][i]))/(2*epsilon);
			}else{
				System.out.print("Time mismatch: Lo=" + "dataLo[0][i]" + " Hi = " + dataHi[0][i]);
			}
		}
		String message1 = "Derivative: d log(delta m)/d tau where tau equiv (h - h_sp)/h_sp.";
		String outFile = params.sget("Output File") + File.separator + "deriv";
		FileUtil.initFile(outFile, params, message1);
		for (int i = 0; i < size; i++){
			FileUtil.printlnToFile(outFile, dataLo[0][i], deriv[i]);
		}
		
		
	}

	
}
