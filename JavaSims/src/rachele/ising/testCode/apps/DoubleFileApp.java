package rachele.ising.testCode.apps;

import rachele.util.FileUtil;
import scikit.jobs.Control;
import scikit.jobs.Job;
import scikit.jobs.Simulation;
import scikit.jobs.params.FileValue;

public class DoubleFileApp extends Simulation{

	public void animate() {

		
	}

	public static void main(String[] args) {
		new Control(new DoubleFileApp(), "Double Phi0 File");
	}
	
	public void clear() {
		
	}

	public void load(Control c) {
		params.add("Input Config",new FileValue("/home/erdomi/data/lraim/configs1d/config"));
		params.add("Output Config",new FileValue("/home/erdomi/data/lraim/configs1d/config2"));
		params.addm("Lp value", 128);
	}

	public void run() {
		int Lp = params.iget("Lp value");
		double [] phi = new double [Lp];
		double [] phiDouble = new double [Lp*2];
		phi = FileUtil.readConfigFromFile(params.sget("Input Config"), Lp);
		for (int i = 0; i < Lp; i++){
			phiDouble[i] = phi [i];
			phiDouble[Lp + i] = phi [i];
		}
		String outFile = params.sget("Output Config");
		FileUtil.deleteFile(outFile);
		FileUtil.writeConfigToFile(outFile, 2*Lp, phiDouble);
		System.out.println("output config length is " + 2*Lp);
		Job.animate();
	}

}
