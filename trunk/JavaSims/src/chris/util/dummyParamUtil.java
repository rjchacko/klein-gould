package chris.util;

import scikit.jobs.params.Parameters;

public class dummyParamUtil {
	
	public static Parameters ofcParams(Parameters params){
		
		/*
		 * 
		 * List of all parameters in ofc2Dfast or classes that extend it:
		 * 
		 *  ("Data Directory");
		 *	("Data File");
		 *	("Random Seed");
		 *	("Interaction Shape");
		 *	("Interaction Radius (R)");
		 *	("Lattice Size");
		 *	("Boundary Condtions");
		 *	("Failure Stress (\u03C3_f)");
		 *	("\u03C3_f width");
		 *	("Residual Stress (\u03C3_r)");
		 *	("\u03C3_r width");
		 *	("Dissipation (\u03B1)");
		 *	("\u03B1 width");	
		 *	("L_cg (must be odd)");
		 *	("t_cg");
		 *	("Heal Time");
		 *	("HT Width");	
		 *	("Number of Lives");
		 *	("NL width");
		 *  ("\u03D5");
		 *
		 *
		 */
		
		Parameters dparms = params;
		
		String[] pstr = {"Data Directory","Data File","Random Seed","Interaction Shape","Interaction Radius (R)","Lattice Size",
						 "Boundary Condtions", "Failure Stress (\u03C3_f)", "\u03C3_f width", "Residual Stress (\u03C3_r)", 
						 "\u03C3_r width", "Dissipation (\u03B1)", "\u03B1 width", "L_cg (must be odd)", 
						 "t_cg", "Heal Time", "HT Width", "Number of Lives", "NL width", "\u03D5"};
		
		
		for (int jj = 0 ; jj < 20 ; jj++){
			if(params.containsKey(pstr[jj])) continue;
			
			switch (jj) {
	
				case 0:
					dparms.add(pstr[0],"/Users/cserino/Documents/BU/Research/default_location");
					continue;

				case 1:
					dparms.add(pstr[1],"default");
					continue;
					
				case 2:
					dparms.add(pstr[2],0);
					continue;
					
				case 3:
					dparms.add(pstr[3],"Circle");
					continue;

				case 4:
					dparms.add(pstr[4],1);
					continue;
					
				case 5:
					dparms.add(pstr[5],10);
					continue;			

				case 6:
					dparms.add(pstr[6],"Open");
					continue;

				case 7:
					dparms.add(pstr[7],2.);
					continue;
						
				case 8:
					dparms.add(pstr[8],0.);
					continue;			
				
				case 9:
					dparms.add(pstr[9],1.);
					continue;

				case 10:
					dparms.add(pstr[10],0.);
					continue;

				case 11:
					dparms.add(pstr[11],0.05);
					continue;		
					
				case 12:
					dparms.add(pstr[12],0.);
					continue;

				case 13:
					dparms.add(pstr[13],5);
					continue;

				case 14:
					dparms.add(pstr[14],10);
					continue;			
				case 15:
					dparms.add(pstr[15],1);
					continue;

				case 16:
					dparms.add(pstr[16],0);
					continue;

				case 17:
					dparms.add(pstr[17],1);
					continue;
					
				case 18:
					dparms.add(pstr[18],0);
					continue;

				case 19:
					dparms.add(pstr[19],1.);
					continue;
					
				default:
					throw new IllegalArgumentException("Index overflow");
			}	
		}

		return dparms;
	}

}
