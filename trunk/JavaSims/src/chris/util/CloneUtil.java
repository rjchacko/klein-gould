package chris.util;

import scikit.jobs.params.Parameters;
import chris.ofcdamage.damage2Dfast;
import chris.ofcdamage.ofc2Dfast;

public class CloneUtil {
	
	public static ofc2Dfast cloneOFC(Parameters params, ofc2Dfast cp, boolean cloneRand){
		
		ofc2Dfast ret = new ofc2Dfast(dummyParamUtil.ofcParams(params), cp.getSr(), cp.getSf(), cp.getStress(), cp.getSbar(), cp.getData(), cp.getGR(), cp.getLastShower());
		if(cloneRand) ret.setRandToClone(cp.getRand());

		return ret;
	}

	public static damage2Dfast cloneDamage(Parameters params, damage2Dfast cp, boolean cloneRand){
		
		damage2Dfast ret = new damage2Dfast(dummyParamUtil.ofcParams(params), cp.getSr(), cp.getSf(), cp.getStress(), cp.getSbar(), cp.getData(), cp.getGR(), cp.getLastShower(), cp.getFailed(), cp.Ndead(), cp.getLN(), cp.getLives());
		if(cloneRand) ret.setRandToClone(cp.getRand());

		return ret;
	}

}
