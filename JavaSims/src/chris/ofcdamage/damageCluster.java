package chris.ofcdamage;

import java.io.File;

import chris.foo.ofc.clusters.ClustersV4;
import chris.util.PrintUtil;
import scikit.jobs.params.Parameters;

public class damageCluster extends damage{

	
	private ClustersV4 clusters; 
	private boolean perc;
	
	public damageCluster(Parameters params) {
		super(params);
		clusters = new ClustersV4(getL(), getBC());
		perc = false;
	}
	
	protected void add2cluster(int site){
		
		perc = clusters.addSite(site);
		return;
	}
	
	public int[] getClusters(){

		int N     = getL()*getL();
		int[] ret = new int[N];

		for (int jj = 0 ; jj < N ; jj++){
			ret[jj] = clusters.getClusterNumber(jj);
		}

		return ret;
	}
	
	protected boolean sysFail(){
		
		return perc;
	}
	
	public void printCdistance(int cycle){
		
		int N    = getL()*getL();
		int[] dX = new int[N];
		int[] dY = new int[N];
		
		for(int jj = 0 ; jj < N ; jj++){
			dX[jj] = clusters.getDist(jj);
		}
		for(int jj = 0 ; jj < N ; jj++){
			dY[jj] = clusters.getDist(jj+N);
		}
		
		PrintUtil.printArrayToFile(getOutdir()+File.separator+"ClustersDX_"+(cycle+1)+".txt", dX,N,1);
		PrintUtil.printArrayToFile(getOutdir()+File.separator+"ClustersDY_"+(cycle+1)+".txt", dY,N,1);
		
		return;
	}
	
	public int pcnPC(){
		
		return clusters.whichCluster();
	}
	
	public int getSize(int cn){
		
		return clusters.getClusterSize(cn);
	}
	
	public int getSizeV2(int cn){
		
		return clusters.getClusterSizeV2(cn);
	}

}
