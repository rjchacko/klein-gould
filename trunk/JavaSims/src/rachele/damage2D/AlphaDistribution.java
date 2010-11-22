package rachele.damage2D;

import rachele.damage2D.multidx.Damage;
import kip.util.Random;

public class AlphaDistribution {
	static int pow;
	static int L;
	static int N;
	static int R;
	static Random random = new Random();
	static String infoFileName;
	static double [] alpha;

	public AlphaDistribution(int power, int range, int randomSeed, String fileName) {
		pow=power;
		L=(int) Math.pow(2, pow);
		N=L*L;
		System.out.println("N = " + N);
		R=range;
		alpha = new double [N];
		setRandomSeed(randomSeed);
		infoFileName=fileName;


	}

	static void setRandomSeed(int rs){
		random.setSeed(rs);	
	}
	
	static public double [] setAlphaArray(String alphaDistribution, double alphaParam){

		System.out.println("Alpha dist = " + alphaDistribution);
		if(alphaDistribution=="FlatRandom"){
			for (int i = 0; i < N; i++) alpha[i]=random.nextDouble();
		}else if(alphaDistribution=="Gaussian"){
			for (int i = 0; i < N; i++){
				double a =random.nextGaussian()/10 + alphaParam;
				if(a > 1.0) alpha[i] = 1.0;
				else if (a < 0) alpha[i] = 0.0;
				else alpha[i] = a;
			}
		}else if(alphaDistribution=="Constant"){
			for (int i = 0; i < N; i++) alpha[i]=alphaParam;
		}else if(alphaDistribution=="Gaussian about zero"){
			for (int i = 0; i < N; i++){
				double a =Math.abs(random.nextGaussian()/5.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/5.0);
				}
				alpha[i] = Math.abs(a);
			}
		}else if(alphaDistribution=="GaussianSplit"){
			for (int i = 0; i < N; i++){
				double a =(random.nextGaussian()/10.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/10.0);
				}
				if(a < 0.0) a = 1.0 + a;
				alpha[i] = a;
			}
		}else if(alphaDistribution=="GaussianAboutHalf"){
			for (int i = 0; i < N; i++){
				double a =(random.nextGaussian()/10.0);
				while(Math.abs(a)>0.5){
					a =(random.nextGaussian()/10.0);
				}
				a += 0.5;
				alpha[i] = a;
			}
		}else if(alphaDistribution=="DeadBlocks"){
			boolean [] aa = Damage.setDamage("Place Dead Blocks", 32,0.0, 32);
			double a;
			for(int i = 0; i < N; i++){
				a=(random.nextGaussian()/10.0);
				while(Math.abs(a)>1.0){
					a =(random.nextGaussian()/10.0);
				}
				if(aa[i]){
					alpha[i] = 1.0- Math.abs(a);
				}else{
					alpha[i] = Math.abs(a);
				}
			}	
		}else if(alphaDistribution=="ManyGaussians"){
			System.out.println("many Gaussians");
			int dx = 16;
			int Lp = L/dx;
			int Np = Lp*Lp;
			boolean [] blockSet = new boolean [Np];
			for(int i = 0; i < Np; i ++) blockSet[i] = false;
			for(int i = 0; i < Np; i ++){
				//choose a block
				int block = (int)(random.nextDouble()*(double)(Np));
				while(blockSet[block]){
					block = (int)(random.nextDouble()*(double)(Np));
				}
				double center = ((double)i+0.5)/(double)(Np);
				setGaussianBlock(center, block, dx, Np);
				blockSet[block] = true;
			}

		}else if(alphaDistribution=="Fractal"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;

			int maxPow = pow-2;  //largest lattice size is 8 x 8
			int [] intervalCount = new int [maxPow+1];
			for (int ip = maxPow; ip >= 0; ip--){
				int dx = (int)Math.pow(2 ,ip);
				int Lp = L/dx;
				int Np = Lp*Lp;
				double center = 1.0-((double)ip+0.5)/(double)(maxPow+1);
				for (int i = 0; i < Np; i++){
					double pr = 0.25;  						//Set one quarter of blocks to alpha distributed about center
					if(alpha[getFirstBlockSite(i, ip)]<0){
						if (random.nextDouble()<pr){
							setGaussianBlock(center, i, dx, maxPow+1);
							intervalCount[ip]+=dx*dx;
						}
					}
				}	
			}
			int [] unsetList = randomizeUnsetSites();
			for (int i = 0; i < unsetList.length; i++){
				int minI = findMinInterval(intervalCount);
				double center = 1.0-((double)minI+0.5)/(double)(maxPow+1);
				setGaussianBlock(center, unsetList[i], 1, maxPow+1);
				intervalCount[minI]+=1;
			}
		}else if (alphaDistribution=="Quarters"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;

			//set one block 1/4 the lattice size
			int ip = pow - 1;
			int dx = (int)Math.pow(2 ,ip);
			double center = 0.5/4.0;
			setGaussianBlock(center, 1, dx, 4);
			//set 4 blocks 1/8 the size of lattice
			ip = pow - 2;
			dx = (int)Math.pow(2 ,ip);
			int Lp = L/dx;
			int Np = Lp*Lp;
			center = 1.5/4.0;
			int noSet = 0;
			while(noSet < 4){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			//set 16
			ip = pow - 3;
			dx = (int)Math.pow(2 ,ip);
			Lp = L/dx;
			Np = Lp*Lp;
			center = 2.5/4.0;
			noSet = 0;
			while(noSet < 16){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			//set 16 more blocks
			center = 3.5/4.0;
			noSet = 0;
			while(noSet < 16){
				//pick a block of size dx
				int randBlock = (int)(random.nextDouble()*Np);
				if(alpha[getFirstBlockSite(randBlock, ip)]<0){
					setGaussianBlock(center, randBlock, dx, 4);
					noSet+=1;
					System.out.println("Set block of size " + dx);
				}
			}
			int [] unsetList = randomizeUnsetSites();
			System.out.println("no unset = " + unsetList.length);
		}else if (alphaDistribution=="Eights"){
			for (int i = 0; i < N; i++)
				alpha[i]=-1.0;

			int maxPow = pow-2;  //largest lattice size is 4 X 4
			for (int ip = maxPow; ip > 0; ip--){
				int dx = (int)Math.pow(2 ,ip);
				int Lp = L/dx;
				int Np = Lp*Lp;
				double center = (7.5-(double)ip)/8.0;
				int noSet = 0;
				int maxToSet = Np/8;
				while(noSet < maxToSet){
					//pick a block of size dx
					int randBlock = (int)(random.nextDouble()*Np);
					if(alpha[getFirstBlockSite(randBlock, ip)]<0){
						setFlatBlock(center, randBlock, dx, 8);
						noSet+=1;
						System.out.println("Set block of size " + dx);
					}
				}

			}
			double center = 7.5/8.0;
			for (int i = 0; i < N; i++){
				if(alpha[i]< 0){
					setFlatBlock(center, i, 1, 8);
				}
			}
		}else{
			System.out.println("No Gaussian Dist");
		}
		
		return alpha;
	}


	static int [] randomizeUnsetSites(){
		int ct = 0;
		for (int i = 0; i < N; i++){
			if(alpha[i]<0){
				ct += 1;
			}
		}
		int [] list = new int [ct];
		for (int i = 0; i < ct; i++){
			list[i] = -1;
		}
		for (int i = 0; i < N; i++){
			if(alpha[i]<0){
				int index = (int)(random.nextDouble()*ct);
				while(list[index]>0){
					index = (int)(random.nextDouble()*ct);
				}
				list[index]=i;
			}
		}
		return list;
	}

	static int findMinInterval(int [] a){
		int minCt = N;
		int minI = -1;
		for (int i = 0; i < a.length; i++){
			if (a[i] < minCt){
				minCt = a[i];
				minI = i;
			}
		}
		return minI;
	}

	static void setGaussianBlock(double center, int block, int dx, int noInts){
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				double rand = random.nextGaussian();
				while(Math.abs(rand)>=1.0){
					rand = random.nextGaussian();
				}
				alpha[site] = rand/(double)(noInts*2)+center;
			}
		}
	}

	static void setFlatBlock(double center, int block, int dx, int noInts){
		int Lp = L/dx;
		int xp = block%Lp;
		int yp = block/Lp;
		double intervalWidth = 1.0/noInts;
		for (int y = yp*dx; y < yp*dx+dx; y++){
			for (int x = xp*dx; x < xp*dx + dx; x++){
				int site = y*L+x; 
				double rand = random.nextDouble()*intervalWidth;

				alpha[site] = rand+center-intervalWidth/2.0;
			}
		}
	}

	static int getFirstBlockSite(int blockNo, int ip){
		int dx = (int)Math.pow(2, ip);
		int Lp = L/dx;
		int xp = blockNo%Lp;
		int yp = blockNo/Lp;
		int x = dx*xp;
		int y = dx*yp;
		int site = y*L+x;
		return site;
	}

}
