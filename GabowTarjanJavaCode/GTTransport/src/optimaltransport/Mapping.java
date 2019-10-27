package optimaltransport;

public class Mapping {
	private long startTime;

	/**
	 * Computes an additive approximation of the optimal transport between two discrete probability distributions, A and B, 
	 * using the Gabow-Tarjan transportation problem algorithm as a subroutine.
	 * We say that the set A contains the 'demand' vertices and the set B contains the 'supply' vertices.
	 * @param n : The number of supply and demand vertices, i.e., n = |A| = |B|.
	 * @param supplies : The probability distribution associated with the supply locations. 
	 * @param demands : The probability distribution associated with the supply locations. 
	 * Requirements: Both supplies and demands must be size n and sum to 1.
	 * @param C : An n x n cost matrix, where C[i][j] gives the cost of transporting 1 unit of mass from the ith demand vertex
	 * to the jth supply vertex.
	 * 
	 * Computes a transport plan with additive error at most delta away from the optimal and stores the transport in the
	 * 'flow' variable, which can be retieved afterwards using the 'getFlow()' method.
	 */
	public Mapping(int n, double[] supplies, double[] demands, double[][] C, double delta) {
		startTime = System.currentTimeMillis();
		
		double max = 0;
		
		for (int i = 0; i < n; i++) {
			max = Double.max(max, supplies[i]);
		}
		
		double maxCost = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				maxCost = Double.max(maxCost, C[i][j]);
			}
		}
		
		//Convert the inputs into an instance of the transportation problem with integer supplies, demands, and costs.
		int[][] scaledC = new int[n][n];
		int[] scaledDemands = new int[n];
		int[] scaledSupplies = new int[n];
		double alpha = 4.*maxCost / delta;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				scaledC[i][j] = (int)(C[i][j] * alpha);
			}
			scaledDemands[i] = (int)(Math.ceil(demands[i] * alpha * n));
			scaledSupplies[i] = (int)(supplies[i] * alpha * n);
		}
	
		//Call the main Gabow-Tarjan algorithm routine to solve the scaled instance.
		//Returns a maximum-size transport plan additive error at most sum(scaledSupplies).
		GTTransport gt = new GTTransport(scaledC, scaledSupplies, scaledDemands, n);
		
		//Record the efficiency-related results of the main routine
		this.iterations = gt.getIterations();
		this.APLengths = gt.APLengths();
		this.mainRoutineTimeInSeconds = gt.timeTakenInSeconds();
		this.timeTakenAugment = gt.timeTakenAugmentInSeconds();
		
		//Get the solution for the scaled instance.
		int[][] scaledFlow = gt.getSolution();
	
		flow = new double[n][n];
		
		//Scale back flows and compute residual (leftover) supplies and demands.
		double[] residualSupply = new double[n];
		double[] residualDemand = new double[n];

		for (int i = 0; i < n; i++) {
			residualSupply[i] = supplies[i];
			residualDemand[i] = demands[i];
		}

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				flow[i][j] = scaledFlow[i + n][j] / (n * alpha); 
				residualSupply[j] -= flow[i][j];
				residualDemand[i] -= flow[i][j];
			}
		}
		
		//Push back some flow incoming to demand constraints that are violated.
		for (int j = 0; j < n; j++) {
			for (int i = 0; residualDemand[j] < 0 && i < n; i++) {
				double reduction = Double.min(-residualDemand[j], flow[j][i]);
				flow[j][i] -= reduction;
				residualDemand[j] += reduction;
				residualSupply[i] += reduction;
			}
		}

		//Arbitrarily match the remaining supplies
		for (int i = 0; i < n; i++) {
			for (int j = 0; residualSupply[i] > 0 && j < n; j++) {
				double increase = Double.min(residualSupply[i], residualDemand[j]);
				flow[j][i] += increase;
				residualDemand[j] -= increase;
				residualSupply[i] -= increase;
			}
		}
		
		this.totalCost = 0;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				totalCost += flow[i][j] * C[i][j];
			}
		}
		
		timeTaken = (System.currentTimeMillis() - startTime)/1000.;
	}
	
	//A method that can be called to verify the contents of the produced flow to ensure it is feasible.
	public void verifyFlow(double[] supplies, double[] demands, double[][] flow) {	
		int n = supplies.length;
		for (int i = 0; i < n; i++) {
			double sumB = 0;//sum for flow outgoing from supply vertex i
			double sumA = 0;//sum for flow incoming to demand vertex i
			for (int j = 0; j < n; j++) {
				sumB += flow[j][i];
				sumA += flow[i][j];
			}
			double residualB = supplies[i] - sumB;
			double residualA = demands[i] - sumA;
			double threshold = 0.00001;
			if (Math.abs(residualB) > threshold) {
				System.err.println("Violation B: " + residualB + " at index " + i);
			}
			if (Math.abs(residualA) > threshold) {
				System.err.println("Violation A: " + residualA + " at index " + i);
			}
		}
	}
	
	private double mainRoutineTimeInSeconds;
	
	public double getMainRoutineTimeInSeconds() {
		return mainRoutineTimeInSeconds;
	}

	private double[][] flow;
	private double totalCost;
	private double timeTaken;
	private double timeTakenAugment;
	private int iterations;
	private int APLengths;
	

	public double[][] getFlow() {
		return flow;
	}


	public double getTotalCost() {
		return totalCost;
	}


	public double getTimeTaken() {
		return timeTaken;
	}


	public int getIterations() {
		return iterations;
	}


	public int getAPLengths() {
		return APLengths;
	}
	
	public double getTimeTakenAugment() {
		return timeTakenAugment;
	}
}
