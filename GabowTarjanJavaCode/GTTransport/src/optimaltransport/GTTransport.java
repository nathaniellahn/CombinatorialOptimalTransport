package optimaltransport;

import java.util.ArrayList;
import java.util.Arrays;

public class GTTransport {
	private int[][] CBA;
	private int[][] CAB;

	private int[] y;
	private int[][] capacityAB;
	private int[][] capacityBA;
	private boolean[] bFree;
	private boolean[] aFree;
	private int n;
	private int[] vertexVisited;
	
	private long timeTaken;
	private long timeTakenAugment;
	
	private int iterations;
	private int APLengths;
	
	//Number of calls to findAP
	private int numAPCalls = 0;
	
	// Preallocated array of size 2*n for findAP function.
	private ArrayList<Integer> path; 
	
	private int[][] finalCapacity;
	
	//Note: Does not actually compute the solution. Rather, returns the 
	//solution that was previously computed by the constructor.
	public int[][] getSolution() {
		return finalCapacity;
	}
	
	public double timeTakenInSeconds() {
		return timeTaken / 1000.;
	}
	
	public double timeTakenAugmentInSeconds() {
		return (double)timeTakenAugment / 1E9;
	}
	
	public int getIterations() {
		return iterations;
	}
	
	public int getNumAPCalls() {
		return this.numAPCalls;
	}
	
	public int APLengths() {
		return APLengths;
	}
	
	@Override
	public String toString() {
		return "GTTransport [timeTaken=" + timeTaken + ", iterations=" + iterations + ", APLengths=" + APLengths + "]";
	}
	
	
	//Solves instance of problem in constructor.
	//Solution is then recovered by calling the above "get solution" method.
	//Other relevant values can be recovered using their respective getter methods.
	//Solving problem in constructor guarantees any input instance is only
	//solved once.
	public GTTransport(int[][] C, int[] supplies, int[] demands, int n) {
		long startTime = System.currentTimeMillis();
		this.n = n;
		
		//Convention: CBA contains costs for edges directed from B to A.
		//CBA[i][j], 0 <= i,j < n, gives the cost of the edge from the
		//ith vertex of B to the jth vertex of A.
		//CAB[j][i] is the cost from the ith vertex of A to the jth vertex of B.
		//These costs are symmetric, but storing both is useful for efficiency reasons.
		//We follow this convention also for slacks, etc.
		//Note that this convention is the opposite of that used for the original Matlab implementation of this code.
		//The difference is because of how both languages store arrays (column vs. row order)
		CBA = transpose(C);
		CAB = C;		
		
		//Whether all supplies or demands of a vertex have been satisfied.
		bFree = new boolean[n];
		aFree = new boolean[n];
		
		for(int i = 0; i < n; i++) {
			bFree[i] = supplies[i] != 0;
			aFree[i] = demands[i] != 0;
		}
		
		//Dual weights for B and A
		y = new int[2*n];
		
		//Remaining "unmatched" supply or demand.
		int[] deficiencyB = supplies;
		int[] deficiencyA = demands;
		capacityAB = new int[n][n];
		capacityBA = new int[n][n];
		
		for(int i = 0; i < n; i++) {
			for(int j = 0; j < n; j++) {
				capacityBA[i][j] = Math.min(deficiencyB[i], deficiencyA[j]); 
			}
		}
		
		//Preallocate maximum AP length of 2*n
		path = new ArrayList<Integer>(2*n);
		
		//Main iteration loop. Continue executing until all supply vertices have no more residual.
		while (any(bFree)) {
			iterations++;
			//Perform Dijkstra's algorithm
			//To identify the distances lv
			//in slack from the nearest free vertex of B.
			//First n vertices are type B. 
			//Second set of n vertices are type A.
			int[] lv = new int[2*n];
			
			//The minimum distance to any free vertex of A.
			int distAF = Integer.MAX_VALUE;
			
			//For each vertex, store whether this vertex has been chosen
			//As minDist before. All such vertices have already been included
			//in the Dijkstra shortest path tree, and are final.
			boolean[] finalDist = new boolean[2*n];
			Arrays.fill(lv, Integer.MAX_VALUE);
			for (int i = 0; i < n; i++) {
				if (bFree[i]) {
					lv[i] = 0;
				}
			}
			
			//Main Dijkstra loop. Each iteration adds one vertex to the shortest
			//path tree. Total time for all iterations is O(n^2).
			//Will break out early if free vertex of A is reached.
			for (int i = 0; i < 2*n; i++) {
				//Find the next vertex to add to the shortest path tree
				int minDist = Integer.MAX_VALUE;
				int minIndex = -1; //Placeholder
				for (int v = 0; v < 2*n; v++) { 
					if (lv[v] < minDist && !finalDist[v]) {
						minDist = lv[v];
						minIndex = v;
					}
				}
				finalDist[minIndex] = true;
				
				//From the for loop condition, the early breaking out upon
				//reaching a free vertex, and the fact that the total supply is 
				//less than the total demand, there should always be an augmenting
				//path found, meaning minDist < Inf. The asserts double check this if enabled.
				assert(minIndex != -1);
				assert(minDist < Integer.MAX_VALUE);
				
				if (minIndex < n) {
					//Add a vertex of type B to the tree
					//Update distances to all neighbors in A
					for (int a = 0; a < n; a++) {
						if (capacityBA[minIndex][a] > 0) {
							int aIndex = a + n;
							int newDist = lv[minIndex] + CBA[minIndex][a] + 1 - y[minIndex] - y[aIndex];
							if (newDist < lv[aIndex]) {
								lv[aIndex] = newDist; 
							}
						}
					}
				}
				else {
					//Add a vertex of type A to the tree
					//Update distances to all neighbors in B
					int a = minIndex - n;
					if (aFree[a]) {
						distAF = lv[minIndex];
						break;
					}
					for (int b = 0; b < n; b++) {
						if (capacityAB[a][b] > 0) {
							int newDist = lv[minIndex] + y[a + n] + y[b] - CAB[a][b];
							if (newDist < lv[b]) {
								lv[b] = newDist;
							}
						}
					}
				}
			}
		
			//Since there is a free vertex of B left, there is always some AP left
			assert(distAF < Integer.MAX_VALUE);

			assert(distAF > 0); //Ensures that a maximal set of vertex-disjoint shortest augmenting paths was found last phase.
			
			//Now, perform dual adjustments 
			for (int i = 0; i < lv.length; i++) {
				int delta = Math.max(distAF - lv[i], 0);
				if (i < n) { 
					//i is a vertex of B; increase dual
					 y[i] += delta;
					 assert(y[i] >= 0);
				}
				else {
					// i is a vertex of A. Decrease dual.
					y[i] -= delta;
					assert(y[i] <= 0);
				}
			}

			//Let the admissible graph be the set of 0 slack edges.
			//Now, we iteratively execute partial DFS searches from each free
			//vertex of B to find a maximal set of vertex-disjoint admissible
			//augmenting paths.
			
			//This is used by the DFS procedure to track the
		    //largest explored neighbor index of every vertex.
		    //Following our convention, 0:n-1 -> B and n:2n-1 -> A.
			//These values persist throughout all partial DFS searches this phase
			vertexVisited = new int[2*n];
			for (int i = 0; i < n; i++) {
				vertexVisited[i] = n-1;
			}
			for (int i = 0; i < n; i++) {
				vertexVisited[i + n] = -1;
			}
			
			for (int vertex = 0; vertex < n; vertex++) {
				if (bFree[vertex]) {
					//For each free vertex, repeatedly find admissible APs
					//until no more can be found.
					while (deficiencyB[vertex] > 0 && vertexVisited[vertex] < 2*n - 1) {
						ArrayList<Integer> ap = null;
						
						ap = findAP(vertex);
						
						//Comment this out to make code faster
						//Uncomment if time taken by augmentations needs to be measured.
						//long startTimeAugment = System.nanoTime();
						if (ap.size() == 0) {
							//No AP was found. Move to next free vertex.
							break;
						}
						APLengths += ap.size() - 1;
						
						//Compute the bottleneck capacity beta: the maximum amount of flow that can be pushed
						//without violating some vertex / edge capacity constraint.
						int beta = Integer.min(deficiencyB[ap.get(0)], deficiencyA[ap.get(ap.size() - 1) - n]);
						for (int j = 0; j < ap.size() - 1; j++) {
							int u = ap.get(j);
							int v = ap.get(j + 1);			
							if (u >= n) {							
								//edge is directed from A to B
								beta = Integer.min(beta, capacityAB[u - n][v]);							
							}
							else {
								//edge is directed from B to A
								beta = Integer.min(beta, capacityBA[u][v - n]);
							}
						}
						
						//Augment along path by value beta
						//First, update the edge capacities.
						for (int j = 0; j < ap.size() - 1; j++) {
							int u = ap.get(j);
							int v = ap.get(j + 1);
							if (u >= n) {
								//edge is directed from A to B
								capacityAB[u - n][v] -= beta;
								capacityBA[v][u - n] += beta;
								if (capacityAB[u - n][v] > 0) {
									//Allow edge to be reused on future augmenting
									//paths this phase:
									vertexVisited[u] = v - 1;
								}
							}
							else {
								//edge is directed from B to A
								capacityBA[u][v - n] -= beta;
								capacityAB[v - n][u] += beta;
								if (capacityBA[u][v - n] > 0) {
									//Allow edge to be reused on future augmenting
									//paths this phase:
									vertexVisited[u] = v - 1;
								}
							}
						}
						
						//Next, update the deficiencies of the endpoints of the path
						int first = ap.get(0);
						deficiencyB[first] -= beta;
						if (deficiencyB[first] == 0) {
							bFree[first] = false;
						}
						
						int last = ap.get(ap.size() - 1) - n;
						deficiencyA[last] -= beta;
						if (deficiencyA[last] == 0) {
							aFree[last] = false;
						}
						
						//Comment this out to make code faster.
						//Uncomment if time taken by augment procedure needs to be measured.
						//timeTakenAugment += System.nanoTime() - startTimeAugment;
					}
				}
			}//End of main for loop for DFS'
			
			
		}//End of while loop (phase loop)
		
		//Finally, convert capacities to a slightly different format used
		//by the rest of the program (i.e., the Mapping class).
		
		finalCapacity = new int[2*n][2*n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				finalCapacity[i + n][j] = capacityAB[i][j];
				finalCapacity[i][j + n] = capacityBA[i][j];
			}
		}
		
		timeTaken = System.currentTimeMillis() - startTime;
	}
	
	//Execute a partial DFS from the given vertex.
	//If an admissible AP is found, return it.
	//else, return an empty list.
	public ArrayList<Integer> findAP(int start) {
		numAPCalls++;
		//Path is a preallocated array of size 2*n
		path.clear(); 
		path.add(start);
		while (!path.isEmpty()) {
			int end = path.get(path.size() - 1);
			if (end >= n && aFree[end - n]) {
				//Found an AP
				return path; 
			}
			//Attempt to expand path
			boolean backtrack = true;
			int rangeStart = vertexVisited[end] + 1;
			int rangeEnd = end < n ? 2*n : n;
			for (int i = rangeStart; i < rangeEnd; i++) { 
				vertexVisited[end] = i;
				if (end < n) {
					//current vertex is type B
					int a = i - n;
					if (CBA[end][a] + 1 - y[end] - y[a + n] == 0 && capacityBA[end][a] > 0) {
						backtrack = false;
						//Add vertex to path
						path.add(i);
						break;
					}
				}
				else {
					//current vertex is type A
					int a = end - n;
					if (capacityAB[a][i] > 0 && y[a + n] + y[i] == CAB[a][i]) {
						backtrack = false;
						//Add vertex to path
						path.add(i);
						break;
					}
				}
			}
			if (backtrack) {
				//No more options to explore from this vertex. Remove
				//last vertex from path.
				path.remove(path.size() - 1);
			}
		
		}
		return path;
	}
	
	//Returns true if there is any instance of true in the input array.
	//Otherwise, returns false.
	static boolean any(boolean[] free) {
		for (int i = 0; i < free.length; i++) {
			if (free[i]) {
				return true;
			}
		}
		return false;
	}
	
	public static int[][] transpose(int[][] matrix){
		int[][] t = new int[matrix.length][matrix[0].length];
		for(int i = 0; i < matrix.length; i++) {
			for(int j = 0; j < matrix[0].length; j++) {
				t[j][i] = matrix[i][j];
			}
		}
		return t;
	}
	
	//Some methods for outputting matrix contents if desired.
	public static String arrToString(int[][] arr) {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < arr.length; i++) {
			for (int j = 0; j < arr[0].length; j++) {
				sb.append(arr[i][j]);
				sb.append(' ');
			}
			sb.append('\n');
		}
		return sb.toString();
	}
	
	public static String arrToString(int[] arr) {
		StringBuilder sb = new StringBuilder();
		for (int j = 0; j < arr.length; j++) {
			sb.append(arr[j]);
			sb.append(' ');
		}
		sb.append('\n');
		
		return sb.toString();
	}
	
	
	
	
	
	
}