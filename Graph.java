import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.stream.Collectors;


public class Graph {

  // Keep a fast index to nodes in the map
  private Map<Integer, Vertex> vertexNames;

  /**
   * Construct an empty Graph with a map. The map's key is the name of a vertex
   * and the map's value is the vertex object.
   */
  public Graph() {
    vertexNames = new HashMap<>();
  }

  /**
   * Adds a vertex to the graph. Throws IllegalArgumentException if two vertices
   * with the same name are added.
   * 
   * @param v
   *          (Vertex) vertex to be added to the graph
   */
  public void addVertex(Vertex v) {
    if (vertexNames.containsKey(v.name))
      throw new IllegalArgumentException("Cannot create new vertex with existing name.");
    vertexNames.put(v.name, v);
  }

  /**
   * Gets a collection of all the vertices in the graph
   * 
   * @return (Collection<Vertex>) collection of all the vertices in the graph
   */
  public Collection<Vertex> getVertices() {
    return vertexNames.values();
  }

  /**
   * Gets the vertex object with the given name
   * 
   * @param name
   *          (String) name of the vertex object requested
   * @return (Vertex) vertex object associated with the name
   */
  public Vertex getVertex(String name) {
    return vertexNames.get(name);
  }

  /**
   * Adds a directed edge from vertex u to vertex v
   * 
   * @param nameU
   *          (String) name of vertex u
   * @param nameV
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addEdge(int nameU, int nameV, Double cost) {
    if (!vertexNames.containsKey(nameU))
      throw new IllegalArgumentException(nameU + " does not exist. Cannot create edge.");
    if (!vertexNames.containsKey(nameV))
      throw new IllegalArgumentException(nameV + " does not exist. Cannot create edge.");
    Vertex sourceVertex = vertexNames.get(nameU);
    Vertex targetVertex = vertexNames.get(nameV);
    Edge newEdge = new Edge(sourceVertex, targetVertex, cost);
    sourceVertex.addEdge(newEdge);
  }

  /**
   * Adds an undirected edge between vertex u and vertex v by adding a directed
   * edge from u to v, then a directed edge from v to u
   * 
   * @param name
   *          (String) name of vertex u
   * @param name2
   *          (String) name of vertex v
   * @param cost
   *          (double) cost of the edge between vertex u and v
   */
  public void addUndirectedEdge(int name, int name2, double cost) {
    addEdge(name, name2, cost);
    addEdge(name2, name, cost);
  }


  /**
   * Computes the euclidean distance between two points as described by their
   * coordinates
   * 
   * @param ux
   *          (double) x coordinate of point u
   * @param uy
   *          (double) y coordinate of point u
   * @param vx
   *          (double) x coordinate of point v
   * @param vy
   *          (double) y coordinate of point v
   * @return (double) distance between the two points
   */
  public double computeEuclideanDistance(double ux, double uy, double vx, double vy) {
    return Math.sqrt(Math.pow(ux - vx, 2) + Math.pow(uy - vy, 2));
  }

  /**
   * Computes euclidean distance between two vertices as described by their
   * coordinates
   * 
   * @param u
   *          (Vertex) vertex u
   * @param v
   *          (Vertex) vertex v
   * @return (double) distance between two vertices
   */
  public double computeEuclideanDistance(Vertex u, Vertex v) {
    return computeEuclideanDistance(u.x, u.y, v.x, v.y);
  }

  /**
   * Calculates the euclidean distance for all edges in the map using the
   * computeEuclideanCost method.
   */
  public void computeAllEuclideanDistances() {
    for (Vertex u : getVertices()){
      for (Edge uv : u.adjacentEdges) {
        Vertex v = uv.target;
        uv.distance = computeEuclideanDistance(u.x, u.y, v.x, v.y);
      }
  	}
  }



  // STUDENT CODE STARTS HERE

  public void generateRandomVertices(int n) {
    vertexNames = new HashMap<>(); // reset the vertex hashmap
   
    //creates points using random generator
    for(int i = 0; i < n ; i++){
    	Random generator = new Random();
    	Vertex x = new Vertex(i,generator.nextInt(100),generator.nextInt(100));
    	vertexNames.put(i,x);
    }
    //calcualtes ditance for all points
    for(int x : vertexNames.keySet()){
    	for(int y : vertexNames.keySet()){
    			Edge add = new Edge(vertexNames.get(x),vertexNames.get(y),100);
        		vertexNames.get(x).adjacentEdges.add(add);
    	}
    }
    computeAllEuclideanDistances(); // compute distances
  }

  //calculates nearestNeighbor, this may or may not produce the best result
  public List<Edge> nearestNeighborTsp() {
	  //takes the shortest distance
	  double closest = Integer.MAX_VALUE;
	  List<Edge> bestPath = new LinkedList<Edge>(); 
	  
	  //goes through all start points
	  for( int i = 0; i < vertexNames.size(); i++){
		  List<Edge> nearest = new LinkedList<Edge>();
		  double distance = 0;
		  
		  //initialize all points
		  for( int x : vertexNames.keySet()){
				vertexNames.get(x).known = false;	
		  }
		  //goes through all points
		  Vertex v = vertexNames.get(i);
		  v.known = true;
		  int known  = 1;
		  
		  while(known < vertexNames.size()){
			  double minimum = Integer.MAX_VALUE;
			  Edge min = v.adjacentEdges.get(0);
			  for(Edge edge : v.adjacentEdges){
				  if(edge.target.known == false && edge.distance < minimum){
					  minimum = edge.distance;
					  min = edge;
				  }
			  }
			  distance = distance + minimum;
			  v = min.target;
			  v.known = true;
			  nearest.add(min);
			  known++;
		  }
		  
		  //end case
		  for(Edge edge : v.adjacentEdges){
			  if(vertexNames.get(i).equals(edge.target)){
				  nearest.add(edge);
				  distance = distance + edge.distance;
			  }
		  }
		  //finding shortest path
		  if(closest > distance){
			  closest = distance;
			  bestPath = nearest;
		  }
	  }
	  
	  return bestPath;
  }

  public List<Edge> bruteForceTsp() {
	  LinkedList<Integer> permutation = new LinkedList<Integer>();
	  for(int i : vertexNames.keySet()){
		  permutation.add(i);
	  }
	  
	  //takes the shortest distance
	  double closest = Integer.MAX_VALUE;
	  List<Edge> bestPath = new LinkedList<Edge>(); 
	  
	  //using permute class to calculate all permutations
	  permute here = new permute(permutation);
	  LinkedList<LinkedList<Integer>> Permutation = here.list();
	  here.print();
	  //here.list()  return a list of all permutation, LinkedList<LinkedList<Integer>>
	  //using a list of edges the program goes through each of the edges to find the correct one and compiles them
	  for(LinkedList<Integer> x : Permutation){
		  double distance = 0;
		  List<Edge> path = new LinkedList<Edge>();
		  Vertex v = vertexNames.get(x.get(0));
		  for( int i = 1; i < x.size(); i++){
			  for(Edge edge : v.adjacentEdges){
				  if(edge.target.equals(vertexNames.get(x.get(i)))){
					  distance = distance + edge.distance;
					  path.add(edge);
					  v = vertexNames.get(x.get(i));
					  break;
				  }
			  }
		  }
		  
		  //last edge connecting to first edge
		  for(Edge edge : v.adjacentEdges){
			  if(vertexNames.get(x.get(0)).equals(edge.target)){
				  distance =  distance + edge.distance;
				  path.add(edge);
			  }
		  }
		  
		  //finding shortest path
		  if(closest > distance){
			  closest = distance;
			  bestPath = path;
		  }
	  }
	  return bestPath;
  }
  
  //permutation class to calculate all possible permutation
  private static class permute{
		LinkedList<Integer> path;
		LinkedList<LinkedList<Integer>> permutation;
		boolean calc = false; 
		
		permute(LinkedList<Integer> combination){
			path = combination;
			permutation = new LinkedList<LinkedList<Integer>>();
		}
		
		public LinkedList<LinkedList<Integer>> list(){
			calc = true;
			return list(path,0);
		}
		
		private LinkedList<LinkedList<Integer>> list(LinkedList<Integer> path, int location){
			//return case, final case
			if(location == path.size()){
				return permutation;
			}
			
			else{
				//base case
				if(location == 0){
					LinkedList<Integer> base = new LinkedList<Integer>();
					base.add(path.get(location));
					permutation.add(base);
					location++;
					return list(path,location);
				}
				//recurrance 
				else{
					LinkedList<LinkedList<Integer>> temp = new LinkedList<LinkedList<Integer>>();
					for(LinkedList<Integer> y: permutation){
						int size = y.size();
						for(int i = 0; i <= size; i++){
							LinkedList<Integer> x = new LinkedList<Integer>();
							x.addAll(y);
							x.add(i,path.get(location));
							temp.add(x);
						}
					}
					permutation = temp;
					location++;
					return list(path,location);
				}
			}
		}
		
		//printing permutation for checking, must be called after list if used
		public void print(){
			if (calc == true){
				System.out.println("\nHere are the permutations: ");
				for(LinkedList<Integer> x: permutation){
					System.out.println(x.toString());
				}
			}else{
				System.out.println("No permuatation calculation have been done yet.");
			}
		}
		
	}

  


  // STUDENT CODE ENDS HERE



  /**
   * Prints out the adjacency list of the graph for debugging
   */
  public void printAdjacencyList() {
    for (int u : vertexNames.keySet()) {
      StringBuilder sb = new StringBuilder();
      sb.append(u + " ");
      sb.append(vertexNames.get(u).known);
      sb.append(" -> [ ");
      for (Edge e : vertexNames.get(u).adjacentEdges) {
        sb.append(e.target.name);
        sb.append("(");
        sb.append(e.distance);
        sb.append(") ");
      }
      sb.append("]");
      System.out.println(sb.toString());
    }
  }
}
