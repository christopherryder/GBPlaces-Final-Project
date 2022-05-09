/*
Author: Christopher Ryder
Email: Christopher.Ryder-2@student.manchester.ac.uk
Created on Saturday December 1 11:32:41 2018

Declare a class for implementing Christofides Heurestic for a TSP problem for a given vector of places.

Christofides heurestic is a heurestic method for solving the tsp, designed to generate solution within 3/2 the optimal solution.
[https://www2.seas.gwu.edu/~simhaweb/champalg/tsp/tsp.html]
[https://en.wikipedia.org/wiki/Christofides_algorithm]

The heuristic algorithm in its most general form is:

The general algorithm for the christofides heuristic is as follows:
	1. Find a Minimum Spanning Tree (MST).
	2. Identify all the odd valency verticies in the 'graph'
	3. Find the minimum matching between all of these odd valency verticies.
		i. this is the same as finding the eulerian cycle.
	4. Traverse the cycle, starting a vertex, V and skip previously visited verticies.

The first link, seas.gwu.edu, presents an interesting idea for further optimising the solution using K-OPT optimisation, this could be an interesting extension.
[https://en.wikipedia.org/wiki/2-opt]
[https://ocw.mit.edu/courses/sloan-school-of-management/15-053-optimization-methods-in-management-science-spring-2013/lecture-notes/MIT15_053S13_lec17.pdf]
[http://pedrohfsd.com/2017/08/09/2opt-part1.html]

This method (un)fortunately had me dive quite deep down a rabbit hole of graph theory...
*/

//header guards stop preprocessor from including files multiple times -> breaks code as redefinitions are usually bad
#ifndef __CPATHANALYSIS_HEADER__
#define __CPATHANALYSIS_HEADER__

#include <vector>

class CPathAnalysis
{
public:

	//Method takes a distance matrix as constructor, D_ij must equal D_ji and the trace must be 0.
	CPathAnalysis(std::vector<std::vector<double>> distanceMatrix);

	//finds a MST via prims algorithm.
	void findMST();

	//finds the minimumMatching of the odd nodes of the MST.
	void minimumMatching();

	//find the eulerian cycle.
	void findEulerian();

	//find the hamiltonian cycle: => skip repeated nodes in the eulerian.
	void findHamiltonCycle();

	//call all the methods.
	std::vector<int> getFinalPath();
	int getVerticies() { return mVerticies; }

	//getter functions, keep them inline as theres no need to take up space in the cpp file.
	double getPathLength() { return mPathLength; }

private:

	//define a struct to hold information about each vertex.
	struct vertex {

		//parent of the vertex.
		int parent;

		//critical value, used for selection of verticies - normally the corresponding distance matrix value.
		double cVal;

		//hold adjacency relations.
		std::vector<int> adjacent;

		bool inMST;

		//is of odd valency.
		bool isOdd = false;
		bool visited = false;
	};

	//define a graph, a collection of verticies.
	std::vector<vertex> mGraph;
	int mVerticies;

	//distance matrix holds the value of the distance from each point to each point.
	std::vector<std::vector<double>> mDistanceMatrix;

	//hold the path in terms of the index of the place, hence a vector of ints.
	std::vector<int> mPath;

	//the final path, to be returned.
	std::vector<int> mFinalPath;
	double mPathLength;

	//list off odd valency verticies.
	std::vector<int> mOddVerticies;
};
#endif