/*Author: Christopher Ryder
Email: Christopher.Ryder-2@student.manchester.ac.uk
Created on Saturday December 1 11:32:41 2018

This source file implements methods declared in the CPathAnalysis.h file. A summary of this class is given in that file.
*/

#include <algorithm>
#include <stack>

#include "CPathAnalysis.h"

CPathAnalysis::CPathAnalysis(std::vector<std::vector<double>> distanceMatrix)
{
	mDistanceMatrix = distanceMatrix;
	mVerticies = mDistanceMatrix[0].size();

	//set up the graph.
	mGraph = std::vector<vertex>(mVerticies);
	mPathLength = 0;
}

void CPathAnalysis::findMST()
{
	/*Prim's Algorithm of A Level Maths Decision 1 Fame.
	[https://www.geeksforgeeks.org/prims-mst-for-adjacency-list-representation-greedy-algo-6/]

	1) Create a set mstSet that keeps track of vertices already included in MST.
	2) Assign a key value to all vertices in the input graph. Initialize all key values as INFINITE. Assign key value as 0 for the first vertex so that it is picked first.
	3) While mstSet doesn’t include all vertices
		….a) Pick a vertex u which is not there in mstSet and has minimum key value.
		….b) Include u to mstSet.
		….c) Update key value of all adjacent vertices of u. To update the key values, iterate through all adjacent vertices. For every adjacent vertex v,
		if weight of edge u-v is less than the previous key value of v, update the key value as weight of u-v

	Interesting psuedo code for this can be found here:	http://www.public.asu.edu/~huanliu/AI04S/project1.htm */

	//initialise the verticies, none are in MST and all have inf cVal.
	for (int i = 0; i < mVerticies; i++)
	{
		mGraph[i].inMST = false;
		mGraph[i].cVal = std::numeric_limits<double>::max();
	}

	//the depot is to be chosen first, so is given a cVal of 0.
	mGraph[0].cVal = 0;
	//the depot has NO parent.
	mGraph[0].parent = -1;

	//consider for each vertex.
	for (int i = 0; i < mVerticies; i++)
	{
		//need the index with the minimum cVal that is NOT in the MST.
		int minIndex;

		double minValue = std::numeric_limits<double>::max();
		for (int j = 0; j < mVerticies; j++)
		{
			if (mGraph[j].inMST == false && mGraph[j].cVal < minValue)
			{
				minValue = mGraph[j].cVal;
				minIndex = j;
			}
		}

		//add it to the MST.
		mGraph[minIndex].inMST = true;

		//consider adjacent nodes
		for (int j = 0; j < mVerticies; j++)
		{
			// mDistMat[minIndex][j] is non zero only for adjacent vertices of j
			// 'inMST' flag for this vertex is false.
			// Update the cVal only if  mDistMat[minIndex][j] is smaller than cVal[j]
			if (mGraph[j].inMST == false && mDistanceMatrix[minIndex][j] < mGraph[j].cVal)
			{
				mGraph[j].parent = minIndex;
				mGraph[j].cVal = mDistanceMatrix[minIndex][j];
			}
		}
	}

	//consider adjacency relations.
	for (int i = 0; i < mVerticies; i++)
	{
		//define the index of the parent vertex.
		int pVertex = mGraph[i].parent;
		//parent of the terminal is -1, therefore we do not consider this. 
		if (pVertex != -1)
		{
			//the vertex, i , is adjacent to its parent.
			mGraph[i].adjacent.push_back(pVertex);
			//the parent is adjacent to its child, i.
			mGraph[pVertex].adjacent.push_back(i);
		}
	}
}

void CPathAnalysis::minimumMatching()
{
	/*I attempted to try and implement the hungarian algorithm or the edmonds blossom algorithm or some other minimal matching
	algorithm but they seem too indepth for me to be able to take advantage of them.

	This is a fairly average approximation of a minimal matching, but what is important is that there is a PERFECT matching.
	To truly say the solution to this TSP is 3/2 optimal we NEED a minimal PERFECT matching. :'(.
	Perhaps we could brute force it but I don't fancy that O(n!) complexity for n=100.*/

	//define a matrix of odd valency verticies.
	for (int i = 0; i < mVerticies; i++)
	{
		if (mGraph[i].adjacent.size() % 2 != 0)
		{
			mOddVerticies.push_back(i);
		}
	}

	//make diagonals max, so that they don't get mixed up in finding the minimum cost, cannot pair node to itself anyway; don't desire cycles.
	//best to do this in a temp distance matrix however, as we do not want to mess up our original.
	std::vector<std::vector<double>> tempDM = mDistanceMatrix;
	for (unsigned int i = 0; i < mDistanceMatrix.size(); i++)
	{
		//set the trace to be infinite, as trace is originally 0 we don't want to match a vertex to itself.
		tempDM[i][i] = std::numeric_limits<double>::max();
	}

	int row, col, matchings = 0;
	double minValue;

	for (unsigned int i = 0; i < mOddVerticies.size(); i++)
	{
		//no sense repeating the loop if we've made all the matchings we can make.
		if (matchings == mOddVerticies.size())
			break;

		minValue = std::numeric_limits<double>::max();
		for (unsigned int j = 0; j < mOddVerticies.size(); j++)
		{
			//consider only points with a smaller cost than the min value and who are not already paired.
			if (tempDM[mOddVerticies[i]][mOddVerticies[j]] < minValue && mGraph[mOddVerticies[i]].isOdd == false && mGraph[mOddVerticies[j]].isOdd == false)
			{
				//slightly confusing syntax, denotes the cost of going between the two odd verticies, i and j.
				minValue = tempDM[mOddVerticies[i]][mOddVerticies[j]];
				row = i;
				col = j;
			}
		}

		//to add a minimum value we need to be sure a vertex has not been previously added.
		if (mGraph[mOddVerticies[row]].isOdd == false && mGraph[mOddVerticies[col]].isOdd == false)
		{
			//tell the graph that these verticies have been included in the odd pair matching.
			mGraph[mOddVerticies[row]].isOdd = true;
			mGraph[mOddVerticies[col]].isOdd = true;

			//add them to the MST.
			mGraph[mOddVerticies[row]].adjacent.push_back(mOddVerticies[col]);
			mGraph[mOddVerticies[col]].adjacent.push_back(mOddVerticies[row]);

			//dont consider these points again.
			tempDM[mOddVerticies[row]][mOddVerticies[col]] = std::numeric_limits<double>::max();
			matchings += 2;
		}
		else {
			//dont consider these points again.
			tempDM[mOddVerticies[row]][mOddVerticies[col]] = std::numeric_limits<double>::max();
		}
	}
}

void CPathAnalysis::findEulerian()
{
	/*
	Eulerian Path is a path in graph that visits every edge exactly once. Eulerian Circuit is an Eulerian Path which starts and ends on the same vertex.
	[http://www.graph-magics.com/articles/euler.php]
	[https://www.geeksforgeeks.org/eulerian-path-and-circuit/]

	- from geeksforgeeks:
	Eulerian Cycle
	An undirected graph has Eulerian cycle if following two conditions are true.
	….a) All vertices with non-zero degree are connected. We don’t care about vertices with zero degree because they don’t belong to Eulerian Cycle or Path (we only consider all edges).
	….b) All vertices have even degree.
	-----------------------------------

	It should be trivial to find the eulerian cycle as we have ensured the graph has only even verticies.

	Fleury's algorithm:
	Start at an arbitrary vertex for even degree graphs.
	For each vertex
		CHOOSE an edge path whose DELETION would not disconned the graph
		IF this isn't possible CHOOSE remaining edge at the current vertex. MOVE to the end of the edge.

		At the end of the algorithm there are no edges left.
	*/

	//stacks are incredibly useful, work on the principal of elements are only added to the top and can only be removed from the top.
	//initialise a stack to hold the verticies, and by implication the edges, to be traversed.
	std::stack<int> vStack;

	//current vertex, next vertex.
	int cVertex = 0;
	int nextVertex = 0;

	while (mGraph[cVertex].adjacent.size() != 0 || !vStack.empty())
	{
		//if the vertex is neighbourless
		if (mGraph[cVertex].adjacent.size() == 0)
		{
			//add vertex to eulerian
			mPath.push_back(cVertex);

			//last in is top of the stack. LIFO.
			int lastVertex = vStack.top();

			//remove this vertex from top of stack.
			vStack.pop();

			//go back to the last vertex as we've encountered a 'dead end'
			cVertex = lastVertex;
		}
		else {

			//if the vertex has neighbours, push this to the top of the stack.
			vStack.push(cVertex);
			//the next vertex to traverse is the one at the END of the adjacency vector.
			int nextVertex = mGraph[cVertex].adjacent.back();
			//remove this from the adjacency vector, as we have saved its index as nextVertex.
			mGraph[cVertex].adjacent.pop_back();

			for (unsigned int i = 0; i < mGraph[nextVertex].adjacent.size(); i++)
			{
				if (mGraph[nextVertex].adjacent[i] == cVertex)
				{
					//remove vertex
					mGraph[nextVertex].adjacent.erase(mGraph[nextVertex].adjacent.begin() + i);
					break;
				}
			}
			cVertex = nextVertex;
		}
	}
	//add the vertex to the path.
	mPath.push_back(cVertex);
}

//same as eulerian but just skip verticies already visited.
void CPathAnalysis::findHamiltonCycle()
{
	//start at the depot.
	mFinalPath.push_back(0);

	//define an iterator to keep track of which vertex we're at.
	int vertexIterator = 0;

	//continue until the N-1 vertex.
	while (vertexIterator != mPath.size() - 1)
	{
		//add the next vertex, if it is NOT in the path.
		if (mGraph[mPath[vertexIterator + 1]].visited == false)
		{
			//aggregate the corresponding distance.
			mPathLength += mDistanceMatrix[mPath[vertexIterator]][mPath[vertexIterator + 1]];

			//increase the iterator.
			vertexIterator += 1;

			//update the path.
			mFinalPath.push_back(mPath[vertexIterator]);
			mGraph[mPath[vertexIterator]].visited = true;
		}
		else {
			//otherwise if its a repeated vertex, delete it from the path.
			mPath.erase(mPath.begin() + vertexIterator);
		}
	}
}

std::vector<int> CPathAnalysis::getFinalPath()
{
	//order of business for the algorithm.
	findMST();
	minimumMatching();
	findEulerian();
	findHamiltonCycle();

	return mFinalPath;
}
