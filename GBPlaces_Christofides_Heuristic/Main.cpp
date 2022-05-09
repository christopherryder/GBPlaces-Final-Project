/*
Author: Christopher Ryder
Email: Christopher.Ryder-2@student.manchester.ac.uk
Created on Wednesday Dec 5 18:01:22 2018

--- GBPlaces.csv FINAL PROJECT ---

SUMMARY:
This project aims to utilise the K++ means algorithm to determine 'good' places for K depots. It should be noted that K++ means does not always find the OPTIMUM depot
but this trade in accuracy is perhaps offset by the ability to have an arbitrary number of depots; though it is nonsensical to have more than the number of points. It should
be noted that if one were to enter k = 100, the algorithm has some problems assigning one point to every depot, as some points are given the same longitude and latitude:
e.g Westminster and London.

After generating K depots, we attempt to find a 'good' route through all of the points, visiting each only only once, using a heuristic solution. It should be noted that the value
given by this solution should be within 3/2 of the optimum solution but it is not /generally/ optimal. Further information on the heuristic can be found in the CPathAnalysis.H file.

After finding the pairings between the depot and cities they are written out to a text file. This file INCLUDES the POSITION of the depot. this is due to the fact that the output
is sometimes so large that on some computers the console window does not allow a user to scroll back so many lines. Similarly, good routes for each depot are also written out to a
text file, detailing the 'move' and the distance added for each movement.
---
Further Information:
I have constructed the K++ means algorithm to seek to minimise distance ONLY. I have decided AGAINST weighting by population as it doesn't really make sense to do when there is
no way to quantify demand. This perhaps could be done by setting up a 'delivery vector' which would describe the demand of each place. it should also be noted that the majority
of the places listed in GBPlaces.csv tend to cluster around Manchester/London therefore for K>1, there is little point weighting by population, although this could be trivially
done by having the 'updateMatchings' method attempt to minimise population*distance rather than just distance.

Ultimately, this code aims to minimise the distance travelled, hence weighting by population will not help in this endeavour. An interesting treatment of delivery vectors and other
such considerations (Truck capacity, etc) can be found in:

The Truck Dispatching Problem, G. B. Dantzig and J. H. Ramser, Management Science, Vol. 6, No. 1 (Oct., 1959), pp. 80-91.

NOTE: an example output of K++ Means depots for K=5 can be found in the attached files; this is purely for visualisation.
*/

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <random>
#include <cmath>
#include <algorithm>

//path analysis class and methods
#include "CPathAnalysis.h"

//Data structure to hold GBPlaces information, and information used for K++-Means.
struct Place {
	std::string name{};
	std::string type{};

	unsigned int population{};
	double latitude{};
	double longitude{};

	//wether assigned to a depot or not.
	bool assigned{ false };
};

//read in the file
void readPlaces(std::ifstream& rFileStream, std::vector<Place>& rPlaceVector)
{

	//store each line of the file as a string.
	std::string line;

	//split the document by line.
	while (std::getline(rFileStream, line, '\n'))
	{
		if (line[0] != '%')
		{
			//split the string into individual values, hold them here.
			std::vector<std::string> values;
			//define a string stream to operate on the string.
			std::istringstream ss(line);

			//split each line by comma into values.
			while (std::getline(ss, line, ','))
			{
				values.push_back(line);
			}

			//create place structure.
			Place temp;
			temp.name = values[0];
			temp.type = values[1];
			temp.population = std::stoi(values[2]);
			temp.latitude = std::stod(values[3]);
			temp.longitude = std::stod(values[4]);

			//append place to vector.
			rPlaceVector.push_back(temp);
		}
	}
}

//generate random Z(integer).
int randZ(int lower, int upper)
{
	//c++ 11 rand generator; one does not need to seed this with time.
	std::random_device r;
	std::default_random_engine engine(r());

	//employ boundary conditions for the random distribution.
	std::uniform_int_distribution<int> uniform(lower, upper);

	//generate the number using the random engine.
	int number = uniform(engine);
	return number;
}

//generate random R(real)
double randR(double lower, double upper)
{
	//c++ 11 rand generator; one does not need to seed this with time.
	std::random_device r;
	std::default_random_engine engine(r());

	//employ boundary conditions for the random distribution.
	std::uniform_real_distribution<double> uniform(lower, upper);

	//generate the number using the random engine.
	double number = uniform(engine);
	return number;
}

//calculate the great circle distance between two 'place' objects using the haversine formula.
double gcDist(const Place base, const Place place)
{
	//radius of earth in km.
	double r = 6.371e3;

	//convert lat long to radians.
	double bLat = base.latitude * (2 * 3.14159265 / 360);
	double bLong = base.longitude * (2 * 3.14159265 / 360);

	double pLat = place.latitude * (2 * 3.14159265 / 360);
	double pLong = place.longitude * (2 * 3.14159265 / 360);

	//define parameters of haversine formula.
	double havLat = 0.5 * (1 - std::cos(bLong - pLong));
	double havLong = 0.5 * (1 - std::cos(bLat - pLat));
	double latProduct = std::cos(bLat) * std::cos(bLong);

	//calculate gc dist
	double gcDistance = 2 * r * std::asin(std::sqrt(havLat + (latProduct * havLong)));
	return gcDistance;
}

void assignCentroids(std::vector<std::vector<Place>>& rAssignments, std::vector<Place>& rPlaces, int k)
{
	//Define a vector representing the squared distance from the starting centroid.
	std::vector<double> DXsq;
	double cumulativeDXsq = 0;

	//define a cumulative distribution function
	std::vector<double> cdf;

	for (unsigned int k = 0; k < rPlaces.size(); k++)
	{
		//calculate distance between kcentre and every place.
		DXsq.push_back(std::pow(gcDist(rAssignments[0][0], rPlaces[k]), 2));

		//define cdf being proportional to the square distance from the starting centroid.
		cumulativeDXsq += DXsq[k];
		cdf.push_back(cumulativeDXsq);
	}

	for (unsigned int k = 0; k < cdf.size(); k++)
	{
		//normalise cdf to unity.
		cdf[k] /= cumulativeDXsq;
	}

	while (rAssignments.size() != k)
	{
		//generate a uniform random real number.
		double p = randR(0, 1);

		for (unsigned int k = 0; k < cdf.size(); k++)
		{
			//search for lowest index such that p < cdf[i].
			if (p < cdf[k] && rPlaces[k].assigned == false)
			{
				//assign this place as a centroid.
				Place centroid = rPlaces[k];
				centroid.name = "DEPOT";
				centroid.type = "DEPOT";
				centroid.population = 100000;

				//add it to the assignments vector
				rAssignments.push_back({ centroid });
				centroid.assigned = true;

				//dont pick this point again.
				cdf[k] = -1;
				break;
			}
		}
	}
}

void matchToCentroid(std::vector<std::vector<Place>>& rAssignments, std::vector<Place>& rPlaces)
{
	for (unsigned int i = 0; i < rPlaces.size(); i++)
	{
		//initialise D to 'inf'.
		double D = std::numeric_limits<double>::max();
		int minIterator{};

		for (unsigned int k = 0; k < rAssignments.size(); k++)
		{
			//we can check if the depot position is equal to the place, if so dont both calculating for the other depots.
			if (rAssignments[k][0].latitude == rPlaces[i].latitude && rAssignments[k][0].longitude == rPlaces[i].longitude)
			{
				minIterator = k;
				break;
			}
			else if (gcDist(rAssignments[k][0], rPlaces[i]) < D)
			{
				//otherwise if the point is closer to another depot, assign it there.
				D = gcDist(rAssignments[k][0], rPlaces[i]);
				minIterator = k;
			}
		}
		if (rPlaces[i].assigned == false)
		{
			//make the assignment
			rAssignments[minIterator].push_back(rPlaces[i]);
			rPlaces[i].assigned = true;
		}
	}
}

void updateCentroids(std::vector<std::vector<Place>>& rAssignments)
{
	//consider each centroid.
	for (unsigned int i = 0; i < rAssignments.size(); i++)
	{
		//cumulative lat/longs
		double cumLat = 0;
		double cumLong = 0;

		//consider each element in the centroid. DO NOT include the parent.
		for (unsigned int j = 1; j < rAssignments[i].size(); j++)
		{
			cumLat += rAssignments[i][j].latitude;
			cumLong += rAssignments[i][j].longitude;
		}

		//as long as there is actually an assignment to the depot, needn't bother otherwise.
		if (rAssignments[i].size() > 1)
		{
			cumLat /= (rAssignments[i].size() - 1);
			cumLong /= (rAssignments[i].size() - 1);
		}

		//update the depots coordinates with the new average.
		rAssignments[i][0].latitude = cumLat;
		rAssignments[i][0].longitude = cumLong;
	}
}

bool updateMatchings(std::vector<std::vector<Place>>& rAssignments)
{
	//assume that the algorithm has reached convergence.
	bool hasUpdated = false;

	//for each depot.
	for (unsigned int i = 0; i < rAssignments.size(); i++)
	{
		//for each child of the depot.
		for (unsigned int j = 1; j < rAssignments[i].size(); j++)
		{
			//check if its better with another parent.
			std::vector<double> dist;
			for (unsigned int k = 0; k < rAssignments.size(); k++)
			{
				double distNew = gcDist(rAssignments[k][0], rAssignments[i][j]);
				dist.push_back(distNew);
			}

			int minDepot = std::distance(dist.begin(), std::min_element(dist.begin(), dist.end()));

			//no use changing the depot if it is already assigned to the optimum one.
			if (i != minDepot)
			{
				//if there has been a change, update this.
				hasUpdated = true;

				//make the new match
				rAssignments[minDepot].push_back(rAssignments[i][j]);

				//delete the old one.
				rAssignments[i].erase(rAssignments[i].begin() + j);
			}
		}
	}
	return hasUpdated;
}

int main()
{
	std::cout << "Welcome!\n\n";
	std::cout << "Reading places...\n\n";

	std::vector<Place> places;

	//input file stream, opens the file.
	std::ifstream readFileStream("GBPlaces.csv");

	//this is why I pass the fileStream as a reference, so I can neatly handle error handling in main().
	if (!readFileStream)
	{
		std::cout << "Could not open file GBPlaces.csv...\nPress enter to exit.";
		std::cin.ignore();
		return -1;
	}

	readPlaces(readFileStream, places);
	readFileStream.close();

	std::cout << "How many depots, k, would you like to initialise? (1 <= k <= 100)\nk = ";

	//take input for the k++ means algortihm
	int k;
	std::cin >> k;
	if (k <= 0 || k >= 100)
	{
		while (k < 1 || k > 100)
		{
			std::cout << "Demand input which statisfies (1 <= k <= 100).\nk = ";
			std::cin >> k;
		}
	}
	std::cout << "\n";
	std::cin.ignore();

	//to begin choose a centre.
	int randInt = randZ(0, places.size() - 1);
	std::cout << "The starting depot was placed in " << places[randInt].name << " [" << randInt << "].\n\n";

	//define a starting centroid.
	std::vector<Place> kCentre = { places[randInt] };
	kCentre[0].name = "DEPOT";
	kCentre[0].type = "DEPOT";
	//depot is assigned a population of 100k for the arbitrary reason that I use my gbplaces python plotter to visualise the K means - it scales point size with population.
	kCentre[0].population = 100000;

	//consider a vector of place vectors, the zeroth entry of each represents the depot, subsequent entries are assigned places.
	std::vector<std::vector<Place>> assignments;

	//introduce the first depot, the starting centroid.
	assignments.push_back(kCentre);

	//assign the other centroids until there are k centroids.
	assignCentroids(assignments, places, k);

	//run the initial matching of each place to each centroid.
	matchToCentroid(assignments, places);

	//repeat this process until convergence.
	bool updating = true;
	int i = 0;
	while (updating)
	{
		updateCentroids(assignments);
		updating = updateMatchings(assignments);
		i++;
	}

	std::cout << "K++ means found in " << i << " iterations, for k = " << k << ".\n";
	std::cout << "\nWriting the matching to 'depot-place-matchings.txt'.\n\n";

	//Write the depot-place matchings to file.
	std::ofstream fileStream;
	fileStream.open("depot-place-matchings.txt");

	if (!fileStream)
	{
		std::cout << "Could not write to file depot-place-matchings.txt...\nPress enter to exit.";
		std::cin.ignore();
		return -2;
	}

	//write depot matching to file, include the information about each place too.
	for (unsigned int i = 0; i < assignments.size(); i++)
	{
		fileStream << "DEPOT [" << i << "]:" << "\n";
		for (unsigned int j = 0; j < assignments[i].size(); j++)
		{
			fileStream << "\t [" << j << "] - " << assignments[i][j].name << "," << assignments[i][j].type << "," << assignments[i][j].population << "," << assignments[i][j].latitude << "," << assignments[i][j].longitude << "\n";
		}
	}

	fileStream.close();

	std::cout << "Generating routes for each depot...\n\n";
	std::cout << "Writing routes for each depot to 'depot-routes.txt'\n\n";

	//prepare a file stream to write to file.
	fileStream.open("depot-routes.txt");

	if (!fileStream)
	{
		std::cout << "Could not write to file depot-routes.txt...\nPress enter to exit.";
		std::cin.ignore();
		return -3;
	}

	//take these matchings and generate 'paths' for them all.
	//generate Distance matricies for each depot.
	double cumPath = 0;
	for (unsigned int k = 0; k < assignments.size(); k++)
	{
		std::vector<std::vector<double>> dM;
		for (unsigned int i = 0; i < assignments[k].size(); i++)
		{
			std::vector<double> row;
			for (unsigned int j = 0; j < assignments[k].size(); j++)
			{
				row.push_back(gcDist(assignments[k][i], assignments[k][j]));
			}
			dM.push_back(row);
		}

		//GENERATE Paths for each depot.
		CPathAnalysis cpa(dM);
		std::vector<int> finalPath = cpa.getFinalPath();

		//write to file.
		fileStream << "\nGenerating a route for depot [" << k << "]. \n\n";

		for (unsigned int i = 0; i < finalPath.size() - 1; i++)
		{
			fileStream << "[" << assignments[k][finalPath[i]].name << "] --> [" << assignments[k][finalPath[i + 1]].name << "]" << " - travelling " << dM[finalPath[i]][finalPath[i + 1]] << "km.\n";
		}

		fileStream << "\nTotal distance travelled for this route = " << cpa.getPathLength() << "km.";
		fileStream << "\n---------------\n\n";

		cumPath += cpa.getPathLength();
	}

	fileStream << "The total distance travelled to visit all points is " << cumPath << "km.\n";
	std::cout << "The total distance travelled to visit all points is " << cumPath << "km.\n";
	fileStream.close();

	std::cin.ignore();
}