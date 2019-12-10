/*
 * electric3Dprinter.h
 *
 *  Created on: Dec 10, 2019
 *      Author: ensab
 */

#ifndef ELECTRIC3DPRINTER_H_
#define ELECTRIC3DPRINTER_H_

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <numeric>
#include <ctime>

#include "finiteDifference.h"

using namespace std;


enum tankShape
{
	cylinder,
	sphere
};


// Comment this line out if you are not using OpenMP (https://www.openmp.org/) to parallelise and speed up the Gauss-Seidel PDE solution.

#define PARALLEL

// Set true for progress reports etc.

bool debug = false;

// List of nodes on the boundary in cyclic order and ascending z. The maximum possible
// number, maxBoundary, is all the points on the X, Y faces.  These are the electrodes.

const int maxElectrodes = 4*nodes*nodes;
int electrodeCount = 0;
int electrodeNodes[maxElectrodes][3];

// The bigger this is, the closer the sigmoid function is to a step function

const double sigPower = 50.0;

// The sigmoid function that decides if a given charge integral will be solid

inline double Sigmoid(const double a, const double sigmoidOffset, const double sigmoidMultiplier)
{
	return 1.0 - exp(sigmoidMultiplier*(a - sigmoidOffset))/(exp(sigmoidMultiplier*(a - sigmoidOffset)) + 1.0);
}

// Empirical relationship between radius and voltage.  See the spreadsheet:
//
// diameter-vs-voltage.ods

inline double RadiusToVoltage(double r)
{
	return r*(r*(-0.0019*r + 0.3766) - 11.7326) + 101.023;
}

#endif /* ELECTRIC3DPRINTER_H_ */
