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

#include "electrodes.h"
#include "finiteDifference.h"


using namespace std;


enum tankShape
{
	cylinder,
	sphere
};

// Set true for progress reports etc.

bool debug = false;

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
