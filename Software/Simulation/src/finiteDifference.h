/*
 * finiteDifference.h
 *
 *  Created on: Dec 10, 2019
 *      Author: ensab
 */

#ifndef FINITEDIFFERENCE_H_
#define FINITEDIFFERENCE_H_

// Comment this line out if you are not using the testing functions

#define TESTING

// Comment this line out if you are not using OpenMP (https://www.openmp.org/) to parallelise and speed up the Gauss-Seidel PDE solution.

#define PARALLEL

// Multiplication is faster than division...

const double oneSixth = 1.0/6.0;

// The source and sink
// NB sources must be 2*N where N is odd.

const int sources = 2;


class FiniteDifference
{

friend class Electrodes;

public:
	FiniteDifference(int nodeNumbers);
	void SetCriteria(double rlx, double conv, int mxI);
	void PrintChargeRange();
	void OutputDisc(const char* fileName, double*** a, const int activeR, const int z);
	void OutputTensor(const char* fileName, double*** a);
	Electrodes ElectrodePattern();

private:
	double PDE(const int x, const int y, const int z);
	double GaussSeidelOnePass();
	void GausSeidelIteration();
	void GradientMagnitudes();
	void SigmoidCharge(const double sigmoidOffset, const double sigmoidMultiplier);
	void ChargeSetUp();
	bool OnBoundary(const int x, const int y, const int z);
	void Initialise();
	void Reset();
	void BoundaryConditions(const int b, const int z, const double v);

#ifdef TESTING
	void TestBoundaryDisc();
	void TestCylinder(const int r, const int z0, const int z1);
#endif

	Electrodes electrodes;

	double*** potential;
	double*** lastPotential;

	// potential[][][] is the potential field, V. lastPotential[][][] is a copy of it from the last iteration.
	// field[][][] is used to store the magnitude of the field vectors, computed at the end of one solution.
	// chargeIntegral[][][] is the accumulated charge that has flowed through each node for all solutions.
	// thresholdedChargeIntegral[][][] is a thresholded version of chargeIntegral[][][] from a sigmoid function.

	double*** field;
	double*** chargeIntegral;
	double*** thresholdedChargeIntegral;

	// True in the active region. This could be computed on the fly; but it's faster
	// to store it.  What's memory for?

	bool*** inside;

	int source[sources][3];

	// The centre and radius of the disc that is a cross section of the cylinder

	int xCentre;
	int yCentre;
	int zCentre;
	int radius;

	// Gauss-Seidel convergence criterion

	const double convergence = 0.0001;

	// Gauss-Seidel relaxation (1.0 to turn off).

	double relax = 1.0;

	// Stop Gauss-Seidel after this many iterations if there's no convergence

	int maxIterations = 3000;

};

inline Electrodes FiniteDifference::ElectrodePattern() { return electrodes; }


#endif /* FINITEDIFFERENCE_H_ */
