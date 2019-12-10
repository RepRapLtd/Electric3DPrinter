/*
 * finiteDifference.h
 *
 *  Created on: Dec 10, 2019
 *      Author: ensab
 */

#ifndef FINITEDIFFERENCE_H_
#define FINITEDIFFERENCE_H_

// Multiplication is faster than division...

const double oneSixth = 1.0/6.0;

// Size of the grid

const int nodes = 50;

// The source and sink
// NB sources must be 2*N where N is odd.

const int sources = 2;


class FiniteDifference
{
public:
	FiniteDifference();
	void PrintChargeRange();
	void OutputDisc(const char* fileName, const double a[nodes+2][nodes+2][nodes+2], const int activeR, const int z);
	void OutputTensor(const char* fileName, const double a[nodes+2][nodes+2][nodes+2]);

private:
	double PDE(const int x, const int y, const int z);
	double GaussSeidelOnePass();
	void GausSeidelIteration();
	void GradientMagnitudes();
	void SigmoidCharge(const double sigmoidOffset, const double sigmoidMultiplier);
	void ChargeSetUp();
	bool OnBoundary(const int x, const int y, const int z);
	void Initialise();
	void BoundaryConditions(const int b, const int z, const double v);

	double potential[nodes+2][nodes+2][nodes+2], lastPotential[nodes+2][nodes+2][nodes+2];

	// potential[][][] is the potential field, V. lastPotential[][][] is a copy of it from the last iteration.
	// field[][][] is used to store the magnitude of the field vectors, computed at the end of one solution.
	// chargeIntegral[][][] is the accumulated charge that has flowed through each node for all solutions.
	// thresholdedChargeIntegral[][][] is a thresholded version of chargeIntegral[][][] from a sigmoid function.

	double field[nodes+2][nodes+2][nodes+2],
	       chargeIntegral[nodes+2][nodes+2][nodes+2], thresholdedChargeIntegral[nodes+2][nodes+2][nodes+2];

	// True in the active region. This could be computed on the fly; but it's faster
	// to store it.  What's memory for?

	bool inside[nodes+2][nodes+2][nodes+2];

	int source[sources][3];

	// The centre and radius of the disc that is a cross section of the cylinder

	const int xCentre = 25;
	const int yCentre = 25;
	const int zCentre = 25;
	const int radius = 22;

	// Gauss-Seidel convergence criterion

	const double convergence = 0.0001;

	// Gauss-Seidel relaxation (1.0 to turn off).

	double relax = 1.0;


	// Stop Gauss-Seidel after this many iterations if there's no convergence

	const int maxIterations = 3000;

};


#endif /* FINITEDIFFERENCE_H_ */
