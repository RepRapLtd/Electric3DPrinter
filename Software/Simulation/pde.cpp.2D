/*
 * pde.cpp
 *
 * Solves Laplace's equation for a potential field given source and sink voltages at points on the periphery.
 *
 * The region of solution is a disk, and the sources and sinks are specified as a list of potentials.
 *
 * The sources and sinks can be moved between several solutions and the cumulated charge moved through each pixel
 * calculated.
 *
 *  Created on: 26 Jul 2019
 *
 *      Author: Adrian Bowyer
 *              RepRap Ltd
 *              https://reprapltd.com
 *
 *     Licence: GPL
 *
 *
 * To plot the results:
 *
 *  $ gnuplot
 *  gnuplot> set hidden3d
 *  gnuplot> splot 'potential.dat' with lines
 *
 * The output files are:
 *
 *   potential.dat - the electric potentials
 *   field.dat - the magnitude of the electric field
 *   charge.dat - the cumulated charge at each node
 *
 */

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Set true for progress reports etc.

bool debug = false;

// Gauss-Seidel convergence criterion

const double convergence = 0.001;

// Size of the grid

const int n = 50;
const int m = 50;

// The centre and radius of the disc

const int xc = 25;
const int yc = 25;
const int radius = 22;

// The source and sink
// NB sources must be 2*N where N is odd.

const int sources = 2;

int source[sources][2];

// List of nodes on the boundary in cyclic order

const int maxBoundary = 2*(n+m);
int boundaryCount = 0;
int boundaryNodes[maxBoundary][2];

// Stop after this many iterations if there's no convergence

const int maxIterations = 3000;

// potential[][] is the potential field, V. lastPotential[][] is a copy of it from the last iteration.
// field[][] is used to store the magnitude of the field vectors, computed at the end of one solution.
// chargeIntegral[][] is the accumulated charge that has flowed through each node for all solutions.
// thresholdedChargeIntegral[][] is a thresholded version of c[][] subject to, for example, a sigmoid function.

double potential[n+2][m+2], lastPotential[n+2][m+2], field[n+2][m+2], chargeIntegral[n+2][m+2], thresholdedChargeIntegral[n+2][m+2];

// Run the simulation this many times, incrementing the angle of the electrodes each time.

const int angles = 40;

// True in the active region. This could be computed on the fly; but it's faster
// to store it.  What's memory for?

bool inside[n+2][m+2];

//**********************************************************************************************

// Solve the PDE at a single node [i][j]

double PDE(int i, int j)
{
	// Do nothing outside the disc

	if(!inside[i][j])
		return potential[i][j];

	// Don't mess with the sources and sinks

	for(int k = 0; k < sources; k++)
	{
		if( (i == source[k][0]) && (j == source[k][1]) )
			return potential[i][j];
	}

	// Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

	double vxm = potential[i-1][j];
	if(!inside[i-1][j])
		vxm = potential[i+1][j];

	double vxp = potential[i+1][j];
	if(!inside[i+1][j])
		vxp = potential[i-1][j];

	double vym = potential[i][j-1];
	if(!inside[i][j-1])
		vym = potential[i][j+1];

	double vyp = potential[i][j+1];
	if(!inside[i][j+1])
		vyp = potential[i][j-1];

	// The actual PDE

	return 0.25*(vxm + vxp + vym + vyp);
}

// Do a single pass of the Gauss-Seidel iteration.

double GaussSeidelOnePass()
{
	double rms = 0.0;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{
			potential[i][j] = PDE(i, j);
			double r = (potential[i][j] - lastPotential[i][j]);
			rms = rms + r*r;
			lastPotential[i][j] = potential[i][j];
		}
	}
	rms = sqrt(rms/( (double)m*(double)n) );
	return rms;
}


// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the constant convergence.

void GausSeidelIteration()
{
	double rms = 100.0*convergence;
	int k = 0;
	while(k < maxIterations && rms > convergence)
	{
		rms = GaussSeidelOnePass();
		if(debug)
			cout << "Iteration: " << k << ", rms:" << rms << endl;
		k++;
	}
	if(k >= maxIterations)
		cout << "No convergence!, rms: " << rms << endl;
}


// Compute the magnitudes of the gradient vectors in e[][].  This is
// the electric field, E.

void GradientMagnitudes()
{
	double xd, yd;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{

			// At the edges use linear gradients; parabolas elsewhere

			if(!inside[i+1][j])
				xd = potential[i][j] - potential[i-1][j];
			else if(!inside[i-1][j])
				xd = potential[i+1][j] - potential[i][j];
			else
				xd = 0.5*(potential[i+1][j] - potential[i-1][j]);

			if(!inside[i][j+1])
				yd = potential[i][j] - potential[i][j-1];
			else if(!inside[i][j-1])
				yd = potential[i][j+1] - potential[i][j];
			else
				yd = 0.5*(potential[i][j+1] - potential[i][j-1]);

			field[i][j] = sqrt(xd*xd + yd*yd);
			chargeIntegral[i][j] += field[i][j];
		}
	}
}

// The sigmoid function that decides if a given charge integral will be solid

double Sigmoid(double a, double sigmoidOffset, double sMultiplier)
{
	return 1.0 - exp(sMultiplier*(a - sigmoidOffset))/(exp(sMultiplier*(a - sigmoidOffset)) + 1.0);
}


// Threshold the charges.  Note - this INVERTS the charges high gives 0; low gives 1.

void SigmoidCharge(double sigmoidOffset, double sMultiplier)
{
	double low = chargeIntegral[xc][yc];
	double high = low;
	for(int i=1;i<n;i++)
	{
		for(int j=1;j<m;j++)
		{
			if(inside[i][j])
			{
				if(chargeIntegral[i][j] < low)
					low = chargeIntegral[i][j];
				if(chargeIntegral[i][j] > high)
					high = chargeIntegral[i][j];
			}
		}
	}

	double scale = 1.0/(high - low);

	for(int i=1;i<n;i++)
	{
		for(int j=1;j<m;j++)
		{
			if(inside[i][j])
			{
				thresholdedChargeIntegral[i][j] = Sigmoid(scale*(chargeIntegral[i][j] - low), sigmoidOffset, sMultiplier);
			} else
				thresholdedChargeIntegral[i][j] = 0.0;

		}
	}
}


// Initialise the charges at the nodes

void ChargeSetUp()
{
	for(int i=1;i<n;i++)
	{
		for(int j=1;j<m;j++)
		{
			chargeIntegral[i][j] = 0.0;
		}
	}
}

bool OnBoundary(int i, int j)
{
	if(!inside[i][j])
		return false;

	for(int ii = -1; ii < 2; ii++)
	{
		for(int jj = -1; jj < 2; jj++)
		{
			if(!inside[i+ii][j+jj])
				return true;
		}
	}

	return false;
}

void FindBoundary()
{
	// Find the first quadrant

	boundaryCount = 0;
	for(int i = 0; i <= xc; i++)
	{
		for(int j = yc; j >= 0; j--)
		{
			if(OnBoundary(i,j))
			{
				if(boundaryCount >= maxBoundary)
				{
					cout << "Maximum boundary node count (Q0) exceeded!" << endl;
					return;
				}

				boundaryNodes[boundaryCount][0] = i;
				boundaryNodes[boundaryCount][1] = j;
				boundaryCount++;
			}
		}
	}

	// Copy that to the other three

	int bc = boundaryCount;
	for(int i = bc - 2; i >= 0; i--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q1) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = 2*xc - boundaryNodes[i][0];
		boundaryNodes[boundaryCount][1] = boundaryNodes[i][1];
		boundaryCount++;
	}

	bc = boundaryCount;
	for(int i = bc - 2; i > 0; i--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q34) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = boundaryNodes[i][0];
		boundaryNodes[boundaryCount][1] = 2*yc - boundaryNodes[i][1];
		boundaryCount++;
	}

	if(boundaryCount%4)
		cout << "Number of boundary nodes is not a multiple of 4! " << boundaryCount << endl;
}


// Set the boundary conditions and initialise one solution

void BoundaryConditions(int b)
{
	// Set up the active area and initialise the solution to 0.

	for(int i = 0; i <= n; i++)
	{
		int xd = i - xc;
		for(int j = 0; j <= m; j++)
		{
			int yd = j - yc;
			inside[i][j] = xd*xd + yd*yd < radius*radius;
			potential[i][j] = 0.0;
			lastPotential[i][j] = 0.0;
		}
	}

	FindBoundary();

	// Sources and sinks

//	source[0][0] = xc + round((double)(radius - 1)*cos(angle));
//	source[0][1] = yc + round((double)(radius - 1)*sin(angle));
//	source[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
//	source[1][1] = yc + round((double)(radius - 1)*sin(angle + M_PI));

	source[0][0] = boundaryNodes[b][0];
	source[0][1] = boundaryNodes[b][1];
	double angle = atan2(yc - source[0][1], xc - source[0][0]);
	int opposite = (b + boundaryCount/2)%boundaryCount;
	source[1][0] = boundaryNodes[opposite][0];
	source[1][1] = boundaryNodes[opposite][1];

//	potential[source[0][0]][source[0][1]] = 2.0 + sin(4.0*angle);
//	potential[source[1][0]][source[1][1]] = 2.0 + sin(4.0*angle + M_PI);

	potential[source[0][0]][source[0][1]] = 1.0;
	potential[source[1][0]][source[1][1]] = -1.0;
}


// Output for GNUplot.  If activeR is positive, just output that
// radius of the disc for close-ups of the middle.

void Output(char* name, double a[n+2][m+2], int activeR)
{
	// Find the most negative value in the mesh and use that
	// as the values outside the disc.

	double negValue = a[xc][yc];
	for(int i = 0; i <= n; i++)
		for(int j = 0; j <= m ; j++)
		{
			if(activeR <= 0)
			{
				if(inside[i][j] && a[i][j] < negValue)
					negValue = a[i][j];
			} else
			{
				int xd = i - xc;
				int yd = j - yc;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					if(a[i][j] < negValue)
						negValue = a[i][j];
				}
			}
		}

	// Stick the data in a file so that GNUPlot can
	// plot it.

	ofstream outputFile;
	outputFile.open(name);
	for(int i = 0; i <= n; i++)
	{
		for(int j = 0; j <= m ; j++)
		{
			double val = negValue;
			if(activeR <= 0)
			{
				if(inside[i][j])
					val = a[i][j];
			} else
			{
				int xd = i - xc;
				int yd = j - yc;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					val = a[i][j];
				}
			}
			outputFile << val << '\n';
		}
		outputFile << '\n';
	}
	outputFile.close();
}


// Function to plot the boundary to test that it's properly set-up.
// Not normally called.

void TestBoundary()
{
	BoundaryConditions(0.0);
	for(int i = 0; i <= n; i++)
	{
		for(int j = 0; j <= m ; j++)
			inside[i][j] = true;
	}
	potential[1][1] = 0.5;
	for(int i = 0; i < boundaryCount; i++)
	{
		potential[boundaryNodes[i][0]][boundaryNodes[i][1]] += 1;
		cout << i << ": (" << boundaryNodes[i][0] << ", " << boundaryNodes[i][1] << ")" << endl;
	}
	Output("boundary.dat", potential, -1);
}


// Self-explanatory, I hope.

int main()
{
//	TestBoundary();

	ChargeSetUp();
	BoundaryConditions(0);

	for(int i = 0; i < boundaryCount/4; i++)
	{
		BoundaryConditions(i);
		GausSeidelIteration();
		GradientMagnitudes();
	}
	for(int i = boundaryCount/2; i < (3*boundaryCount)/4; i++)
	{
		BoundaryConditions(i);
		GausSeidelIteration();
		GradientMagnitudes();
	}

	SigmoidCharge(0.1, 50);
	Output("threshold.dat", thresholdedChargeIntegral, -1);
	Output("potential.dat", potential, -1);
	Output("field.dat", field, -1);
	Output("charge.dat", chargeIntegral, -1);
}





