/*
 * finiteDifference.cpp
 *
 *  Created on: Dec 10, 2019
 *      Author: ensab
 */

#include "electric3Dprinter.h"

// Initialise all the tensors for potential, field etc to 0, and also
// set up the active region.

void FiniteDifference::Initialise()
{
	// Set up the active area and initialise the solution to 0.

	for(int x = 0; x <= nodes; x++)
	{
		int xd = x - xCentre;
		for(int y = 0; y <= nodes; y++)
		{
			int yd = y - yCentre;
			bool in = xd*xd + yd*yd < radius*radius;
			for(int z = 0; z <= nodes; z++)
			{
				inside[x][y][z] = in;
				potential[x][y][z] = 0.0;
				lastPotential[x][y][z] = 0.0;
			}

			// Bottom and top

			inside[x][y][0] = false;
			inside[x][y][nodes] = false;
		}
	}
}

// Set the boundary conditions and initialise one solution with potentials applied at at z.
// b is the index into the boundary array that decides how far round the circle voltages will
// be applied.  v is the potential.

void FiniteDifference::BoundaryConditions(const int b, const int z, const double v)
{
	double halfDiagonal = 19.0;
	double angle = 2.0*M_PI*(double)b/(double)electrodeCount;
	if(angle > 0.5*M_PI)
		angle = M_PI - angle;
	double radius = halfDiagonal*sin(M_PI*0.25)/sin(M_PI*0.75 - angle); // Triangle sine rule
	if(z < 7 || z > 43)
		radius = 5;
	double voltage = RadiusToVoltage(radius);

	// Sources and sinks

	source[0][0] = electrodeNodes[b][0];
	source[0][1] = electrodeNodes[b][1];

	int opposite = (b + electrodeCount/2)%electrodeCount;
	source[1][0] = electrodeNodes[opposite][0];
	source[1][1] = electrodeNodes[opposite][1];

	potential[source[0][0]][source[0][1]][z] = voltage;
	potential[source[1][0]][source[1][1]][z] = -voltage;



	//double angle = atan2(yCentre - source[0][1], xCentre - source[0][0]);
	//	potential[source[0][0]][source[0][1]] = 2.0 + sin(4.0*angle);
	//	potential[source[1][0]][source[1][1]] = 2.0 + sin(4.0*angle + M_PI);
	//	source[0][0] = xc + round((double)(radius - 1)*cos(angle));
	//	source[0][1] = yc + round((double)(radius - 1)*sin(angle));
	//	source[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
	//	source[1][1] = yCentre + round((double)(radius - 1)*sin(angle + M_PI));
}


// Solve the PDE at a single node [x][y][z]

double FiniteDifference::PDE(const int x, const int y, const int z)
{
	// Do nothing outside the active region

	if(!inside[x][y][z])
		return potential[x][y][z];

	// Don't mess with the sources and sinks

	for(int l = 0; l < sources; l++)
	{
		if( (x == source[l][0]) && (y == source[l][1]) && (z == source[l][2]))
			return potential[x][y][z];
	}

	// Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

	double vxm = potential[x-1][y][z];
	if(!inside[x-1][y][z])
		vxm = potential[x+1][y][z];

	double vxp = potential[x+1][y][z];
	if(!inside[x+1][y][z])
		vxp = potential[x-1][y][z];

	double vym = potential[x][y-1][z];
	if(!inside[x][y-1][z])
		vym = potential[x][y+1][z];

	double vyp = potential[x][y+1][z];
	if(!inside[x][y+1][z])
		vyp = potential[x][y-1][z];

	double vzm = potential[x][y][z-1];
	if(!inside[x][y][z-1])
		vzm = potential[x][y][z+1];

	double vzp = potential[x][y][z+1];
	if(!inside[x][y][z+1])
		vzp = potential[x][y][z-1];

	// The actual PDE solution

	double solution = (vxm + vxp + vym + vyp + vzm + vzp)*oneSixth;

	// Apply the relaxation

	return relax*solution + (1.0 - relax)*lastPotential[x][y][z];
}

// Do a single pass of the Gauss-Seidel iteration.  Returns the root-mean
// square difference between this solution and the last, the size of which
// is the test of convergence.  This is optionally parallelised for speed using OpenMP.

#ifdef PARALLEL

double FiniteDifference::GaussSeidelOnePass()
{
	// Distribute thread results across an array
	double rms_x[nodes] = { 0.0 };
#pragma omp parallel for
	for (int x = 1; x < nodes; x++)
	{
		for (int y = 1; y < nodes; y++)
		{
			for (int z = 1; z < nodes; z++)
			{
				potential[x][y][z] = PDE(x, y, z);
				double r = (potential[x][y][z] - lastPotential[x][y][z]);
				rms_x[x] += r * r;
				lastPotential[x][y][z] = potential[x][y][z];
			}
		}
	}

	// Accumulate results
	double rms = accumulate(begin(rms_x), end(rms_x), (double)0.0);
	rms = sqrt(rms / (double)(nodes * nodes * nodes));
	return rms;
}

#else

double FiniteDifference::GaussSeidelOnePass()
{
	double rms = 0.0;
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				potential[x][y][z] = PDE(x, y, z);
				double r = (potential[x][y][z] - lastPotential[x][y][z]);
				rms += r*r;
				lastPotential[x][y][z] = potential[x][y][z];
			}
		}
	}
	rms = sqrt( rms/(double)(nodes*nodes*nodes) );
	return rms;
}

#endif



// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the constant convergence.

void FiniteDifference::GausSeidelIteration()
{
	double rms = 100.0*convergence;
	int i = 0;
	while(i < maxIterations && rms > convergence)
	{
		rms = GaussSeidelOnePass();
		if(debug)
			cout << "Iteration: " << i << ", rms:" << rms << endl;
		i++;
	}
	if(i >= maxIterations)
		cout << "No convergence!, rms: " << rms << endl;
}

// Compute the magnitudes of the gradient vectors in field[][][].  This is
// the electric field, E.  Also add the values to the accumulating charge integral.

void FiniteDifference::GradientMagnitudes()
{
	double xd, yd, zd;
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{

				// At the edges use linear gradients; parabolas elsewhere

				if(!inside[x+1][y][z])
					xd = potential[x][y][z] - potential[x-1][y][z];
				else if(!inside[x-1][y][z])
					xd = potential[x+1][y][z] - potential[x][y][z];
				else
					xd = 0.5*(potential[x+1][y][z] - potential[x-1][y][z]);

				if(!inside[x][y+1][z])
					yd = potential[x][y][z] - potential[x][y-1][z];
				else if(!inside[x][y-1][z])
					yd = potential[x][y+1][z] - potential[x][y][z];
				else
					yd = 0.5*(potential[x][y+1][z] - potential[x][y-1][z]);

				if(!inside[x][y][z+1])
					zd = potential[x][y][z] - potential[x][y][z-1];
				else if(!inside[x][y][z-1])
					zd = potential[x][y][z+1] - potential[x][y][z];
				else
					zd = 0.5*(potential[x][y][z+1] - potential[x][y][z-1]);

				field[x][y][z] = sqrt(xd*xd + yd*yd + zd*zd);
				chargeIntegral[x][y][z] += field[x][y][z];
			}
		}
	}
}

// Find and print the maximum and minimum values of the charge integral.

void FiniteDifference::PrintChargeRange()
{
	// Assume the very mid point is not outside the active region...

	double low = chargeIntegral[xCentre][yCentre][nodes/2];
	double high = low;

	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				if(inside[x][y][z])
				{
					if(chargeIntegral[x][y][z] < low)
						low = chargeIntegral[x][y][z];
					if(chargeIntegral[x][y][z] > high)
						high = chargeIntegral[x][y][z];
				}
			}
		}
	}

	cout << "Charge range: " << low << " to " << high << endl;
}

// Threshold the charges.  Note - this INVERTS the charges - high gives 0; low gives 1.
// The assumption is that we have a solid material that is made more liquid the greater
// the integral of the current flowing through it

void FiniteDifference::SigmoidCharge(const double sigmoidOffset, const double sigmoidMultiplier)
{
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				if(inside[x][y][z])
				{
					thresholdedChargeIntegral[x][y][z] = Sigmoid(chargeIntegral[x][y][z], sigmoidOffset, sigmoidMultiplier);
				} else
					thresholdedChargeIntegral[x][y][z] = 0.0;
			}
		}
	}
}


// Initialise the charges at the nodes.  Set them all 0 inside and outside the active region.

void FiniteDifference::ChargeSetUp()
{
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int z = 0; z <= nodes; z++)
			{
				chargeIntegral[x][y][z] = 0.0;
			}
		}
	}
}

// Is the node [x][y][z] on the boundary of the active region? I.e. is it inside neighbouring some outside.
// Tiny inefficiency - it checks [x][y][z] twice. 1/27th unnecessary...

bool FiniteDifference::OnBoundary(const int x, const int y, const int z)
{
	if(!inside[x][y][z])
		return false;

	for(int xDel = -1; xDel < 2; xDel++)
	{
		for(int yDel = -1; yDel < 2; yDel++)
		{
			for(int zDel = -1; zDel < 2; zDel++)
			{
				if(!inside[x+xDel][y+yDel][z+zDel])
					return true;
			}
		}
	}

	return false;
}

// Output one disc at z for gnuplot into file fileName.  If activeR is positive, just output that
// radius of the disc for close-ups of the middle.

void FiniteDifference::OutputDisc(const char* fileName, const double a[nodes+2][nodes+2][nodes+2], const int activeR, const int z)
{
	// Find the most negative value in the mesh and use that
	// as the values outside the disc.

	double negValue = a[xCentre][yCentre][z];
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
		{
			if(activeR <= 0)
			{
				if(inside[x][y][z] && a[x][y][z] < negValue)
					negValue = a[x][y][z];
			} else
			{
				int xd = x - xCentre;
				int yd = y - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					if(a[x][y][z] < negValue)
						negValue = a[x][y][z];
				}
			}
		}
	}

	// Stick the data in a file so that GNUPlot can
	// plot it.  Note messing about with blank lines that gnuplot needs.

	ofstream outputFile;
	outputFile.open(fileName);
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
		{
			double val = negValue;
			if(activeR <= 0)
			{
				if(inside[x][y][z])
					val = a[x][y][z];
			} else
			{
				int xd = x - xCentre;
				int yd = y - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					val = a[x][y][z];
				}
			}
			outputFile << val << '\n';
		}
		outputFile << '\n';
	}
	outputFile.close();
}

// Output the tensor to file fileName that is the entire mesh so that a 3D iso-surface STL file
// can be generated from it.

void FiniteDifference::OutputTensor(const char* fileName, const double a[nodes+2][nodes+2][nodes+2])
{
	// Find the maximum and minimum values.
	// Assume the central point is active...

	double minValue = a[xCentre][yCentre][nodes/2];
	double maxValue = minValue;

	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int z = 0; z <= nodes; z++)
			{
				if(inside[x][y][z])
				{
					if(a[x][y][z] < minValue)
						minValue = a[x][y][z];
					if(a[x][y][z] > maxValue)
						maxValue = a[x][y][z];
				}
			}
		}
	}

	cout << "Tensor minimum and maximum: " << minValue << ", " << maxValue << endl;

	ofstream outputFile;
	outputFile.open(fileName);

	// Need an extra blank layer so the marching cubes can see the top.

	outputFile << nodes+1 << ' ' << nodes+1 << ' ' << nodes+2 << ' ' << minValue << ' ' << maxValue;

	for(int z = 0; z <= nodes; z++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int x = 0; x <= nodes; x++)
			{
				outputFile << ' ';
				if(!inside[x][y][z])
				{
					outputFile << minValue;
				} else
				{
					outputFile << a[x][y][z];
				}
			}
		}
	}

	// Blank layer

	for(int y = 0; y <= nodes; y++)
	{
		for(int x = 0; x <= nodes; x++)
		{
			outputFile << ' ';
			outputFile << minValue;
		}
	}
	outputFile.close();
}

