/*
 * electric3Dprinter.cpp
 *
 * Solves Laplace's equation for a potential field given source and sink voltages at points on the periphery.
 *
 * The region of solution is a cylinder, and the sources and sinks are specified as a list of potentials.
 *
 * The sources and sinks can be moved between several solutions and the cumulated charge moved through each pixel
 * calculated.
 *
 *  Created on: 26 July 2019
 *
 *      Author: Adrian Bowyer
 *              RepRap Ltd
 *              https://reprapltd.com
 *
 *      Licence: GPL
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
 *
 * To convert output to STL:
 *
 *  $ dmc -tensor myoutputdata.tns -iso 0.5 -out myobject.stl
 *
 */

#include "electric3Dprinter.h"

FiniteDifference tank();

//**********************************************************************************************

// For timing comparisons

struct timespec diff(struct timespec start, struct timespec end)
{
	struct timespec temp;

	if(end.tv_sec - start.tv_sec == 0)
	{
		temp.tv_nsec = end.tv_nsec - start.tv_nsec;
	} else
	{
		temp.tv_nsec = ((end.tv_sec - start.tv_sec)*1000000000) + end.tv_nsec - start.tv_nsec;
	}

	return temp;
}

// Find the 2D boundary pattern.  This is a disc in the middle of the Z range,
// and is the same for all Z values.  This uses the 3D OnBoundary() function above, rather
// than bother with an extra 2D version.  It relies on the order of scanning the nodes
// to get the order round the boundary right.

void FindBoundary()
{
	// Find the first quadrant

	electrodeCount = 0;

	// Assume that half way up is immune from top and bottom boundary interference

	int z = nodes/2;

	for(int x = 0; x <= xCentre; x++)
	{
		for(int y = yCentre; y >= 0; y--)
		{
			if(OnBoundary(x, y, z))
			{
				if(electrodeCount >= maxElectrodes)
				{
					cout << "Maximum boundary node count (Q0) exceeded!" << endl;
					return;
				}

				electrodeNodes[electrodeCount][0] = x;
				electrodeNodes[electrodeCount][1] = y;
				electrodeCount++;
			}
		}
	}

	// Copy that quadrant to the other three

	int bc = electrodeCount;
	for(int x = bc - 2; x >= 0; x--)
	{
		if(electrodeCount >= maxElectrodes)
		{
			cout << "Maximum boundary node count (Q1) exceeded!" << endl;
			return;
		}
		electrodeNodes[electrodeCount][0] = 2*xCentre - electrodeNodes[x][0];
		electrodeNodes[electrodeCount][1] = electrodeNodes[x][1];
		electrodeCount++;
	}

	bc = electrodeCount;
	for(int x = bc - 2; x > 0; x--)
	{
		if(electrodeCount >= maxElectrodes)
		{
			cout << "Maximum boundary node count (Q34) exceeded!" << endl;
			return;
		}
		electrodeNodes[electrodeCount][0] = electrodeNodes[x][0];
		electrodeNodes[electrodeCount][1] = 2*yCentre - electrodeNodes[x][1];
		electrodeCount++;
	}

	if(electrodeCount%4)
		cout << "Number of boundary nodes is not a multiple of 4! " << electrodeCount << endl;
}

// Create an quarter of a circle radius r anticlockwise from the X axis.
// It returns the count of the "pixels".  Horn's algorithm.

int QuarterCircleDDA(int r, int circle[nodes][2])
{
	int x = r;
	int y = 0;
	int s = -r;

	// Compute one eighth of the circle. y is both coordinate and counter.

	while(x >= y)
	{
		circle[y][0] = x;
		circle[y][1] = y; // More confusing not to...
		s += 2*y + 1;
		y++;
		if(s > 0)
		{
			s -= 2*x - 2;
			x--;
		}
	}

	int yTop = y;
	y *= 2;

	// If x == y at the eighth-circle point, don't repeat that point.

	if(circle[yTop - 1][0] == circle[yTop - 1][1])
	{
		yTop--;
		y--;
	}

	for(int c = 1; c < yTop; c++)
	{
		circle[y - c - 1][0] = circle[c][1];
		circle[y - c - 1][1] = circle[c][0];
	}

	return y - 1;
}


// Add a circular ring of electrodes of radius r at height z.

void ElectrodeRing(int r, int z)
{
	int circle[nodes][2];

	int count = QuarterCircleDDA(r, circle);

	if(count <= 1)
	{
		electrodeNodes[electrodeCount][0] = xCentre;
		electrodeNodes[electrodeCount][1] = yCentre;
		electrodeNodes[electrodeCount][2] = z;
		electrodeCount++;
		return;
	}


	for(int c = 0; c < count; c++)
	{
		electrodeNodes[electrodeCount][0] = xCentre + circle[c][0];
		electrodeNodes[electrodeCount][1] = yCentre + circle[c][1];
		electrodeNodes[electrodeCount][2] = z;

		electrodeNodes[electrodeCount + count][0] = xCentre - circle[c][1];
		electrodeNodes[electrodeCount + count][1] = yCentre + circle[c][0];
		electrodeNodes[electrodeCount + count][2] = z;

		electrodeNodes[electrodeCount + 2*count][0] = xCentre - circle[c][0];
		electrodeNodes[electrodeCount + 2*count][1] = yCentre - circle[c][1];
		electrodeNodes[electrodeCount + 2*count][2] = z;

		electrodeNodes[electrodeCount + 3*count][0] = xCentre + circle[c][1];
		electrodeNodes[electrodeCount + 3*count][1] = yCentre - circle[c][0];
		electrodeNodes[electrodeCount + 3*count][2] = z;

		electrodeCount++;
	}

	electrodeCount += 3*count;
}


// Construct a spherical shell of electrodes.

void ElectrodeSphere()
{
	int circle[nodes][2];
	int halfNumberOfRings = QuarterCircleDDA(radius, circle);

	int zLow = zCentre - radius;
	int r, z;
	for(int h = 0; h < halfNumberOfRings; h++)
	{
		r = circle[h][1];
		z = zLow + radius - circle[h][0];
		cout << "r = " << r << ", z = " << z << endl;
		ElectrodeRing(r, z);
	}
	for(int h = halfNumberOfRings - 1; h > 0 ; h--)
	{
		r = circle[h][1];
		z = zCentre + circle[h][0] - 1;
		cout << "r = " << r << ", z = " << z << endl;
		ElectrodeRing(r, z);
	}
}




// Save the electrode pattern as a .asc point cloud

void OutputElectrodes()
{
	ofstream outputFile;
	outputFile.open("electrodes.asc");
	for(int ec = 0; ec < electrodeCount; ec++)
	{
		for(int coord = 0; coord < 3; coord++)
		{
			outputFile << electrodeNodes[ec][coord] << " ";
		}
		outputFile << endl;
	}
	outputFile.close();
}


// Function to plot the boundary to test that it's properly set-up.
// Not normally called.

void TestBoundaryDisc()
{
	int z = nodes/2;
	BoundaryConditions(0.0, z, 1.0);
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
			inside[x][y][z] = true;
	}
	potential[1][1][z] = 0.5;
	for(int x = 0; x < electrodeCount; x++)
	{
		potential[electrodeNodes[x][0]][electrodeNodes[x][1]][z] += 1;
		cout << x << ": (" << electrodeNodes[x][0] << ", " << electrodeNodes[x][1] << ")" << endl;
	}
	OutputDisc("boundary.dat", potential, -1, z);
}

// Create a cylinder pattern in thresholdedChargeIntegral[][][] and write it out as
// a tensor to test tensor output.  Not normally called.

void TestCylinder(const int r, const int z0, const int z1)
{
	Initialise();
	for(int x = 0; x <= nodes; x++)
	{
		int xd = x - xCentre;
		for(int y = 0; y <= nodes; y++)
		{
			int yd = y - yCentre;
			for(int z = 0; z <= nodes; z++)
			{
				if(z < z0 || z > z1 || (xd*xd + yd*yd) > r*r)
				{
					thresholdedChargeIntegral[x][y][z] = 0.0;
				} else
				{
					thresholdedChargeIntegral[x][y][z] = 1.0;
				}
			}
		}
	}
	OutputTensor("testCylinder.tns", thresholdedChargeIntegral);
}

// Test the circle DDA.  Not normally called.

void TestCircle(int r)
{
	electrodeCount = 0;
	ElectrodeRing(r, 0);
	//cout << electrodeCount << endl;
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
			potential[x][y][0] = 0.0;
	}
	for(int e = 0; e < electrodeCount; e++)
	{
		//cout << e << ' ' << electrodeNodes[e][0] - xCentre << ' ' << electrodeNodes[e][1] - yCentre << endl;
		potential [electrodeNodes[e][0]] [electrodeNodes[e][1]] [0] += 1.0;
	}

	for(int y = nodes; y >= 0 ; y--)
	{
		for(int x = 0; x <= nodes; x++)
		{
			if(potential[x][y][0] < 0.5)
				cout << '.';
			else
				cout << round(potential[x][y][0]);
		}
		cout << endl;
	}
	cout << endl;
}

// Output the pattern of electrodes as a test.
// Not normally called.

void TestElectrodePattern()
{
	Initialise();
	ElectrodeSphere();
	OutputElectrodes();
}

// (Re)set everything to the initial state

void Reset()
{
	Initialise();
	FindBoundary();
	BoundaryConditions(0, nodes/2, 1.0);
}

// Remind the user what they can do.

void Prompt()
{
	cout << endl << "Commands:" << endl;
	cout << " s: - set sigmoid threshold and output a slice" << endl;
	cout << " r: - reset to the initial state" << endl;
	cout << " t: - create the tensor file" << endl;
	cout << " h: - print this list" << endl;
	cout << " q: - quit" << endl<< endl;
}

// Decide how to process the results.

void Control()
{
	double s;
	int z;
	string fName;

	cout << "Type h for help." << endl;
	while(1)
	{
		cout << "Command: ";
		char c;
		cin >> c;

		switch(c)
		{
		case 't':
			OutputTensor("thresholdTensor.tns", thresholdedChargeIntegral);
			break;

		case 's':
			cout << "Sigmoid value for 0.5 point: ";
			cin >> s;
			SigmoidCharge(s, sigPower);
			cout << "Z value for slice [0, " << nodes << "]: ";
			cin >> z;
			cout << "Output file name: ";
			cin >> fName;
			OutputDisc(fName.c_str(), thresholdedChargeIntegral, -1, z);
			break;

		case 'r':
			Reset();
			break;

		case 'q':
			return;

		default:
			cout << endl << "Unrecognised command - " << c << endl;
		case 'h':
			Prompt();
		}
	}
}


// Self-explanatory, I hope.

int main()
{
//	struct timespec t_start, t_end;
//	clock_gettime(CLOCK_MONOTONIC, &t_start);
//
//	//	TestBoundary();
	TestElectrodePattern();
//	for(int r = 0; r <= radius; r++)
//	{
//		cout << endl << "Radius: " << r << endl;
//		TestCircle(r);
//	}
//	//  TestCylinder(15, 10, 40);
//  //  TestCircle(24);
//
//  Reset();
//
//
//	//for(int vv = 1; vv < 21; vv++)
//	//{
//	//int vv = 3;
//		ChargeSetUp();
//		//cout << "Z for " << vv << " volts: ";
//		for(int z = 1; z < nodes; z++)
//		{
//			cout << z << ' ';
//			cout.flush();
//
//			double v = 1.0;
//
//			for(int angle = 0; angle < electrodeCount/2; angle++)
//			{
//				BoundaryConditions(angle, z, v);
//				GausSeidelIteration();
//				GradientMagnitudes();
//			}
//		}
//		cout << endl;
//		SigmoidCharge(0.0, 50);
//		string fileName = "rectangleAttemptParallel.tns";
//		//char str[20];
//		//sprintf(str,"-%d.tns",vv);
//		//fileName += str;
//		OutputTensor(fileName.c_str(), thresholdedChargeIntegral);
//	//}
//
////	Output("potential.dat", potential, -1, nodes/2);
////	Output("field.dat", field, -1, nodes/2);
////	Output("charge.dat", chargeIntegral, -1, nodes/2);
////
////	PrintChargeRange();
////
////	Control();
//
//		clock_gettime(CLOCK_MONOTONIC, &t_end);
//		cout << "Execution time: " << (double)(diff(t_start, t_end).tv_nsec / 1000000000.0) << "s" << endl;

}






