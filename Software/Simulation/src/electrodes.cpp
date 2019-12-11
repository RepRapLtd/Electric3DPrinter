/*
 * electrodes.cpp
 *
 *  Created on: Dec 11, 2019
 *      Author: ensab
 */

#include "electric3Dprinter.h"

// Find the 2D boundary pattern.  This is a disc in the middle of the Z range,
// and is the same for all Z values.  This uses the 3D OnBoundary() function above, rather
// than bother with an extra 2D version.  It relies on the order of scanning the nodes
// to get the order round the boundary right.

void Electrodes::FindBoundary()
{
	// Find the first quadrant

	electrodeCount = 0;

	// Assume that half way up is immune from top and bottom boundary interference

	int z = nodes/2;

	for(int x = 0; x <= tank.xCentre; x++)
	{
		for(int y = tank.yCentre; y >= 0; y--)
		{
			if(tank.OnBoundary(x, y, z))
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
		electrodeNodes[electrodeCount][0] = 2*tank.xCentre - electrodeNodes[x][0];
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
		electrodeNodes[electrodeCount][1] = 2*tank.yCentre - electrodeNodes[x][1];
		electrodeCount++;
	}

	if(electrodeCount%4)
		cout << "Number of boundary nodes is not a multiple of 4! " << electrodeCount << endl;
}

// Create an quarter of a circle radius r anticlockwise from the X axis.
// It returns the count of the "pixels".  Horn's algorithm.

int Electrodes::QuarterCircleDDA(int r, int circle[nodes][2])
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

void Electrodes::ElectrodeRing(int r, int z)
{
	int circle[nodes][2];

	int count = QuarterCircleDDA(r, circle);

	if(count <= 1)
	{
		electrodeNodes[electrodeCount][0] = tank.xCentre;
		electrodeNodes[electrodeCount][1] = tank.yCentre;
		electrodeNodes[electrodeCount][2] = z;
		electrodeCount++;
		return;
	}


	for(int c = 0; c < count; c++)
	{
		electrodeNodes[electrodeCount][0] = tank.xCentre + circle[c][0];
		electrodeNodes[electrodeCount][1] = tank.yCentre + circle[c][1];
		electrodeNodes[electrodeCount][2] = z;

		electrodeNodes[electrodeCount + count][0] = tank.xCentre - circle[c][1];
		electrodeNodes[electrodeCount + count][1] = tank.yCentre + circle[c][0];
		electrodeNodes[electrodeCount + count][2] = z;

		electrodeNodes[electrodeCount + 2*count][0] = tank.xCentre - circle[c][0];
		electrodeNodes[electrodeCount + 2*count][1] = tank.yCentre - circle[c][1];
		electrodeNodes[electrodeCount + 2*count][2] = z;

		electrodeNodes[electrodeCount + 3*count][0] = tank.xCentre + circle[c][1];
		electrodeNodes[electrodeCount + 3*count][1] = tank.yCentre - circle[c][0];
		electrodeNodes[electrodeCount + 3*count][2] = z;

		electrodeCount++;
	}

	electrodeCount += 3*count;
}


// Construct a spherical shell of electrodes.

void Electrodes::ElectrodeSphere()
{
	int circle[nodes][2];
	int halfNumberOfRings = QuarterCircleDDA(tank.radius, circle);

	int zLow = tank.zCentre - tank.radius;
	int r, z;
	for(int h = 0; h < halfNumberOfRings; h++)
	{
		r = circle[h][1];
		z = zLow + tank.radius - circle[h][0];
		cout << "r = " << r << ", z = " << z << endl;
		ElectrodeRing(r, z);
	}
	for(int h = halfNumberOfRings - 1; h > 0 ; h--)
	{
		r = circle[h][1];
		z = tank.zCentre + circle[h][0] - 1;
		cout << "r = " << r << ", z = " << z << endl;
		ElectrodeRing(r, z);
	}
}




// Save the electrode pattern as a .asc point cloud

void Electrodes::OutputElectrodes()
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

#ifdef TESTING

// Test the circle DDA.  Not normally called.

void Electrodes::TestCircle(int r)
{
	electrodeCount = 0;
	ElectrodeRing(r, 0);
	//cout << electrodeCount << endl;
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
			tank.potential[x][y][0] = 0.0;
	}
	for(int e = 0; e < electrodeCount; e++)
	{
		//cout << e << ' ' << electrodeNodes[e][0] - xCentre << ' ' << electrodeNodes[e][1] - yCentre << endl;
		tank.potential [electrodeNodes[e][0]] [electrodeNodes[e][1]] [0] += 1.0;
	}

	for(int y = nodes; y >= 0 ; y--)
	{
		for(int x = 0; x <= nodes; x++)
		{
			if(tank.potential[x][y][0] < 0.5)
				cout << '.';
			else
				cout << round(tank.potential[x][y][0]);
		}
		cout << endl;
	}
	cout << endl;
}

// Output the pattern of electrodes as a test.
// Not normally called.

void Electrodes::TestElectrodePattern()
{
	tank.Initialise();
	ElectrodeSphere();
	OutputElectrodes();
}

#endif

