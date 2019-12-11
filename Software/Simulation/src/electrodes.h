/*
 * electrodes.h
 *
 *  Created on: Dec 11, 2019
 *      Author: ensab
 */

#ifndef ELECTRODES_H_
#define ELECTRODES_H_

// Comment this line out if you are not using the testing functions

#define TESTING

class FiniteDifference;
extern int nodes;

class Electrodes
{
friend class FiniteDifference;

public:

private:
	FiniteDifference tank;

	// List of nodes on the boundary in cyclic order and ascending z. The maximum possible
	// number, maxBoundary, is all the points on the X, Y faces.  These are the electrodes.

	const int maxElectrodes = 4*nodes*nodes;
	int electrodeCount = 0;
	int electrodeNodes[maxElectrodes][3];

	void FindBoundary();
	int QuarterCircleDDA(int r, int circle[nodes][2]);
	void ElectrodeRing(int r, int z);
	void ElectrodeSphere();
	void OutputElectrodes();

#ifdef TESTING
	void TestCircle(int r);
	void TestElectrodePattern();
#endif
};



#endif /* ELECTRODES_H_ */
