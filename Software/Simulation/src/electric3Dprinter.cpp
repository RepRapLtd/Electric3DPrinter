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






