# E3Dprinter.py
#
# Solves Laplace's equation for a potential field given source and sink voltages at points on the periphery.
#
# The region of solution is a cylinder, and the sources and sinks are specified as a list of potentials.
#
# The sources and sinks can be moved between several solutions and the cumulated charge moved through each pixel
# calculated.
#
#  Created on: 19 April 2021
#  This is a Python version of the previous C++ program electric3Dprinter.cpp
#
#      Author: Adrian Bowyer
#              RepRap Ltd
#              https://reprapltd.com
#
#      Licence: GPL
#
#
# To plot the results:
#
#  $ gnuplot
#  gnuplot> set hidden3d
#  gnuplot> splot 'potential.dat' with lines
#
# The output files are:
#
#   potential.dat - the electric potentials
#   field.dat - the magnitude of the electric field
#   charge.dat - the cumulated charge at each node
#
#
# To convert output to STL:
#
#  $ dmc -tensor myoutputdata.tns -iso 0.5 -out myobject.stl
#
# **********************************************************************************************

import math as maths
import time

# Set true for progress reports etc.

debug = False

# Gauss-Seidel convergence criterion

convergence = 0.0001

# Gauss-Seidel relaxation (1.0 to turn off).

relax = 1.0

# Size of the grid

nodes = 50

# The centre and radius of the disc that is a cross section of the cylinder

xCentre = 25
yCentre = 25
zCentre = 25
radius = 22

# The source and sink
# NB sources must be 2*N where N is odd.

sources = 2

source = [[0.0 for i in range(3)] for j in range(sources)]

# List of nodes on the 2D boundary in cyclic order. These are the same for all Z.

maxBoundary = 4 * nodes
boundaryCount = 0
boundaryNodes = [[0.0 for i in range(2)] for j in range(maxBoundary)]

# Stop Gauss-Seidel after this many iterations if there's no convergence

maxIterations = 3000

# potential[][][] is the potential field, V. lastPotential[][][] is a copy of it from the last iteration.
# field[][][] is used to store the magnitude of the field vectors, computed at the end of one solution.
# chargeIntegral[][][] is the accumulated charge that has flowed through each node for all solutions.
# thresholdedChargeIntegral[][][] is a thresholded version of chargeIntegral[][][] from a sigmoid function.

potential = [[[0.0 for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]
lastPotential = [[[0.0 for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]
field = [[[0.0 for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]
chargeIntegral = [[[0.0 for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]
thresholdedChargeIntegral = [[[0.0 for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]

# True in the active region. This could be computed on the fly; but it's faster
# to store it.  What's memory for?

inside = [[[True for i in range(nodes+2)] for j in range(nodes+2)] for k in range(nodes+2)]

# Multiplication is faster than division...

oneSixth = 1.0 / 6.0

# The bigger this is, the closer the sigmoid function is to a step function

sigPower = 50.0


# **********************************************************************************************

def TimeDifference(start, end):
    #	if end.tv_sec - start.tv_sec == 0:
    #		temp.tv_nsec = end.tv_nsec - start.tv_nsec
    #    else:
    #		temp.tv_nsec = ((end.tv_sec - start.tv_sec)*1000000000) + end.tv_nsec - start.tv_nsec
    return end - start


def PDE(x, y, z):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Do nothing outside the active region

    if not inside[x][y][z]:
        return potential[x][y][z]

    # Don't mess with the sources and sinks

    for l in range(sources):
        if (x == source[l][0]) and (y == source[l][1]) and (z == source[l][2]):
            return potential[x][y][z]

    # Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

    vxm = potential[x - 1][y][z]
    if not inside[x - 1][y][z]:
        vxm = potential[x + 1][y][z]

    vxp = potential[x + 1][y][z]
    if not inside[x + 1][y][z]:
        vxp = potential[x - 1][y][z]

    vym = potential[x][y - 1][z]
    if not inside[x][y - 1][z]:
        vym = potential[x][y + 1][z]

    vyp = potential[x][y + 1][z]
    if not inside[x][y + 1][z]:
        vyp = potential[x][y - 1][z]

    vzm = potential[x][y][z - 1]
    if not inside[x][y][z - 1]:
        vzm = potential[x][y][z + 1]

    vzp = potential[x][y][z + 1]
    if not inside[x][y][z + 1]:
        vzp = potential[x][y][z - 1]

    # The actual PDE solution

    solution = (vxm + vxp + vym + vyp + vzm + vzp) * oneSixth

    # Apply the relaxation

    return relax * solution + (1.0 - relax) * lastPotential[x][y][z]


# Do a single pass of the Gauss-Seidel iteration.  Returns the root-mean
# square difference between this solution and the last, the size of which
# is the test of convergence.

def GaussSeidelOnePass():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower
    rms = 0.0
    for x in range(1, nodes):
        for y in range(1, nodes):
            for z in range(1, nodes):
                potential[x][y][z] = PDE(x, y, z)
                r = (potential[x][y][z] - lastPotential[x][y][z])
                rms += r * r
    rms = maths.sqrt(rms / (nodes * nodes * nodes))
    return rms


# Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
# of the differences between the last pass and the current one is less than the value
# of the constant convergence.

def GausSeidelIteration():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower
    rms = 100.0 * convergence;
    i = 0;
    while i < maxIterations and rms > convergence:
        rms = GaussSeidelOnePass()
        if debug:
            print("Iteration: ", i, ", rms:", rms)
        i += 1
    if i >= maxIterations:
        print("No convergence!, rms: ", rms)


# Compute the magnitudes of the gradient vectors in field[][][].  This is
# the electric field, E.  Also add the values to the accumulating charge integral.

def GradientMagnitudes():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower
    for x in range(1, nodes):
        for y in range(1, nodes):
            for z in range(1, nodes):

                # At the edges use linear gradients; parabolas elsewhere

                if not inside[x + 1][y][z]:
                    xd = potential[x][y][z] - potential[x - 1][y][z]
                elif not inside[x - 1][y][z]:
                    xd = potential[x + 1][y][z] - potential[x][y][z]
                else:
                    xd = 0.5 * (potential[x + 1][y][z] - potential[x - 1][y][z])

                if not inside[x][y + 1][z]:
                    yd = potential[x][y][z] - potential[x][y - 1][z]
                elif not inside[x][y - 1][z]:
                    yd = potential[x][y + 1][z] - potential[x][y][z]
                else:
                    yd = 0.5 * (potential[x][y + 1][z] - potential[x][y - 1][z])

                if not inside[x][y][z + 1]:
                    zd = potential[x][y][z] - potential[x][y][z - 1]
                elif not inside[x][y][z - 1]:
                    zd = potential[x][y][z + 1] - potential[x][y][z]
                else:
                    zd = 0.5 * (potential[x][y][z + 1] - potential[x][y][z - 1]);

                field[x][y][z] = maths.sqrt(xd * xd + yd * yd + zd * zd)
                chargeIntegral[x][y][z] += field[x][y][z]


# Find and print the maximum and minimum values of the charge integral.

def PrintChargeRange():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Assume the very mid point is not outside the active region...

    low = chargeIntegral[xCentre][yCentre][zCentre]
    high = low


    for x in range(1, nodes):
        for y in range(1, nodes):
            for z in range(1, nodes):
                if inside[x][y][z]:
                    if chargeIntegral[x][y][z] < low:
                        low = chargeIntegral[x][y][z]
                    if chargeIntegral[x][y][z] > high:
                        high = chargeIntegral[x][y][z]

    print("Charge range: ", low, " to ", high)


# The sigmoid function that decides if a given charge integral will be solid

def Sigmoid(a, sigmoidOffset, sigmoidMultiplier):
    return 1.0 - maths.exp(sigmoidMultiplier * (a - sigmoidOffset)) / (
            maths.exp(sigmoidMultiplier * (a - sigmoidOffset)) + 1.0)


# Empirical relationship between radius and voltage.  See the spreadsheet:
#
# diameter-vs-voltage.ods

def RadiusToVoltage(r):
    return -0.0019 * r * r * r + 0.3766 * r * r - 11.7326 * r + 101.023


# Threshold the charges.  Note - this INVERTS the charges - high gives 0; low gives 1.
# The assumption is that we have a solid material that is made more liquid the greater
# the integral of the current flowing through it.

def SigmoidCharge(sigmoidOffset, sigmoidMultiplier):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower
    for x in range(1, nodes):
        for y in range(1, nodes):
            for z in range(1, nodes):
                if inside[x][y][z]:
                    thresholdedChargeIntegral[x][y][z] = Sigmoid(chargeIntegral[x][y][z], sigmoidOffset,
                                                                 sigmoidMultiplier)
                else:
                    thresholdedChargeIntegral[x][y][z] = 0.0


# Initialise the charges at the nodes.  Set them all 0 inside and outside the active region.

def ChargeSetUp():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    for x in range(nodes):
        for y in range(nodes):
            for z in range(nodes):
                chargeIntegral[x][y][z] = 0.0


# Is the node [x][y][z] on the boundary of the active region? I.e. is it inside neighbouring some outside.
# Tiny inefficiency - it checks [x][y][z] twice. 1/27th unnecessary...

def OnBoundary(x, y, z):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    if not inside[x][y][z]:
        return False

    for xDel in range(-1, 2):
        for yDel in range(-1, 2):
            for zDel in range(-1, 2):
                if not inside[x + xDel][y + yDel][z + zDel]:
                    return True
    return False


# Find the 2D boundary pattern.  This is a disc in the middle of the Z range,
# and is the same for all Z values.  This uses the 3D OnBoundary() function above, rather
# than bother with an extra 2D version.  It relies on the order of scanning the nodes
# to get the order round the boundary right.

def FindBoundary():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Find the first quadrant

    boundaryCount = 0

    # Assume that half way up is immune from top and bottom boundary interference

    z = zCentre

    for x in range(xCentre + 1):
        for y in range(yCentre, -1, -1):
            if OnBoundary(x, y, z):
                if boundaryCount >= maxBoundary:
                    print("Maximum boundary node count (Q0) exceeded!")
                    return
                boundaryNodes[boundaryCount][0] = x;
                boundaryNodes[boundaryCount][1] = y;
                boundaryCount += 1

    # Copy that quadrant to the other three

    bc = boundaryCount;
    for x in range(bc - 2, -1, -1):
        if boundaryCount >= maxBoundary:
            print("Maximum boundary node count (Q1) exceeded!")
            return
        boundaryNodes[boundaryCount][0] = 2 * xCentre - boundaryNodes[x][0]
        boundaryNodes[boundaryCount][1] = boundaryNodes[x][1]
        boundaryCount + 1

    bc = boundaryCount;
    for x in range(bc - 2, 0, -1):
        if boundaryCount >= maxBoundary:
            print("Maximum boundary node count (Q34) exceeded!")
            return
        boundaryNodes[boundaryCount][0] = boundaryNodes[x][0]
        boundaryNodes[boundaryCount][1] = 2 * yCentre - boundaryNodes[x][1]
        boundaryCount += 1

    if boundaryCount % 4 != 0:
        print("Number of boundary nodes is not a multiple of 4! ", boundaryCount)

# Initialise all the tensors for potential, field etc to 0, and also
# set up the active region.

def Initialise():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Set up the active area and initialise the solution to 0.

    for x in range(nodes + 1):
        xd = x - xCentre
        for y in range(nodes + 1):
            yd = y - yCentre
            within = xd * xd + yd * yd < radius * radius
            for z in range(nodes + 1):
                inside[x][y][z] = within
                potential[x][y][z] = 0.0
                lastPotential[x][y][z] = 0.0

            # Bottom and top

            inside[x][y][0] = False
            inside[x][y][nodes] = False


# Set the boundary conditions and initialise one solution with potentials applied at at z.
# b is the index into the boundary array that decides how far round the circle voltages will
# be applied.  v is the potential.

def BoundaryConditions(b, z, v):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    Initialise()
    FindBoundary()

    halfDiagonal = 19.0
    angle = 2.0 * maths.pi * b / boundaryCount
    if angle > 0.5 * maths.pi:
        angle = maths.pi - angle;
    radius = halfDiagonal * maths.sin(maths.pi * 0.25) / maths.sin(maths.pi * 0.75 - angle)  # Triangle sine rule
    if z < 7 or z > 43:
        radius = 5;
    voltage = RadiusToVoltage(radius)

    # Sources and sinks

    source[0][0] = boundaryNodes[b][0]
    source[0][1] = boundaryNodes[b][1]

    opposite = (b + (boundaryCount//2)) % boundaryCount
    opposite = int(opposite+0.1)
    source[1][0] = boundaryNodes[opposite][0];
    source[1][1] = boundaryNodes[opposite][1];

    potential[source[0][0]][source[0][1]][z] = voltage;
    potential[source[1][0]][source[1][1]][z] = -voltage;


# double angle = atan2(yCentre - source[0][1], xCentre - source[0][0]);
#	potential[source[0][0]][source[0][1]] = 2.0 + sin(4.0*angle);
#	potential[source[1][0]][source[1][1]] = 2.0 + sin(4.0*angle + M_PI);
#	source[0][0] = xc + round((double)(radius - 1)*cos(angle));
#	source[0][1] = yc + round((double)(radius - 1)*sin(angle));
#	source[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
#	source[1][1] = yCentre + round((double)(radius - 1)*sin(angle + M_PI));

# Output one disc at z for gnuplot into file fileName.  If activeR is positive, just output that
# radius of the disc for close-ups of the middle.

def Output(fileName, a, activeR, z):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Find the most negative value in the mesh and use that
    # as the values outside the disc.

    negValue = a[xCentre][yCentre][z]
    for x in range(nodes + 1):
        for y in range(nodes + 1):
            if activeR <= 0:
                if inside[x][y][z] and a[x][y][z] < negValue:
                    negValue = a[x][y][z]
    else:
        xd = x - xCentre
        yd = y - yCentre
        if xd * xd + yd * yd < activeR * activeR:
            if a[x][y][z] < negValue:
                negValue = a[x][y][z]

    # Stick the data in a file so that GNUPlot can
    # plot it.  Note messing about with blank lines that gnuplot needs.

    outputFile = open(fileName, "w+")
    for x in range(nodes + 1):
        for y in range(nodes + 1):
            val = negValue
            if activeR <= 0:
                if inside[x][y][z]:
                    val = a[x][y][z]
    else:
        xd = x - xCentre
        yd = y - yCentre
        if xd * xd + yd * yd < activeR * activeR:
            val = a[x][y][z]
    outputFile.write(val, '\n')
    outputFile.write('\n')
    outputFile.close()


# Output the tensor to file fileName that is the entire mesh so that a 3D iso-surface STL file
# can be generated from it.

def OutputTensor(fileName, a):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    # Find the maximum and minimum values.
    # Assume the central point is active...

    minValue = a[xCentre][yCentre][zCentre]
    maxValue = minValue;

    for x in range(nodes + 1):
        for y in range(nodes + 1):
            for z in range(nodes + 1):
                if inside[x][y][z]:
                    if a[x][y][z] < minValue:
                        minValue = a[x][y][z]
                    if a[x][y][z] > maxValue:
                        maxValue = a[x][y][z]

    print("Tensor minimum and maximum: ", minValue, ", ", maxValue)

    outputFile = open(fileName, "w+")

    # Need an extra blank layer so the marching cubes can see the top.

    outputFile.write(nodes + 1, ' ', nodes + 1, ' ', nodes + 2, ' ', minValue, ' ', maxValue)

    for z in range(nodes + 1):
        for y in range(nodes + 1):
            for x in range(nodes + 1):
                outputFile.write(' ')
                if not inside[x][y][z]:
                    outputFile.write(minValue)
        else:
            outputFile.write(a[x][y][z])

    # Blank layer

    for y in range(nodes + 1):
        for x in range(nodes + 1):
            outputFile.write(' ')
            outputFile.write(minValue)
    outputFile.close();


# Function to plot the boundary to test that it's properly set-up.
# Not normally called.

def TestBoundary():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    z = nodes//2
    BoundaryConditions(0.0, z, 1.0)
    for x in range(nodes + 1):
        for y in range(nodes + 1):
            inside[x][y][z] = True

    potential[1][1][z] = 0.5;
    for x in range(boundaryCount):
        potential[boundaryNodes[x][0]][boundaryNodes[x][1]][z] += 1
        print(x, ": (", boundaryNodes[x][0], ", ", boundaryNodes[x][1], ")")
    Output("boundary.dat", potential, -1, z)


# Create a cylinder pattern in thresholdedChargeIntegral[][][] and write it out as
# a tensor to test tensor output.  Not normally called.

def TestCylinder(r, z0, z1):
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    Initialise();
    for x in range(nodes + 1):
        xd = x - xCentre
        for y in range(nodes + 1):
            yd = y - yCentre
            for z in range(nodes + 1):
                if z < z0 or z > z1 or (xd * xd + yd * yd) > r * r:
                    thresholdedChargeIntegral[x][y][z] = 0.0
        else:
            thresholdedChargeIntegral[x][y][z] = 1.0
    OutputTensor("testCylinder.tns", thresholdedChargeIntegral)


# Remind the user what they can do.

def Prompt():
    print("\nCommands:")
    print(" s: - set sigmoid threshold and output a slice")
    print(" t: - create the tensor file")
    print(" h: - print this list")
    print(" q: - quit")


# Decide how to process the results.

def Control():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    print("Type h for help.")
    while (True):
        c = input("Command: ")
        c = c[0]

        if c == 't':
            OutputTensor("thresholdTensor.tns", thresholdedChargeIntegral)
        elif c == 's':
            s = input("Sigmoid value for 0.5 point: ")
            SigmoidCharge(s, sigPower)
            Output("threshold.dat", thresholdedChargeIntegral, -1, nodes//2)
        elif c == 'q':
            return
        elif c == 'h':
            Prompt()
        else:
            print("\nUnrecognised command - ", c)
            Prompt();


def Run():
    global debug,convergence,relax,nodes,xCentre,yCentre,zCentre,radius,sources,source,maxBoundary,boundaryCount,boundaryNodes,\
        maxIterations,potential,lastPotential,field,chargeIntegral,thresholdedChargeIntegral,inside,oneSixth,sigPower

    startTime = time.time()

    # TestBoundary();
    #  TestCylinder(15, 10, 40);

    BoundaryConditions(0, nodes//2, 1.0);

    # for(int vv = 1; vv < 21; vv++)
    # {
    # nt vv = 3;
    ChargeSetUp();
    # cout << "Z for " << vv << " volts: ";
    for z in range(nodes):
        print(z, ' ', end='', flush=True)

        v = 1.0;

        for angle in range(boundaryCount//2):
            BoundaryConditions(angle, z, v)
            GausSeidelIteration()
            GradientMagnitudes()

    print()
    SigmoidCharge(0.0, 50)
    fileName = "rectangleAttemptParallel.tns"
# char str[20];
# sprintf(str,"-%d.tns",vv);
# fileName += str;
    OutputTensor(fileName.c_str(), thresholdedChargeIntegral)
# }

#	Output("potential.dat", potential, -1, nodes/2);
#	Output("field.dat", field, -1, nodes/2);
#	Output("charge.dat", chargeIntegral, -1, nodes/2);
#
#	PrintChargeRange();
#
#	Control();

    print("Execution time: ", time.time() - startTime)

Run()
