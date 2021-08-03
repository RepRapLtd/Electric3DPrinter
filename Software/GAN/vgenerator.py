# Pattern: 5 16 1 51.5375 45 34 1 -51.5375
# X Y Z V X Y Z V
# NB X Y must be on the circular perimeter

import math as maths
import numpy as np
import mcubes
import random
import subprocess

def ReadTNS(fileName):
 file = open(fileName, "r")
 tensorText = file.read()
 tensorText = tensorText.split()
 xMax = int(tensorText[0])
 yMax = int(tensorText[1])
 zMax = int(tensorText[2])
 vMin = float(tensorText[3])
 vMax = float(tensorText[4])
 count = 5
 tensor = np.ndarray(shape=(xMax, yMax, zMax), dtype=float, order='F')
 for z in range(zMax):
  for y in range(yMax):
   for x in range(xMax):
    v = float(tensorText[count])
    tensor[x][y][z] = (v - vMin)/(vMax - vMin)
    count += 1
 return tensor


def ExportMesh(tnsFileName, meshFileName):
 tensor = ReadTNS(tnsFileName)
# Extract the isosurface
 vertices, triangles = mcubes.marching_cubes(tensor, 0.5)
# Export the result
 mcubes.export_obj(vertices, triangles, meshFileName)
 print("Mesh file exported")


def LoadBoundary(fileName):
 boundary = []
 file = open(fileName, "r")
 boundaryText = file.read()
 boundaryText = boundaryText.split()
 boundaryCount = int(boundaryText[0])
 count = 1;
 for bc in range(boundaryCount):
  x = int(boundaryText[count])
  count += 1
  y = int(boundaryText[count])
  count += 1
  boundary.append([x, y])
 return boundary

def RandomBoundaryPoint(boundary):
 return random.choice(boundary)

def RoughlyOppositeBoundaryPair(boundary):
 len = boundary.__len__()
 a = random.randint(0, len-1)
 b = random.randint(0, round(len/5))
 b = b - round(len/10)
 b = a + round(len/2) + b
 b = b%len
 return (boundary[a], boundary[b])
 #return(RandomBoundaryPoint(boundary), RandomBoundaryPoint(boundary))

# This attempts to put more electrodes at the Z extremes, to try to reduce end
# effects and make a blob in the middle. Sort of works...

def RandomZ(bottom, span):
 y = random.random()
 if y < 0.5:
  x = 2*y*y
 else:
  x = 0.5 + maths.sqrt((y - 0.5)*0.5)
 return bottom + int(x*span)

def RandomVoltage():
 return random.random()*50.0

def TopAndBottomFill(file, boundary, endCount, endDepth):
 for l in range(endCount):
  b = RoughlyOppositeBoundaryPair(boundary)
  b0 = b[0]
  b1 = b[1]
  s = str(b0[0]) + ' ' + str(b0[1]) + ' ' + str(RandomZ(1, endDepth)) + ' ' + str(RandomVoltage()) + ' '
  s = s + str(b1[0]) + ' ' + str(b1[1]) + ' ' + str(RandomZ(1, endDepth)) + ' ' + str(-RandomVoltage()) + '\n'
  file.write(s)
 for l in range(endCount):
  b = RoughlyOppositeBoundaryPair(boundary)
  b0 = b[0]
  b1 = b[1]
  s = str(b0[0]) + ' ' + str(b0[1]) + ' ' + str(RandomZ(51 - endDepth, endDepth)) + ' ' + str(RandomVoltage()) + ' '
  b = RandomBoundaryPoint(boundary)
  s = s + str(b1[0]) + ' ' + str(b1[1]) + ' ' + str(RandomZ(51 - endDepth, endDepth)) + ' ' + str(-RandomVoltage()) + '\n'
  file.write(s)

def CreateVoltages(fileName, voltages, endCount, endDepth):
 boundary = LoadBoundary("boundaryNodes.txt")
 file = open(fileName, "w")
 TopAndBottomFill(file, boundary, endCount, endDepth)
 for l in range(voltages):
  b = RandomBoundaryPoint(boundary)
  s = str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ(1, 51)) + ' ' + str(RandomVoltage()) + ' '
  b = RandomBoundaryPoint(boundary)
  s = s + str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ(1, 51)) + ' ' + str(-RandomVoltage()) + '\n'
  file.write(s)
 file.write("-1 -1 -1 -1 -1 -1 -1 -1\n")
 print("Voltage file created")
 file.close()

def RunASimulation(name, voltages, endCount, endDepth):
 CreateVoltages(name + ".v", voltages, endCount, endDepth)
 subprocess.run(["./Electric3DPrinter", "-i", name + ".v", "-o", name + ".tns", "-so", "0.0", "-sm", "50.0"])
 print("C++ simulation run")
 ExportMesh(name + ".tns", name + ".obj")


random.seed(a=None, version=2)
RunASimulation("t2", 1, 100, 1)
