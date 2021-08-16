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

def RandomBoundaryIndex(boundary):
 return random.randrange(0, boundary.__len__())

# Return two indices into the boundary array that are roughly opposite
# TODO - remove constants

def RoughlyOppositeBoundaryIndices(boundary):
 len = boundary.__len__()
 a = random.randrange(0, len)
 b = random.randint(0, round(len/5))
 b = b - round(len/10)
 b = a + round(len/2) + b
 b = b%len
 return (a, b)

# This attempts to use more electrodes at the Z extremes, to try to reduce end
# effects and make a blob in the middle. Sort of works...

def RandomZ(bottom, span):
# y = random.random()
# if y < 0.5:
#  x = 2*y*y
# else:
#  x = 0.5 + maths.sqrt((y - 0.5)*0.5)
# return bottom + int(x*span)
 return random.randint(bottom, bottom+span)

def RandomVoltage():
 return 50 + random.random()*50.0

def ElectrodePairToFile(file, ends):
 s = ""
 for end in range(2):
  for number in range(4):
   sep = ' '
   if end == 1 and number == 3:
    sep = '\n'
   s += str(ends[end][number]) + sep
 file.write(s)

def TopAndBottomFill(file, boundary, endCount, endDepth):
 for l in range(endCount):
  oppositeBoundaryIndices = RoughlyOppositeBoundaryIndices(boundary)
  boundaryPoint0 = boundary[oppositeBoundaryIndices[0]]
  boundatryPoint1 = boundary[oppositeBoundaryIndices[1]]
  end = []
  ends = []
  end.append(boundaryPoint0[0])
  end.append(boundaryPoint0[1])
  end.append(RandomZ(1, endDepth))
  end.append(RandomVoltage())
  ends.append(end)
  end = []
  end.append(boundatryPoint1[0])
  end.append(boundatryPoint1[1])
  end.append(RandomZ(1, endDepth))
  end.append(-RandomVoltage())
  ends.append(end)
  ElectrodePairToFile(file, ends)
 for l in range(endCount):
  oppositeBoundaryIndices = RoughlyOppositeBoundaryIndices(boundary)
  boundaryPoint0 = boundary[oppositeBoundaryIndices[0]]
  boundatryPoint1 = boundary[oppositeBoundaryIndices[1]]
  end = []
  ends = []
  end.append(boundaryPoint0[0])
  end.append(boundaryPoint0[1])
  end.append(RandomZ(51 - endDepth, endDepth))
  end.append(RandomVoltage())
  ends.append(end)
  end = []
  end.append(boundatryPoint1[0])
  end.append(boundatryPoint1[1])
  end.append(RandomZ(51 - endDepth, endDepth))
  end.append(-RandomVoltage())
  ends.append(end)
  ElectrodePairToFile(file, ends)

def MakeHole(file, boundary, voltages, diameter):
 len = boundary.__len__()
 b = RoughlyOppositeBoundaryIndices(boundary)
 z0 = RandomZ(1, 51)
 z1 = RandomZ(1, 51)
 for v in range(voltages):
  b0 = (b[0] + random.randrange(-diameter, diameter)) % len
  b1 = (b[1] + random.randrange(-diameter, diameter)) % len
  z00 = (z0 + random.randrange(-diameter, diameter)) % 51
  z11 = (z1 + random.randrange(-diameter, diameter)) % 51
  end = []
  ends = []
  end.append(boundary[b0])
  end.append(boundary[b1])
  end.append(RandomZ(1, 51))
  end.append(RandomVoltage())
  ends.append(end)
  b = RandomBoundaryIndex(boundary)
  end = []
  end.append(b[0])
  end.append(b[1])
  end.append(RandomZ(1, 51))
  end.append(-RandomVoltage())
  ends.append(end)
  ElectrodePairToFile(file, ends)



def CreateVoltages(fileName, voltages, endCount, endDepth):
 boundary = LoadBoundary("boundaryNodes.txt")
 file = open(fileName, "w")
 TopAndBottomFill(file, boundary, endCount, endDepth)
 for l in range(voltages):
  boundaryIndex = RandomBoundaryIndex(boundary)
  end = []
  ends = []
  end.append(boundary[boundaryIndex][0])
  end.append(boundary[boundaryIndex][1])
  end.append(RandomZ(1, 51))
  end.append(RandomVoltage())
  ends.append(end)
  boundaryIndex = RandomBoundaryIndex(boundary)
  end = []
  end.append(boundary[boundaryIndex][0])
  end.append(boundary[boundaryIndex][1])
  end.append(RandomZ(1, 51))
  end.append(-RandomVoltage())
  ends.append(end)
  ElectrodePairToFile(file, ends)
 file.write("-1 -1 -1 -1 -1 -1 -1 -1\n")
 print("Voltage file created")
 file.close()

def RunASimulation(name, voltages, endCount, endDepth):
 CreateVoltages(name + ".v", voltages, endCount, endDepth)
 subprocess.run(["./Electric3DPrinter", "-i", name + ".v", "-o", name + ".tns", "-so", "0.0", "-sm", "50.0"])
 print("C++ simulation run")
 ExportMesh(name + ".tns", name + ".obj")


random.seed(a=97)
RunASimulation("t2", 200, 1000, 5)
