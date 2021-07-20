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

def MakeDAE(tnsFileName, daeFileName):
 tensor = ReadTNS(tnsFileName)
# Extract the isosurface
 vertices, triangles = mcubes.marching_cubes(tensor, 0.5)
# Export the result
 mcubes.export_mesh(vertices, triangles, daeFileName, "print")
 print("DAE file exported")


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

def RandomZ():
 z = random.random()
 z = 50.0*random.random()*(0.5 + 0.5*maths.cos(z*2.0*maths.pi))
 z = 1 + int(z)
 return z

def RandomVoltage():
 return random.random()*50.0


def CreateVoltages(fileName, voltages):
 boundary = LoadBoundary("boundaryNodes.txt")
 file = open(fileName, "w")
 for l in range(voltages):
  b = RandomBoundaryPoint(boundary)
  s = str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ()) + ' ' + str(RandomVoltage()) + ' '
  b = RandomBoundaryPoint(boundary)
  s = s + str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ()) + ' ' + str(-RandomVoltage()) + '\n'
  file.write(s)
 file.write("-1 -1 -1 -1 -1 -1 -1 -1\n")
 print("Voltage file created")
 file.close()

def RunASimulation(name, voltages):
 CreateVoltages(name + ".v", voltages)
 subprocess.run(["./Electric3DPrinter", (name + ".v"), (name + ".tns")])
 print("C++ simulation run")
 MakeDAE(name + ".tns", name + ".dae")


random.seed(a=None, version=2)
RunASimulation("t2", 200)
