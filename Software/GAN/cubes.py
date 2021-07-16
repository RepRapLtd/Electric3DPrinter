import numpy as np
import mcubes


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


tensor = ReadTNS("randomShape.tns")

# Extract the isosurface
vertices, triangles = mcubes.marching_cubes(tensor, 0.5)

# Export the result
mcubes.export_mesh(vertices, triangles, "random1.dae", "print")


