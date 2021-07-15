import numpy as np
import mcubes


def ReadTNS(fileName):
 file = open(fileName, "r")
 tensorText = file.read()
 tensorText = tensorText.split()
 xMax = int(tensorText[0])
 yMax = int(tensorText[1])
 zMax = int(tensorText[2])
 count = 5
 tensor = np.ndarray(shape=(xMax, yMax, zMax), dtype=float, order='F')
 for z in range(zMax):
  for y in range(yMax):
   for x in range(xMax):
    tensor[x][y][z] = float(tensorText[count])
    count += 1
 return tensor


tensor = ReadTNS("../Simulation/file-mediated.tns")

# Extract the isosurface
vertices, triangles = mcubes.marching_cubes(tensor, 0.1)

# Export the result
mcubes.export_mesh(vertices, triangles, "print.dae", "print")


