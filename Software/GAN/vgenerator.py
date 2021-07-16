# Pattern: 5 16 1 51.5375 45 34 1 -51.5375
# X Y Z V X Y Z V
# NB X Y must be on the circular perimeter

import random

lines = 4118

boundary = []

def LoadBoundary(fileName):
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

def RandomBoundaryPoint():
 return random.choice(boundary)

def RandomZ():
 return random.randrange(1, 50)

def RandomVoltage():
 return random.random()*50.0




LoadBoundary("boundaryNodes.txt")

random.seed(a=None, version=2)
file = open("randomVoltages.v", "w")
for l in range(lines):
 b = RandomBoundaryPoint()
 s = str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ()) + ' ' + str(RandomVoltage()) + ' '
 b = RandomBoundaryPoint()
 s = s + str(b[0]) + ' ' + str(b[1]) + ' ' + str(RandomZ()) + ' ' + str(-RandomVoltage()) + '\n'
 file.write(s)
file.write("-1 -1 -1 -1 -1 -1 -1 -1\n")
file.close()
