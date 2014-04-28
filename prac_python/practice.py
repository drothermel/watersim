#!/usr/bin/python

import sys
import math
import matplotlib.pyplot as plt

xc = float(sys.argv[1])
yc = float(sys.argv[2])
theta = math.radians(float(sys.argv[3]))
phi = math.radians(109.47)
docm = (1.0/9)*math.cos(phi/2)

oloc = (xc + docm*math.cos(theta), yc + docm*math.sin(theta))
haloc = (oloc[0]-math.cos(phi/2-theta), oloc[1]+math.sin(phi/2-theta))
hbloc = (oloc[0]-math.cos(phi/2+theta), oloc[1]-math.sin(phi/2+theta))

fig = plt.figure()
ax = fig.add_subplot(111)
x_points = [haloc[0], oloc[0], xc, oloc[0], hbloc[0]]
y_points = [haloc[1], oloc[1], yc, oloc[1], hbloc[1]]
p = ax.plot(x_points, y_points, 'b')
plt.ylim([-1,1])
plt.xlim([-1,1])
ax.set_xlabel('x-points')
ax.set_ylabel('y-points')
ax.set_title('Simple XY point plot')
fig.show()

raw_input()