import matplotlib.pyplot as plt
import numpy as np
import re
import pdb
from NeosLab import NeosLab as NL
fileName = '/Users/krishnabalasubramanian/Downloads/Gb_LA.mat'
xLim = 55E-6
yLim = 60E-6
xSteps = 35
ySteps = 35
fitType = '2Gauss'
xStart = 360
xStop = 440
p1 = {"amp":5,"cen": 380,"sigma":5}
p2 = {"amp":5,"cen": 400,"sigma":5}
skipRows = 18
xAxis,inData = NL.readWitecMappedData(fileName,skipRows,xSteps=xSteps,ySteps=ySteps)
xVal=np.zeros((xSteps,ySteps))
yVal=np.zeros((xSteps,ySteps))
zVal=np.zeros((xSteps,ySteps))

for i in range(xSteps):
	for j in range(ySteps):
		xVal[i,j] = i*xLim/xSteps
		yVal[i,j] = j*yLim/ySteps
		xStartIx=np.searchsorted(xAxis,xStart)
		xStopIx=np.searchsorted(xAxis,xStop)
		popt,pconv = NL.customFit(xAxis[xStartIx:xStopIx],inData[i,j,xStartIx:xStopIx],fitType,[p1,p2])
		############### plot only the first peak position
		zVal[i,j] = popt[0] 
		############### plot peak seperation
		#zVal[i,j] = popt[4] - popt[1] 

fig,ax = plt.subplots()
pc=ax.pcolormesh(zVal,vmin = 0,vmax = 50)
fig.colorbar(pc)
plt.show()
