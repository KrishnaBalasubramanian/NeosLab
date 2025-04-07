#################################################
###			 NEOSLAB Helper Codes ###########
#################################################
from scipy.optimize import curve_fit
import numpy as np
import pdb
import matplotlib.pyplot as plt
import os
import re
class NeosLab:
	fitCounter=0
	def customFit(xData,yData,fitType,sp,BG=0,debug=False):
		"""
			Codes to fit known functions. 
				xData: Send only the x points that need fitting
				yData: Send only the y points to the extent of fitting. Keep it closest to the curve
				fitType: Should be one of the following:
						'1Gauss' - Single gaussian, 
						'2Gauss' - Two Gaussians,
						'1Lor' - Single Lorentzian
						'2Lor' - Two Lorentzian
				sp: Should be a dictionary array. Each Dict element should have the following 
					"cen","amp","sigma"

				returns an array of fitting values and error
		"""
		def background(x,m,c):
			return m*x + c
		def _1gaussian(x, amp1,cen1,sigma1):
			return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2)))
		def _2gaussian(x, amp1,cen1,sigma1, amp2,cen2,sigma2):
			return amp1*(1/(sigma1*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen1)/sigma1)**2))) + \
				amp2*(1/(sigma2*(np.sqrt(2*np.pi))))*(np.exp((-1.0/2.0)*(((x-cen2)/sigma2)**2)))
		def _2Lorentzian(x, amp1, cen1, wid1, amp2,cen2,wid2):
			return (amp1*wid1**2/((x-cen1)**2+wid1**2)) +\
				(amp2*wid2**2/((x-cen2)**2+wid2**2))
		def _1Lorentzian(x, amp1, cen1, wid1):
			return (amp1*wid1**2/((x-cen1)**2+wid1**2)) 
		try:
			############### first fit a simple background ##############
			popt, pcov = curve_fit(background, xData, yData, p0=[0,BG])
			yCorr = yData - background(xData,*popt)
			if fitType =='1Gauss':
				popt =[0,0,0]
				pcov = []
				popt, pcov = curve_fit(_1gaussian, xData, yCorr, p0=[sp[0]["amp"], sp[0]["cen"], sp[0]["sigma"]])
			elif fitType =='2Gauss':
				popt =[0,0,0,0,0,0]
				pcov = []
				popt, pcov = curve_fit(_2gaussian, xData, yCorr, p0=[sp[0]["amp"], sp[0]["cen"], sp[0]["sigma"], sp[1]["amp"], sp[1]["cen"], sp[1]["sigma"]])
			elif fitType =='1Lorentzian':
				popt =[0,0,0]
				pcov = []
				popt, pcov = curve_fit(_1Lorentzian, xData, yCorr, p0=[sp[0]["amp"], sp[0]["cen"], sp[0]["sigma"], sp[1]["amp"]])
			elif fitType =='2Lorentzian':
				popt =[0,0,0,0,0,0]
				pcov = []
				popt, pcov = curve_fit(_2Lorentian, xData, yCorr, p0=[sp[0]["amp"], sp[0]["cen"], sp[0]["sigma"], sp[1]["amp"], sp[1]["cen"], sp[1]["sigma"]])
			else:
				popt =[0,0,0,0,0,0]
				pcov = []
				popt, pcov = curve_fit(_2gaussian, xData, yCorr, p0=[sp[0]["amp"], sp[0]["cen"], sp[0]["sigma"], sp[1]["amp"], sp[1]["cen"], sp[1]["sigma"]])
			NeosLab.fitCounter +=1
			if debug:
				fig,ax = plt.subplots()
				ax.plot(xData,yCorr)
				ax.plot(xData,_2gaussian(xData,*popt),'r')
				os.makedirs('debug',exist_ok=True)
				fig.savefig('debug/fitting'+str(NeosLab.fitCounter) + '.png')
				plt.close(fig)
		except Exception as E:
			print(E)
			print('Fitting errors')
		return [popt,pcov]

	def readWitecMappedData(fileName,skipRows,xSteps,ySteps):
		"""
			A function to read Witec mapped data. 
			filename: 
			skip rows: 
			xSteps: number of steps in the x direction
			ySteps: NUmber of steps in the y direction
			
			returns the scan positions in 1/cm and the the data in a 3D matrix data[x,y,:] 


		"""
		file = open(fileName,'r', encoding='utf-8', errors='ignore')
		lines = file.readlines()
		rows = []
		rc=0
		for line in lines:
			if rc > skipRows:
				cols = []
				nums = re.findall('([-+]?\d+\.\d+[eE]?[-+]?\d+)',line)
				if not nums == []:
					for n in nums:
						cols.append(float(n))
					rows.append(cols)
			rc+=1 ### count rows
		inData = np.zeros((xSteps,ySteps,len(rows)))
		xAxis = np.zeros(len(rows))
		for i in range(len(rows)):
			xAxis[i]= rows[i][0]
			c=1
			for y in range(ySteps):
				for x in range(xSteps):
					try:
						inData[x,y,i] = rows[i][c]
					except Exception as E:
						inData[x,y,i] = 0
						#print(E)
						#print('mostly missing data')
					c+=1
		return [xAxis,inData]


    def getTheta(rin,win):
        phi = np.zeros(len(win))
        phi_err = np.zeros(len(win))
        rFun = interp1d(win,rin,bounds_error = False, fill_value=(rin[0],rin[-1]))
        Integrand = lambda w,wp: np.log(rFun(w))/(w + wp)
        #pdb.set_trace()
        for i in range(len(win)):
            phi[i],phi_err = integrate.quad(Integrand,0,10, args = win[i],weight = 'cauchy',wvar = win[i],limit = 100,epsabs = 1E-20)#*2*win[i]/np.pi
            print(phi_err)
            #n[i] = (1 - rin[i]**2)/(1 - 2*rin[i]*np.cos(phi[i]) + rin[i]**2)
            #k[i] = (-2*rin[i]*np.sin(phi[i]))/(1 - 2*rin[i]*np.cos(phi[i]) + rin[i]**2)
        return phi*2*win/np.pi 



    def getnk(rin,theta):
        n = np.zeros(len(rin))
        k = np.zeros(len(rin))
        for i in range(len(rin)):
            n[i] = (1 - rin[i]**2)/(1 - 2*rin[i]*np.cos(phi[i]) + rin[i]**2)
            k[i] = (-2*rin[i]*np.sin(phi[i]))/(1 - 2*rin[i]*np.cos(phi[i]) + rin[i]**2)
        return [n,k]

