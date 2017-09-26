#important note in this code I use this way of error propergation 
#dN=sqrt(N),dt=10ms/sqrt(12)
#dcount_rate=rate*sqrt((dN/N)^2+(dt/t)^2)
#av=(x1+x2)/2
#dav=0.5*sqrt(dx1^2+dx2^2)



import numpy as np
from numpy import cos,sqrt,pi,exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
print "##########################################################################################################"
print "##########################################################################################################"
print "This code is created by Jompoj Wongphechauxsorn, master student in Astrophysics Bonn 5th May 2017"
print "If you found any error feel free to contact me with jompoj.bjstp@gmail.com,jompoj@mpifr-bonn.mpg.de,s6jowong@uni-bonn.de"
print "##########################################################################################################"
print "##########################################################################################################"
E=-14.4*10**(3+9) #E_0=14 keV in neV
c_0=2.997*10**(8+3) # Speed of light in mm/s
#Warning!!!!!! makeing your initial value too precise will never fit with data
r=[0,12,24,33,40,51,65] #change these value
A=[11,5,4,5,7,15] #change these value
C=[-5,-3,-0.5,1,3,5]#change these value
W=1#change these value
S=135#change these value
color=['-r','-g','-b','-r','-g','-b']#change these value for difference colur 

#def lo(x, amp, cen, wid,shift):
#    return shift-(amp/(sqrt(pi/2)*wid)) * exp(-2*(x-cen)**2 /(wid**2))

def lo(x, amp, cen, wid,shift):
    return shift-((1.0/(np.pi*wid))*(amp/(1+((x-cen)/wid)**2)))
print "##########################################################################################################"
print "Please make sure that your file have columm 0 is velociy in mm,3 is count rate,6 is count rate's error"
print "##########################################################################################################"
data=np.loadtxt('data.txt') #load file data.txt that columm 0 is velociy in mm,3 is count rate,6 is count rate error

data=data.transpose() #for sake of simple when you call matrix


popt1, pcov1 = curve_fit(lo,data[0][r[0]:r[1]],data[3][r[0]:r[1]],p0=[A[0],C[0],W,S])
popt2, pcov2 = curve_fit(lo,data[0][r[1]:r[2]],data[3][r[1]:r[2]],p0=[A[1],C[1],W,S]) 
popt3, pcov3 = curve_fit(lo,data[0][r[2]:r[3]],data[3][r[2]:r[3]],p0=[A[2],C[2],W,S])
popt4, pcov4 = curve_fit(lo,data[0][r[3]:r[4]],data[3][r[3]:r[4]],p0=[A[3],C[3],W,S]) 
popt5, pcov5 = curve_fit(lo,data[0][r[4]:r[5]],data[3][r[4]:r[5]],p0=[A[4],C[4],W,S])
popt6, pcov6 = curve_fit(lo,data[0][r[5]:r[6]],data[3][r[5]:r[6]],p0=[A[5],C[5],W,S])
 #curve fitting 
#note that  we just put whatever number that make optimization finish

perr1 = np.sqrt(np.diag(pcov1))
perr2 = np.sqrt(np.diag(pcov2))
perr3 = np.sqrt(np.diag(pcov3))
perr4 = np.sqrt(np.diag(pcov4))
perr5 = np.sqrt(np.diag(pcov5))
perr6 = np.sqrt(np.diag(pcov6))
#calculate error for each parameter 

M=np.copy([popt1,popt2,popt3,popt4,popt5,popt6])  #make array of parameter 
dM=np.copy([perr1,perr2,perr3,perr4,perr5,perr6]) #make array of parameter's error
M=M.astype(float) #Change to float
dM=dM.astype(float) #Change to float
tM=M.transpose() #For sake of simple
tdM=dM.transpose() #For sake of simple
plt.legend()
plt.xlabel('speed(mm/s)') #x axis name
plt.ylabel('Counts') #y axis name

plt.errorbar(data[0],data[3],yerr=data[6],fmt='--o') #error bars 
plt.plot(data[0],data[3],"o",label='experiment') #data
for i in range(0,6):
	plt.plot(data[0][r[i]:r[i+1]+1],lo(data[0][r[i]:r[i+1]+1],M[i][0],M[i][1],M[i][2],M[i][3]),color[i],label="")
#fitting

plt.legend()
plt.show()
print "##########################################################################################################"
print "###############################################Parameter##################################################"
print "##########################################################################################################"
print "amp	damp	cen	dcen	width	dwidth	shift	dshift"

for i in range(0,5):
		print M[i][0],dM[i][0],M[i][1],dM[i][1],M[i][2],dM[i][2],M[i][3],dM[i][3]
print "##########################################################################################################"
print "##########################################################################################################"

iso=np.average(tM[1])
diso=(1/6.0)*np.sqrt(np.sum(tdM[1]**2))
g1_2=((tM[1][3]-tM[1][1])+(tM[1][4]-tM[1][2]))/2.0
dg1_2=0.5*np.sqrt(0.5*np.sqrt(tdM[1][3]**2+tdM[1][1]**2))**2+(0.5*np.sqrt(tdM[1][4]**2+tdM[1][2]**2))
g3_2=((tM[1][1]-tM[1][0])+(tM[1][2]-tM[1][1])+(tM[1][4]-tM[1][3])+(tM[1][5]-tM[1][4]))/4.0
dg3_2=0.5*np.sqrt(0.5*np.sqrt(tdM[1][1]**2+tdM[1][0]**2))**2+(0.5*np.sqrt(tdM[1][2]**2+tdM[1][1]**2)+(0.5*np.sqrt(tdM[1][4]**2+tdM[1][3]**2))**2+(0.5*np.sqrt(tdM[1][5]**2+tdM[1][4]**2)))

print "Isomeric shift=",iso*E/c_0,"$\pm$",diso*E/c_0,"neV"
print "$g_{1/2}$=",g1_2*E/c_0,"$\pm$",dg1_2*E/c_0,"neV"
print "$g_{3/2}$=",g3_2*E/c_0,"$\pm$",dg3_2*E/c_0,"neV"


