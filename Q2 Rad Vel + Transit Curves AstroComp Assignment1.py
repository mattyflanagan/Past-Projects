# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 18:09:08 2024

@author: matth_yxw8uld
"""
#Q2
import numpy as np
import matplotlib.pyplot as plt

#extracting data
d = np.load('JSAstroLab2024_transit_data_21365396(1).npz')

#print('time',d['time'],'flux=',d['flux'],'P=',d['P'])
print(list(d.keys()))

t,f,P=d['time'],d['flux'],d['P']
print('no. of days for orbit=',P)


#plotting flux vs time
plt.figure()
plt.plot(t,f)
plt.xlabel('time (BJD)')
plt.ylabel('flux')
plt.title('Transit Data: Flux vs Time')
plt.show()

#this was printed to find the time stamps of the transit and the out of transit phases
'''for i in range(len(f)):
    if f[i]<0.9875:
        print('tnought part',t[i], i)
        
for i in range(len(f)):
    if 1.0010600892923953>f[i]>0.9875:
        print('limb part',f[i], i)'''

#Central transit time found as middle point of time in transit (from what was printed previously but now commented out)
Tnought=t[161]
print('Tnought=',Tnought,'BJD')
print('amount of time stamps:',len(f))

#calculating Delta F and rho:
#average of flux without transit
topavg=0
for i in range(123):
    topavg+=f[i]
for i in range(102):
    topavg+=f[198+i]

topavg=topavg/225
print('topavg=',topavg)

#average of flux during transit
bottomavg=0
for i in range(61):
    bottomavg+=f[130+i]

bottomavg=bottomavg/61
print('bottomavg=',bottomavg)

DelF=topavg-bottomavg         #change in flux during transit
rho=np.sqrt(np.abs(DelF))     #rho=square root of DelF
print('Delta Flux=',DelF)
print('rho=',rho)


parameters = (Tnought, 6.7, 0.12581726133458684, topavg, (np.radians(90)))  #Corrected variables (this took forever)

#zvalues list created and plotted, just to help debug previously (Tnought was originally taken as the length of time in transit t[191]-t[131])
zvalues=[]

def flux(par, time, P):
    #get parameters
    Tnought, a_Rstar,rho, foot, inc = parameters
    
    #calculate phase angle
    phi=((2 * np.pi / P) * (time - Tnought))
    
    #normalised separation
    b= a_Rstar * np.cos(inc)
    z=np.sqrt(((( (a_Rstar) * np.sin(phi)) )**2) + ( (b * np.cos(phi))**2));
   
    
    #calculate flux
    flux = np.ones(t.size)*foot         #out of transit: foot=topavg
    flux[z<=(1-rho)] = foot-(rho**2)    #flux valeus when fully in transit
    ind = (z > (1-rho)) * (z<=(1+rho))  #when in ingress/egress
    k0 = np.arccos(((rho**2) + (z[ind]**2) -1)/(2 * rho * z[ind])) # kappa nought
    k1 = np.arccos((1 - (rho**2) + (z[ind]**2)) / (2 * z[ind])) # kappa one
    flux[ind] = foot - ((1/np.pi) * ((k0 * (rho**2)) + k1 - np.sqrt((4*(z[ind]**2) -(1+(z[ind]**2) - rho**2)**2)/4 )))
    
    zvalues.append(z)
    return flux
flux_new = flux(parameters,t,P)

#plotting zvalues to help see what was wrong
plt.figure()
plt.title('Z-Values vs Time')
plt.plot(t,np.array(zvalues).flatten())
plt.xlabel('Time (BJD)')
plt.ylabel('Z-Values (a/R_*)')
plt.show()

#And finally, plotting the fit over the flux
plt.figure()
plt.plot(t, f, label='Observed Flux')
plt.axhline(y=bottomavg, linestyle='dotted',label='Flux in Full Transit')
plt.axhline(y=topavg, linestyle='dotted',label='Flux out of Transit',color='navy')
plt.plot(t, flux_new, label='Calculated Flux')
plt.xlabel('Time (BJD)')
plt.ylabel('Flux')
plt.title(' Fitted Transit Data: Flux vs Time')
plt.legend(loc=6)
plt.show()
