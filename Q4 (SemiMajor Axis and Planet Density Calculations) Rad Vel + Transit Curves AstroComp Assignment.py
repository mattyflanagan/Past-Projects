# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 17:13:16 2024

@author: matth_yxw8uld
"""

import numpy as np

#radius and mass of our sun in SI
R_sunm = 6.957e8    #radius in metres
M_sunkg = 1.989e30  # Solar mass in kg

R_star=0.669*R_sunm
earthmass=5.972e+24 #earth mass in kg
G = 6.67430e-11 
P = 1.2410229                 #period of the planet (days).
M_star = 0.69*M_sunkg
Mpsini=2.0130239385819713*earthmass
topavg=0

#finding planet radius from rho and R_star
d = np.load('JSAstroLab2024_transit_data_21365396(1).npz')
f=d['flux']
for i in range(125):
    topavg+=f[i]
for i in range(100):
    topavg+=f[200+i]

topavg=topavg/225
print('topavg=',topavg)

#average of flux during transit
bottomavg=0
for i in range(75):
    bottomavg+=f[125+i]

bottomavg=bottomavg/75
print('bottomavg=',bottomavg)

DelF=topavg-bottomavg #change in flux during transit
rho=np.sqrt(DelF)

R_planet=R_star/rho

#calculating density and semimajor axis

density = Mpsini / ((4/3) * np.pi * R_planet**3)

#Calculated the semi-major axis of the orbit.
smaxis = ((G * (M_star + Mpsini)) / (4 * np.pi**2))**(1/3) * P**(2/3)

#Printed results.
print("Density of planet=", density, "kg/m^3")
print("Semi-major axis=", smaxis, "m")

