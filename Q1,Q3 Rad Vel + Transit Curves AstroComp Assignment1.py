# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 20:55:43 2024

@author: matth_yxw8uld
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import correlate
from scipy.optimize import curve_fit

#Downloaded the RV Data.
d = np.load("JSAstroLab2024_rv_data_21365396.npz")

#To find how many spectra are in the data set.
#print(list(d.keys()))
no_spectra = 40
print("Number of spectra=", no_spectra)

#Extracting relevant info.
wavelength = d['wavelength']
v = d['velocity']
spectemp = d['spec_template']

times = []
spectra = []

for i in range(40):
    times_key = f'time_BJD_{i+1}'
    ts = d[times_key]
    spectrum_key = f'spectrum_{i+1}'
    spectrum = d[spectrum_key]

    times.append(ts)
    spectra.append(spectrum)
    
def cross_correlate(spectemp, spectrum):
    correlresult = correlate(spectrum - np.mean(spectrum), spectemp - np.mean(spectemp), mode='same')
    return correlresult

def gaussian(x, a, b, data):
    return a * np.exp(-(x - b)**2 / (2 * data**2))

rvcurve = np.zeros(no_spectra)
def calc_rv_curve(spectemp, spectra, wavelength):
    for i in range(40):
        spectra1 = spectra[i]
        correlresult = cross_correlate(spectemp, spectra1)

        #Using curve_fit to fit the Gaussian model to the cross-correlation function.
        try:
            popt, pcov = curve_fit(gaussian, wavelength, correlresult, p0=[np.max(correlresult), wavelength[np.argmax(correlresult)], 1.0])
        except RuntimeError:
            popt = [0, 0, 1]
        a_fit, b_fit, data_fit = popt
        rv_shift = b_fit - wavelength[len(wavelength) // 2]
        rvcurve[i] = rv_shift #Employs the velocity shifts to the rv curve.

    return rvcurve - np.mean(rvcurve)

gauss_rv = calc_rv_curve(spectemp, spectra, wavelength)
#print("Radial Velocity Shifts:",gauss_rv)

t_list = [float(time) for time in times]

#plotting spectrum template
plt.figure()
plt.plot(wavelength, spectemp, label = 'Spectral Template')
plt.title("Spectrum Template")
plt.ylabel("Intensity")
plt.xlabel("Wavelength (Angstrom)")
plt.legend()
plt.show()

#Plotting RV curve
plt.figure()
plt.plot(t_list, gauss_rv, marker='o', linestyle='-')
plt.xlabel('Time (BJD)')
plt.ylabel('Radial Velocity Shift (km/s)')
plt.title('RV Curve Data (RV Shift vs Time)')
plt.grid(True)
plt.show()

#Now calculating velocity shift...
def calc_vel_shift(vel, spectrum):
    correl = correlate(vel - vel.mean(), spectrum - spectrum.mean(), mode='same')
    shift = v[np.argmax(correl)]
    return shift

#velocity shifts for all spectra
velocity_shifts = [calc_vel_shift(spectemp, spectrum) for spectrum in spectra]

#plotting after this cross correlation.
plt.figure()
for spectrum in spectra[:40]:
    plt.plot(v, correlate(spectrum - spectrum.mean(), spectemp - spectemp.mean(), 'same'))

plt.title("Zoomed Cross Correlation vs Velocity Shifts")
plt.xlabel("δV (km/s) ")
plt.ylabel("Cross Correlation Value")
plt.xlim([-0.5e+6,0.5e+6])
plt.show()

plt.figure()
plt.title("Cross Correlation vs Velocity Shifts")
plt.xlabel("δV (km/s) ")
plt.ylabel("Cross Correlation Value")
for spectrum in spectra[:40]:
    plt.plot(v, correlate(spectrum - spectrum.mean(), spectemp - spectemp.mean(), 'same'))
plt.show()

#######################################################################################
#Part 3
#Defined parameters:
Ms = 0.69           #Stellar mass in solar masses
Rs = 0.669          #Stellar radius in solar radii
Ts = 4300           #Stellar temperature in Kelvin
P = 1.2410229       #Orbital period in days.
G = 6.67430e-11     #Gravitational constant (m^3 kg^-1 s^-2).
Ms = 0.69           #Stellar mass in solar masses.

#defining new variables for calculations (Unit conversions needed as G is in Nm^2/kg^2)
SolarMasses=1.989e+30 #Solar mass in Kg
MsKG=Ms*SolarMasses   #Star mass in Kg
inc = np.radians(90)  #orbital inc from part 2 
Psec= 107224.37856    #orbital period in seconds 
inc=np.radians(90)    #angle of inclination from part 2
print('MsKG=',MsKG)

def sinfunc(Kstar, t, offset,phi,A):
    return Kstar*np.sin(-(A * np.pi * t / Psec + phi)) + offset

Sin_Fit=[]
for i in range(len(times)):
    phi=((2 * np.pi / P) * (times[i] - 2459639.5417019445))
    Sin_Fit.append(sinfunc(0.0013,times[i],-0.0001, phi, 2.005))
    
plt.figure()
plt.plot(t_list, gauss_rv,marker='o', linestyle='-', label='RV Curve')
plt.plot(t_list, Sin_Fit, label='Sin Fit',color='red')
plt.legend()
plt.xlabel('Time (BJD)')
plt.ylabel('Radial Velocity Shift(km/s)')
plt.title('RV Curve Data (RV Shift vs Time)')
plt.grid(True)
plt.show()


#Plotting phase-folded RV curve.
phase = np.mod(t_list, P) / P  #Normalized phase.
phase.sort()

plt.figure(figsize=(10, 6))
plt.scatter(phase, gauss_rv, label='Data')            #plotting phase folded data
plt.plot(phase,Sin_Fit, color='red',label='Sin Fit')  #plotting fit
plt.xlabel('Phase')
plt.legend()
plt.ylabel('Radial Velocity Shift(km/s)')
plt.title('Phase-folded Radial Velocity Curve')
plt.grid(True)
plt.show()

#Finding Minimum Planet Mass
Kstar=0.0013      #from our sin fit
Mp = (Kstar * ((MsKG)**(2/3)) * P**(1/3)) / ((2*np.pi*G)**1/3) #assumed that Mp+Ms roughly equals Ms, as Ms>>Mp
print('Minimum Mass of the Planet=', Mp,'Kg')
  

