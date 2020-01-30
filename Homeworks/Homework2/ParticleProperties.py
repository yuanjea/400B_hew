import numpy as np #import NumPy module
import astropy.units as u #import AstroPy module
from ReadFile import Read #import module to call Read from Readfile

def ParticleInfo(filename,ptype,pnumber): #Define a funciton that takes filename,ptype=particle type, and pnumber=particle number as inputs
    
     data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)#read file again via ReadFile store remainder of the file starting from the 4th line and onwards

     index = np.where(data['type'] == ptype) # create an index which takes data that has type = 2 (DiskStars)

     #calculating 3D distance
     r = np.sqrt(data[index]['x'][pnumber]**2 + data[index]['y'][pnumber]**2 + data[index]['z'][pnumber]**2)*u.kpc #r = sqrt(x^2 + y^2 + z^2) * kpc
     #calculating 3D velocity
     v = np.sqrt(data[index]['vx'][pnumber]**2 + data[index]['vy'][pnumber]**2 + data[index]['vz'][pnumber]**2)*u.km/u.s #v = sqrt(vx^2 + vy^2 +vz^2) * km/s
     #Retrieving Mass data
     m = data[index]['m'][pnumber]*1e10*u.Msun #units of 10^10 mass of the sun
    
     return r, v, m #returning distance,r ; velocity,v ; mass,m 

r, v, m = ParticleInfo('MW_000.txt',2,99) #calling ParticleInfo function with inputs as filename=MW_000.txt ; ptype=2 ; pnumber=99 (0 to 99)

R = r.to(u.lyr) #assigning distance,r units to lightyears using "to" function

print("3D distance: ", np.around(r,3)) #print 3D distance value in 3 decimal places and kpc units to screen
print("3D distance: ", np.around(R,3)) #print 3D distance value in 3 decimal places and lyr units to screen
print("3D velocity: ", np.around(v,3)) #print 3D velocity value in 3 decimal places and km/s units to screen
print("Mass: ", m) #print Mass values to screen 

    

    
