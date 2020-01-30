import numpy as np #import NumPy module
import astropy.units as u #import AstroPy module

def Read(filename): #Define a function that takes a file as input
    
    file = open(filename,'r') #opening file

    line1 = file.readline() #reading the first line of MW_000.txt
    label, value = line1.split()
    time = float(value)*u.Myr #store time in units of Myr 

    line2 = file.readline() #reading the second  line of MW_000.txt
    label, value = line2.split()
    total = float(value)# total = total number of particles

    file.close() #close the file

    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)#store remainder of the file starting from the 4th line and onwards
    
    print("Particle Type: ",data['type'][1]) #check: printing 1st column 2nd entry
    
    return time, total, data #returning time, total number of particles, and data array


Read(filename = "MW_000.txt") #calling the function with filename = MW_000.txt 

    
    
