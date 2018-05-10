from droplet import *
import numpy as np
#T is used to set temperature of the droplets, mdot is used to
#calculate mean droplet sizes
def entrainer(T, mdot):
    #assume everytime the kernel entrains 10 droplets, based on the fuel mass flow rate
    m10 = mdot*1e-6
    num = 10
    m = m10/num
    fuel = 'A2'
    drop_tempo = droplet(1, fuel, T, 1) #a "useless" droplet to get fuel info
    #calculate r based on mass put into the the kernel, right now, no statistically distribution of r
    r = np.power(3.0/4.0/drop_tempo.rho/np.pi*m,1.0/3.0)
    frac = 1
    
    d_arr = []
   # print(r)
    #create 10 droplets of the same size that give the proper mass flow rate
    for i in range(0,num):
        d_arr.append(droplet(r, fuel, T, frac))
    return d_arr
