import numpy as np
class droplet(object):
    #temperature is in K, SI units are used, droplet only can consist one kind of fuel
    #will use average properties if fuel blends are used
    def __init__(self, radius=50e-6, fuel_name='A2', temp=300, fraction=1):
        self.r = radius
        self.fuel = fuel_name
        self.T = temp #in kelvin
        self.frac = fraction

        #definition of fuel properties
        if self.fuel == 'A2':
            self.IBP = 159.2 + 273
            self.FBP = 270.5 + 273
            self.ABP = (self.IBP + self.FBP)/2
            self.sigma = 23.3*1e-3 #surface tension, J/m^2
            self.rho = 1018.1 - 0.7*self.T #linear interpolation, kg/m^3
            self.m = 4.0/3.0*self.rho*np.pi*(self.r)**3.0
            self.hfg = 0.428e6 #J/kg
            self.k = 0.08 #W/m-K
            self.cp = 2.6e3 #J/kg-K
            self.MW = 158.6
            self.boiled = False
            #why is viscosity not in here??
        if self.fuel == 'C1':
            
            pass

        if self.fuel == 'C5':
            pass
        
