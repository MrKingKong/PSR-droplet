import cantera as ct
import numpy as np
import matplotlib.pyplot as plt
from droplet import *
from entrainer import *
from vaporizer_temp import *
ct.suppress_thermo_warnings()

 

def execute():
    #==============================================================================
#     set up the cantera gases and reactors
#==============================================================================
    drop = droplet(50e-6, 'A2', 300, 1) #drop is defined only to use the drop properties
    #environment air
    gas_env = ct.Solution('a2.cti')
    gas_env.TPX = 300, ct.one_atm, 'O2:29,N2:71'

    gas_fuel = ct.Solution('a2.cti')
    gas_fuel.TPX = drop.ABP, ct.one_atm, 'POSF10325:1' #gaseous fuel at average boiling temperature

    gas_kernel = ct.Solution('a2.cti')
    gas_kernel.TPX = 4000, ct.one_atm, 'O2:29,N2:71'
    gas_kernel.equilibrate('HP')
    
    res_air = ct.Reservoir(gas_env)
    res_fuel = ct.Reservoir(gas_fuel)
    kernel = ct.IdealGasConstPressureReactor(gas_kernel)
    kernel.volume = 2.4e-8 #m^3

    mfc1 = ct.MassFlowController(res_air, kernel) #connect air reservoir to kernel
    mfc2 = ct.MassFlowController(res_fuel, kernel) #connect fuel reservoir to kernel

    mfc1.set_mass_flow_rate(6e-5)
    mfc2.set_mass_flow_rate(0)

    w = ct.Wall(kernel, res_fuel)
    w.set_heat_flux(0)
    w.expansion_rate_coeff = 1e8

    net = ct.ReactorNet({kernel})
    net.atol = 1e-10 #absolute tolerance
    net.rtol = 1e-10 #relative tolerance
#==============================================================================
#   endtime, time steps, and data to be stored    
#==============================================================================
    dt = 1e-6 #change time step for advancing the reactor
    endtime = 1000e-6
    mdot_fuel = 1.98e-6 #this number should be calculated based on equivalence ratio, kg/s
    
    droplets = []
    T = []
    time = []
    m = []
    m_dot = []
    q_dot = []
    droplet_count = []
    X_fuel = []
    X_O = []
    X_CO2 = []
    num_d = []
    
#==============================================================================
#     advance reaction
#==============================================================================
    print("Running simulation")
    for i in range(1,int(endtime/dt)):
        print (int(endtime/dt) - i)
        #entrain droplets every 1 microsecond
        droplets.extend(entrainer(400, 5.*mdot_fuel))  
        mdot, qdot = vaporize(gas_kernel.T, droplets, dt)
        
#        print(gas_kernel.T, mdot, qdot, len(droplets))
#        print(str(mdot) + str(qdot) + str(gas_kernel.T)+'  '+str(len(droplets)))
        mfc2.set_mass_flow_rate(mdot)
        w.set_heat_flux(qdot) #heat required to heat up the droplets during that time step
        
        net.advance(i*dt)
        
        num_d.append(len(droplets))
        #storing variables
        T.append(gas_kernel.T)
        time.append(net.time)
        m.append(kernel.mass)
        m_dot.append(mdot)
        q_dot.append(qdot)
        droplet_count.append(len(droplets))
        if gas_kernel.X[0] < 0:
            gas_kernel.X[0] = 0
        X_fuel.append(gas_kernel.X[gas_kernel.species_index("POSF10325")])
        X_O.append(gas_kernel.X[gas_kernel.species_index("O")])
        X_CO2.append(gas_kernel.X[gas_kernel.species_index("CO2")])
        
#==============================================================================
#   plotting important information      
#==============================================================================
    plt.plot([t*1e3 for t in time],T, linewidth=3)
    plt.xlabel('time(ms)', FontSize=15)
    plt.ylabel('Temperature (K)', FontSize=15)
    plt.show()
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot([t*1e3 for t in time], X_fuel, color = 'red', linewidth = 3)
    ax1.set_xlabel('time (ms)', FontSize=15)
    ax1.set_ylabel('X_fuel', color='red', FontSize=15)
#    ax2.plot([t*1e3 for t in time], X_O, color = 'green', linewidth = 3)
    ax2.set_ylabel('X_CO2', color='green', FontSize=15)
    ax2.plot([t*1e3 for t in time], X_CO2, color = 'blue', linewidth = 3)
    plt.show()
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot([t*1e3 for t in time], m_dot, color = 'red', linewidth = 3)
    ax1.set_xlabel('time (ms)', FontSize=15)
    ax1.set_ylabel('m_dot', color='red', FontSize=15)
    ax2.plot([t*1e3 for t in time], q_dot, color = 'green', linewidth = 3)
    ax2.set_ylabel('q_dot', color='green', FontSize=15)
    plt.show()
    
    fig,ax1 = plt.subplots()
    ax1.plot([t*1e3 for t in time], num_d, linewidth=3)
    ax1.set_xlabel('time(ms)', FontSize=15)
    ax1.set_ylabel('#droplets', FontSize=15)
    plt.show()
   
execute()   