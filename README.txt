vaporizer_temp fix mdot oscillation problem, by using correct logics. This file contains droplet vaporization prior to reaching wetbulb temperature. An dT/dt ode is solved. 

vaporize_simple contains simple algebraic equation for heating. No ODE is solved as not needed. Less accurate. 

Main is the execution file. Droplet contains droplet class. entrainer generates droplets. Other files are used for testing and verification. Make a file folder for the correct codes. 

Still working on this. 