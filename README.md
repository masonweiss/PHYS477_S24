# PHYS477_S24
final project for PHYS 477 - Physics of Nuclear Energy

This project explores the use of the Bateman Equations to model isotopes known as "Neutron Poisons," particularly Xe-135, in a sample nuclear reactor. The relationship between I-135 and its decay product Xe-135 is described using a system of ODEs, which, when approximated using the Forward Euler Method or computed analytically, can demonstrate the impact of isotopes with high absorption cross-sections on reactor operation.

The main.py file serves as an entry point and driver for the program, and references the constants.py file which determines constants not associated with plotting or the numerical method. 

The plots generated by the program include:
* I-135 and Xe-135 concentrations for the following operations:
  * Reactor start-up from 0% to 100% power (analytic solution)
  * Reactor shut-down from 100% to 0% power (analytic solution)
  * Reduction from 100% to 50% power (numerical solution)
  * Increase from 50% to 100% power (numerical solution)
 
* Poison reactivity plots for Xe-135 for the following operations:
  * Reactor shut-down from 100% to 0% power (analytic solution)
  * Reduction from 100% to 50% power (numerical solution)
  * Increase from 50% to 100% power (numerical solution)
 
The Bateman_Nuclear_Poisons.pdf file includes a comprehensive description of the methods used in the programs described above, and a description of the workflow used to estimate the constants. 


Last Updated: April 26, 2024
