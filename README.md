# SILL_DEGASSING
source files for numerical simulations presented in:
"Rapid degassing in basaltic sills as a source of Deep Long Period volcanic earthquakes" by Oleg Melnik1,Vladimir Lyakhovsky, and Nikolai M. Shapiro

**fort_source.f90** - source file for the main program
Please note that the Fortan code was developed for IBM AIX 7.1 system and should be compiled using the 
IBM compiler XLF 15.1.2 using open multiprocessing (omp) option. Please contact Dr. Vladimir Lyakhovsky
 _vladimir.lyakhovsky@gmail.com_ 

**input_files.zip** - input files for the main code

**output_potency_seismograms** - elements of the mesh grid

**POTENCY.zip** - output files containing values of seismic potency tensor

**output_potency_seismograms.py** - python script to plot final results

**lgne.sac** - example of real data shown in the figures
