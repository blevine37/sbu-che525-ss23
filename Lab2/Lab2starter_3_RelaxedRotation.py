"""
Lab2starter_3_RelaxedRotation.py

This code uses Psi4 to perform a relaxed scan of the
torsional PES of hydrogen peroxide at the Hartree 
Fock level.  By relaxed we mean that we will scan several 
dihedral angles, and all other internal coordinate (bond 
lengths, angles) are optimized at each angle.

This is a starter code for Lab2.

Origin CHE 525, S21 Problem Development
Author: Tom Allison, Ben Levine
"""
#%%
# Import modules =======================================================
import psi4
import numpy as np  

#%%
# Set up Psi4 ==========================================================

psi4.core.clean()
psi4.core.clean_options()
psi4.set_memory('4000 MB')  # Can make this much larger on Seawulf, each compute node has more than 100 GB RAM
psi4.set_num_threads(4)    # Can make this much larger on Seawulf, each compute node can support 28 threads.
                            # But it doesn't help much for small molecules...
psi4.core.set_output_file('Lab2starter_3_RelaxedRotation.dat', False) #this command sets psi4's output to a file. Comment this line out if you want to see the output on the terminal.

#%%
# Define QC method to use and angles to scan ===========================

method = 'scf/6-31G'

#%%
# Set up Z-matrix string ==============
# This geometry was taken from the previous geometry optimization.
# Specifically, we opened the previous molden file in molden, opened
# molden's Z-matrix editor, and read the values.
HOOH_string= """
    o
 o   1 oo2     
 h    1 ho3         2 hoo3      
 h    2 ho4         1 hoo4          3 dih4   
 
oo2=        1.462461
ho3=        0.954335
hoo3=       101.178
ho4=        0.954364
hoo4=       101.183
dih4=       {0}"""

# Here is a list of angles that we want to compute, in degrees
angles = [180, 165, 150, 135, 120, 105, 90, 75, 60, 45, 30, 15, 0]

# This creates an empty 2-D array (matrix) that will store the energies 
# that come from our calculations.  The argument (0,2) tells empty that 
# we would request a 2 x 0 array, meaning that each row will have 2 
# elements, but initially there are zero rows.  We will add rows below.
Etheta = np.empty((0,2))

# This for loop will repeat for every element of the array angles, 
# defined on the previous line.  For each iteration, the variable theta
# will be set to a different angle (first 180, then 165, etc...)
for theta in angles:

    # Define the molecules.  Here the format function replaces 
    # placeholder, {0}, with the value of theta
    HOOH = psi4.geometry(HOOH_string.format(theta)) 
    
    # By default, many electronic structure codes will take advantage of
    # the symmetry of the molecules to speed up the calculations.  However,
    # as we scan the PES, the symmetry of HOOH will change.  So, here we 
    # turn off symmetry by setting the point group to C1
    HOOH.reset_point_group('c1')
 
    # This tells the geometry optimizer that we want to fix the dihedral 
    # between atoms 4, 2, 1, and 3.  (Keep in mind that the order you 
    # list the atoms in matters!  The + operator and str() function 
    # are useful tools to manipulate strings in python.  See how we define
    # molden file names below (or just google) for more details.
    psi4.set_module_options('optking', {'frozen_dihedral': '4 2 1 3 '})

    # Run a calculation of the energy, optimizing the energy using all the 
    #coordinates except the previously defined frozen dihedral.
    Ei, wfni = psi4.optimize(method, molecule = HOOH, return_wfn = True)

    # This adds the angle (theta) and energy (E) to the array (Etheta)
    # that we created above.  Each row will contain 2 values: theta and Ei
    Etheta = np.append(Etheta, [[theta,Ei]], axis=0)

    # The object wfn0 now contains the optimized geometry and 
    # associated wave function.  You can output this information 
    # to a molden file as follows.  Here str converts theta from a 
    # number to a string, and the + sign concatenates several
    # strings into one long sting.  For example, if theta is 180, the 
    # filename will be "Lab2starter_3_RelaxedRotation.180.molden").  This 
    # allows us to save a new molden file for each iteration of the loop.
    psi4.molden(wfni,"Lab2starter_3_RelaxedRotation."+str(theta)+".molden")

# Here we print the results array (Etheta) to a CSV file that can be
# read later by a python plotting program or the spreadsheet of your
# choice.  CSV files are also human readable.
np.savetxt("Lab2starter_3_RelaxedRotation.csv", Etheta, delimiter=",",header='H-O-O-H Dihedral Angle,E_relaxed',comments='')
