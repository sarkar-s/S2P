# Run gel modeling code
import sys, os
from optparse import OptionParser
from PolymerFluctuationSSA import PolymerFluctuationSSA
#from multiprocessing import Pool
import multiprocessing
from itertools import repeat
import platform as pl
import warnings
import numpy as np

def create_empty_files(casename):
    os.chdir("..")
    ofile = open(casename+'_inputs.txt','w')
    
    print >> ofile, '# This file should contain the directory name and the specific input file names for this job\n'
    
    print >> ofile, 'Case folder = '+casename
    print >> ofile, 'Monomers inputs = '+casename+'_monomers_inputs.txt'
    print >> ofile, 'Initiator inputs = '+casename+'_initiator_inputs.txt'
    print >> ofile, 'Curing inputs = '+casename+'_curing_inputs.txt'
    print >> ofile, 'Simulation inputs = '+casename+'_simulation_inputs.txt'
    
    ofile.close()
    
    if os.path.isdir(casename)==False:
        os.mkdir(casename)

    os.chdir(casename)
    
    # Create the sample monomer file
    monomers_file = open(casename+'_monomers_inputs.txt','w')
    
    print >> monomers_file, '# This file contains all the physical and chemical details of the photopolymerizing system\n'
    print >> monomers_file, '# Polymerization parameters'
    print >> monomers_file, 'Propagation rate constant = '
    print >> monomers_file, 'Cyclization(or cross-linking) factor = '
    print >> monomers_file, 'Termination rate constant = '
    print >> monomers_file, 'Mean neighboring functional groups = '
    print >> monomers_file, 'Functionality = '
    
    print >> monomers_file, '\n# Monomer properties'
    print >> monomers_file, 'Molecular volume (in Angstroms^3) = '
    print >> monomers_file, 'Molecular weight (in Da) = '
    print >> monomers_file, 'Viscosity (in Pa.s) = '
    print >> monomers_file, 'Reaction enthalpy (in kJ/mole) = '
    
    print >> monomers_file, '\n# Bond energy specification'
    print >> monomers_file, 'Bond length stiffness (in eV) = '
    print >> monomers_file, 'Bond angle stiffness (in eV) = '
    print >> monomers_file, 'Interaction energy with initiators (in eV) = '
    print >> monomers_file, 'Reaction enthalpy (in kJ/mole) = '
    print >> monomers_file, 'Density (in g/cm^3) = '
    print >> monomers_file, 'Speed of sound (in m/s) = '
    print >> monomers_file, 'Boiling point (in K) = '
    print >> monomers_file, 'Coefficient of thermal expansion (/K) = '
    print >> monomers_file, 'Infinite frequency shear modulus (Pa) = '
    
    monomers_file.close()
    
    # Create the sample initiator file
    initiators_file = open(casename+'_initiator_inputs.txt','w')
    
    print >> initiators_file, '# This file contains the specifics of the initiator/co-initiator system\n'
    
    print >> initiators_file, 'Initiation type (unimolecular or bimolecular) = '
    print >> initiators_file, 'Percentage of initiators (mole fraction) = '
    print >> initiators_file, '\n# The following inputs are only for a bimolecular initiation system'
    print >> initiators_file, 'Percentage of co-initiators (mole fraction) = '
    print >> initiators_file, 'Initiator radius = '
    print >> initiators_file, 'Co-initiator radius = '
    print >> initiators_file, 'Initiator molecular weight = '
    print >> initiators_file, 'Co-initiator molecular weight = '
    print >> initiators_file, 'Excitation lifetime of the initiator = '
    
    initiators_file.close()
    
    # Create the sample curing condition file
    curing_file = open(casename+'_curing_inputs.txt','w')
    
    print >> curing_file, '# This file contains the details of the curing light\n'
    
    print >> curing_file, 'Lighting mode (Linear ramp or Pulse) = '
    print >> curing_file, 'Excitation rate constant = '
    print >> curing_file, 'Curing time (in seconds) = '
    print >> curing_file, 'Temperature (in K) = '
    
    print >> curing_file, '\n# Additional input for linear ramp'
    print >> curing_file, 'Ramp time (in seconds) = '
    
    print >> curing_file, '\n# Additional input for pulse'
    print >> curing_file, 'Pulse width (in seconds) = '
    print >> curing_file, 'Pulse gap (in seconds) = '
    
    curing_file.close()

    # Create the sample algorithm input file
    algo_file = open(casename+'_simulation_inputs.txt','w')
    
    print >> algo_file, '# This file contains the inputs for the BFM lattice and the algorithm\n'
    
    print >> algo_file, 'Lattice half size = '
    print >> algo_file, 'Simulation time = '
    print >> algo_file, 'Output interval = '
    print >> algo_file, 'Output casename = '+casename
    print >> algo_file, 'Write video output = '
    print >> algo_file, 'Debye function values = '
    
    algo_file.close()
    
    sys.exit()
    
def process_options():
    usage = "usage: %prog [option1] arg1 [option2] arg2"
    parser = OptionParser(usage=usage)
    parser.add_option("-i", dest="inputfile", help="Input filename for a photopolymerization simulation.", default="")
    parser.add_option("-n", dest="sample", help="Number of stochastic trajectories to be simulated.", default="")
    parser.add_option("-v", dest="verbose", help="Verbosity of display outputs.", default="1")
    
    [options, args] = parser.parse_args()
    
    if len(options.inputfile) != 0:
        try:
            InputFile = open(options.inputfile, 'r')
        except IOError:
            print >> sys.stderr , "ERROR : Cannot open inputfile. Check inputfile name."
            
        InputLines = InputFile.readlines()
        
        removeSet = []
        for l in InputLines:
            if l[0] == '#' or l[0] == '\n':
                removeSet.append(l)
        
        for rem in removeSet:
            InputLines.remove(rem)
            
        if options.verbose=="1":
            verbosity = True
        else:
            verbosity = False

        return InputLines, int(options.sample), verbosity

    else:
        print >> sys.stderr, "\n Improper usage. Please enter ' python gelmaker.py -h ' for instructions."
        sys.stdout.flush()
        sys.exit()
        
Input_Data = []

def simulation_job(input_lines,proc_no,verbosity):
    gelMaker = PolymerFluctuationSSA(input_lines,proc_no,verbosity)
    gelMaker.initialize_structs()
    gelMaker.time_loop()

def f(x):
    return x*x

if __name__ == '__main__':
    if pl.system()=='Windows':
        warnings.filterwarnings("ignore", category=np.VisibleDeprecationWarning) 

    input_lines, pool_size, verbosity = process_options()
    Input_Data.append(input_lines)
    
    if pl.system()=='Windows':
        simulation_job(input_lines,0,verbosity)
    else:
        #pool = Pool(processes=pool_size)
        
        if input_lines:
            procs = []
            
            for i in range(pool_size):
                p = multiprocessing.Process(target=simulation_job,args=(Input_Data[0],i,verbosity))
                
                procs.append(p)
                p.start()
            
            for p in procs:
                p.join()
                
            #    pool.apply_async(simulation_job,args=(Input_Data[0],i,verbosity,))
            #    
            #pool.close()
            #pool.join()