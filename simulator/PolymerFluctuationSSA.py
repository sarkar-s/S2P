""" A framework for simulating a cross-linked polymer network created during photopolymerization
..moduleauthor:: Swarnavo Sarkar <swarnavo.sarkar@nist.gov>
    ..platform: Unix, Windows
    ..synopsis: This code is built upon 4 main techniques from polymer and chemical kinetics theory.
    (1) Bond Fluctuation Model (BFM) to represent a polymer network.
    (2) Stochastic Simulation Algorithm to sample reaction and fluctuation events.
    (3) Lattice diffusion theory to fluctuate the polymer network.
    (4) Bond conduction theory to compute diffusion coefficient.
    (5) Coherent Potential Approximation to compute phonon properties.
"""

import os, sys
import numpy as np
#import scipy.signal as signal
from scipy.optimize import curve_fit
import random as rand
import math
from LightSource import *
from EffectiveDiffusion import *
from ThermalMotion import *
from NetworkDisorder import *
from Monomers import *
import IndexUtils
import vtk_writer
import string
import time

class PolymerFluctuationSSA(object):
    def __init__(self,master_input_lines,proc_no,_verbosity):
	"""Constructor for the PolymerFluctuationSSA instance
	This function reads all the chemical and simulation inputs
	and the number of samples to simulate and the verbosity of display.
	"""
	master_input_dict = {}
	self.verbosity = _verbosity
	
	for l in master_input_lines:
	    if string.find(l,' = ') != -1:
		thisSet = l.split(' = ')
		thisSet[0] = thisSet[0]
		thisSet[1] = thisSet[1].rstrip('\r\n')
		master_input_dict[thisSet[0]] = thisSet[1]
	
	self.case_directory = master_input_dict['Case folder']
	
	self.code_directory = os.getcwd()
	
	os.chdir('../'+self.case_directory)
	
	self.main_input_dict = {}
	
	monomers_file = master_input_dict['Monomers inputs']
	data_lines = self.read_input_lines(monomers_file)
	
	initiator_file = master_input_dict['Initiator inputs']
	data_lines += self.read_input_lines(initiator_file)
	
	curing_file = master_input_dict['Curing inputs']
	data_lines += self.read_input_lines(curing_file)
	
	simulation_file = master_input_dict['Simulation inputs']
	data_lines += self.read_input_lines(simulation_file)
	
	# Create input dictionary
	input_dict = {}
	
	for l in data_lines:
	    if string.find(l,' = ') != -1:
		this_set = l.split(' = ')
		this_set[0] = this_set[0]
		this_set[1] = this_set[1].rstrip('\r\n')
		input_dict[this_set[0]] = this_set[1]
		
	self.read_all_inputs(input_dict,proc_no)

    def read_input_lines(self,filename):
	"""This function reads all the lines in a file.
	"""
	file_handle = open(filename, 'r')
	file_lines = file_handle.readlines()
	file_handle.close()
	
	return file_lines

    def read_all_inputs(self,input_dict,proc_no):
	""" This function reads the lines from the input file and creates the dictionary of input quantities.
	Args:
	input_lines (list): The list of input data
	"""
	self.proc_no = proc_no

	# Monomer inputs
	self.Nm0 = float(input_dict['Mean neighboring functional groups'])
	self.temp = float(input_dict['Temperature (in K)'])
	self.eta = float(input_dict['Viscosity (in Pa.s)'])
	self.bond_angle_stiffness = float(input_dict['Bond angle stiffness (in eV)'])
	self.bond_length_stiffness = float(input_dict['Bond length stiffness (in eV)'])
	self.interaction_energy = float(input_dict['Interaction energy with initiators (in eV)'])
	self.reaction_enthalpy = float(input_dict['Reaction enthalpy (in kJ/mole)'])
	self.density = float(input_dict['Density (in g/cm^3)'])
	
	try:
	    self.max_length = int(input_dict['Chain length'])
	except KeyError:
	    self.max_length = None

	#self.thermal_alpha = float(input_dict['Coefficient of thermal expansion (/K)'])
	
	# Simulation inputs
	self.lattice_size_half = float(input_dict['Lattice half size'])
	self.simulation_time = float(input_dict['Simulation time'])
	self.writeInt = float(input_dict['Output interval'])
	self.outName = input_dict['Output casename']
	self.write_video = input_dict['Write video output']
	self.equilibrium_tolerance = float(input_dict['Equilibrium tolerance'])
	self.minimal_sample_size = float(input_dict['Minimum sample size'])
	
	try:
	    self.step_limit = int(float(input_dict['Termination step']))
	except KeyError:
	    self.step_limit = -1
	
	# Initiator inputs
	self.nI_pt = float(input_dict['Percentage of initiators (mole fraction)'])
	self.curing_time = float(input_dict['Curing time (in seconds)'])
	self.initiation_type = input_dict['Initiation type (unimolecular or bimolecular)']
	# Interaction energy between the monomer units
        #self.eps_self = float(input_dict['Self interaction energy (in eV)'])

	self.total_bonds = 0

	if self.initiation_type=='bimolecular':
	    self.nI_co_pt = float(input_dict['Percentage of co-initiators (mole fraction)'])
	    
	    # Interaction energy between initiators and monomers
	    self.eps_I = float(input_dict['Interaction energy with initiators (in eV)'])
	    
	    # Interaction energy between co-initiators and monomers
	    self.eps_co_I = float(input_dict['Interaction energy with initiators (in eV)'])
	    
	self.all_monomers = Monomers(input_dict)
	self.lightObject = LightSource(input_dict)
	self.diffusion = EffectiveDiffusion(input_dict)
	#self.diffusion.set_energy_levels(self.bond_length_stiffness)
	
	# Object to calculate specific heat capacity
	self.directionBonds = np.zeros(shape=(3))
	#self.MaterialModel = MaxwellModel(input_dict)
	
	self.create_output_directory()
	
	#self.compute_shrinkage_factor()
	
    def create_output_directory(self):
	#Change working directory
	self.oldcwd = os.getcwd()
	os.chdir("..")
	tempcwd = os.getcwd()
	os.system("mkdir "+self.case_directory)
	os.chdir(tempcwd+"/"+self.case_directory)
	os.system("mkdir "+self.outName)
	os.chdir(self.outName)
	
	sys.stdout = open(self.outName+'_output_'+str(self.proc_no)+'.txt','w')
	sys.stderr = open(self.outName+'_error_'+str(self.proc_no)+'.txt','w')
	
    def initialize_structs(self):
	""" Initializes the data structures for the simulation
	"""
	# kB - Boltzmann factor
        self.kB = 1.38065E-23

	# Full dimension of the BFM lattice
	self.lattice_size = 2*int(self.lattice_size_half)
	
	# Limiting lattice distance for the pair correlation function
	self.R_lim = 10 #int(self.lattice_size/8)
	
	IndexUtils.make_interaction_set(self.all_monomers.interaction_range)

	# Lattice size for the vertices of the BFM lattice
	self.super_lattice_size = self.lattice_size
	
	# Total sites in the lattice
	self.total_lattice_sites = self.lattice_size**3
	
	# Initialize the upper limit in the IndexUtils module
	IndexUtils.get_max_index(self.lattice_size)
	self.super_lattice = np.zeros(shape=(self.super_lattice_size,self.super_lattice_size,self.super_lattice_size))
	self.junction_super_lattice = np.zeros(shape=(self.super_lattice_size,self.super_lattice_size,self.super_lattice_size))

	# Matrices to note which chain is taking up a lattice location
	self.location_chain = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	self.location_boxno = -1*np.ones(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Matrix to note which chain is occupying a junction point
	self.junction_chain = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	self.junction_boxno = -1*np.ones(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Matrix to note which chain is occupying a junction point
	self.junction_neighbor_sites = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Matrix to note secondary junction neighbor
	self.junction_second_neighbors = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Matrix to store the type of monomer at the location
	self.location_monomers = -1*np.ones(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Matrix to store the type of functionality at the location
	self.location_functionalities = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# Diffusion tracking objects
	#self.blocked_matrix = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	self.neighbor_sites = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	self.second_neighbor_matrix = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	self.total_sites = int(self.lattice_size)**3
	
	# Dict for counting junctions in each type of chains
	self.junctions_in_chains = {}

        # Number of available sites for the a node to move around (ranges from 1 to 6)
        #self.movability = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	# This dict has a key which refers to a lattice site and the value is a list of size 6 entries 0 or 1
        self.movable_nodes = []
        self.total_movable_nodes = 0
	self.highly_movable_nodes = []
        self.total_highly_movable_nodes = 0
	self.highly_movable_propensities = 0.0
	self.high_propensities = []

	# Moving probability values
	self.movable_dict = np.zeros(shape=(self.lattice_size**3,6))
	self.moving_probs = np.zeros(shape=(self.lattice_size**3,6))
	self.bond_motion_update = np.zeros(shape=(self.lattice_size**3,6))
	self.max_moving_probs = np.zeros(shape=(self.lattice_size**3))

	# Count of interfaces available for a molecule to move around
	self.total_diffusible_interfaces = 6*(self.lattice_size**3)
	self.transient_diffusible_interfaces = self.total_diffusible_interfaces
	self.max_interface_update = 3*(3**2)*(3+1)

	# Energy lattice - stores the energy associate with a given node in the polymer network
	#self.associated_energy = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	
	# Total energy of the network
	self.total_network_energy = 0.0

	# All thermal propensities
	self.all_thermal_propensities = {}
	self.thermal_propensities = np.zeros(shape=(self.lattice_size**3))
	
	self.mean_diffusion_factor = 1.0
	self.d_factor_oligomer = 1.0
	self.oligomer_eps_self = self.all_monomers.mean_eps_self
	
	# Convergence parameters
	self.counted_steps = 0
	self.current_mean_fluctuation = 0.0
	self.last_mean_fluctuation = 0.0
	self.total_fluctuation = 0.0
	self.fluctuation_values = []

	# Counter to keep track of the different type of events in the system
	self.event_counter = {}
	self.event_counter['polymerize'] = 0
	self.event_counter['fluctuate'] = 0
	self.event_counter['initiate'] = 0
	self.event_counter['radical termination'] = 0
	
	if self.initiation_type=='bimolecular':
	    self.event_counter['excite'] = 0
	    self.event_counter['decay'] = 0
	    
	# Statistics of leaping
	self.fluctuation_epr = 0

	# Radical neighborhood data
	self.problist = {}
	self.totalProbSet = {}
	self.neighbor_monomers = {}
	self.neighbor_groups = {}
	self.max_A = {}
	self.fluctuation_quanta = {}
	self.max_fluc_trials = 0
	
	for i in xrange(0,7):
	    self.fluctuation_quanta[i] = 0
	
	# Sites that can form links to the radicals
	self.non_zero_sites = {}
	
	self.junction_list = []

        # Lattice box that went through the latest update
        self.last_update_box = None

	# Step number for the last polymerization event
	self.last_non_fluctuation_step = 0
	self.last_non_fluctuation_time = 0.0

	# Stores the last 10 values of configurational entropy
	self.entropy_list = []
	self.max_entropy_rate = 0.0

        # Terminated pairs
        self.terminated_pairs = {}
	self.early_terminations = []
	
	# Connectivity among chains
	self.chain_connectivity = {}

        # List of radical chains hanging around without nearby monomers
        self.no_neighbors = []

        # Degree of conversion for the polymerization process
        self.DoC = 0.0
        self.oldDoCfact = 0
                
        # Functional radical matrix
        # Stores the end point indices of the active macromolecule radicals
        self.Nr = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
        
	# Total number of sites occupied in the lattice
	self.blocked_site_count = 0
	# Total number of sites neighboring to a polymer site
	#self.neighbor_site_count = 0
	# Total number of sites which are a surface site
	self.surface_site_count = 0
	# Total number of junction surface
	self.junction_surface_count = 0
	self.critical_phi_0 = 0.0

	# Total no. of nodes in the system
	self.total_nodes = 0
	self.last_node_size = 0
	self.totalJunctionPoints = 0
	self.realJunctionPoints = 0
	self.max_connected_unit = 0
	self.last_unit = 0
	
	# Energy population counter
	self.energy_level_population = np.zeros(shape=(13))
	self.energy_counter = {}
	self.bond_populations = {}
	self.bond_populations['high'] = 0
	self.bond_populations['low'] = 0
	self.total_energy_histogram = {}

	self.last_polymer_propensity = 0.0
	
	# Polymer chains - dictionary of coordinates making up the polymer chains
	self.all_chains = {}
	self.active_chains = []
	# Dictionary of monomer types making up the polymer chains
	self.chain_monomers = {}
	
	self.all_monomers.compute_domain_volume(self.density,self.lattice_size_half)

	self.unit_vol, self.sysVol, self.total_effective_units, self.molWt = self.all_monomers.get_volumes()
	
	# Create initiators
	self.create_initiators()
	self.all_monomers.compute_monomer_populations(self.density,self.lattice_size_half,self.nI)
	
	self.TotalMonomers = self.all_monomers.total_monomers       
        self.initialFuncs = self.all_monomers.total_functionality
        self.transientFuncs = self.initialFuncs
        self.transientMonomers = self.TotalMonomers

	# Real length of a lattice box
	self.lattice_h = (self.sysVol**(1.0/3.0))/self.lattice_size
	
	# Heat data
	self.reaction_enthalpy *= (1000.0/(6.022E+23))
	self.heat_production = 0.0
	self.heat_increment = 0.0
	
	# Current time
	self.current_time = 0

	# Initialize Thermal Motion object
	self.fluctuator = ThermalMotion(self.temp,self.sysVol,self.bond_length_stiffness,self.bond_angle_stiffness,self.lattice_size,self.lattice_h)
	self.fluctuator.compute_bond_angle_energies(IndexUtils.all_angles)
	self.fluctuator.compute_monomer_self_diffusion(self.unit_vol,self.all_monomers.molecular_vol,self.eta)
	self.fluctuator.get_monomer_self_diffusion(self.all_monomers.hopping_rate)
	
	self.fluctuator.compute_mean_diffusivity(self.all_monomers.monomer_fraction)
	
	# Initialize phonon properties object
	self.disorder_calculator = NetworkDisorder()
	#self.disorder_calculator.setup_scaling_factors(self.eps_self,self.lattice_h)
	
	# Temperature variables
	self.temperature_update = 0.0
	self.current_temperature = self.temp
	
	# Dict to store the population of different types of bond lengths
	self.bond_pops = {}
	self.bond_pops[4.0], self.bond_pops[5.0], self.bond_pops[6.0] = 0, 0, 0
	self.bond_pops[9.0], self.bond_pops[10.0] = 0, 0
	
	self.bond_pop_update = {}
	
	# Fluctuation trial counter
	self.fluctuation_trials, self.trial_updates = {}, {}
	self.fluctuation_trials['success'] = 0
	self.fluctuation_trials['failure'] = 0

	self.trial_updates['success'] = 0
	self.trial_updates['failure'] = 0
	
	# Quantities to keep track of equilibriation
	self.fluctuation_samples = []
	#self.fluctuation_times = []
	self.fluctuation_sum = 0.0
	self.trajectory_entropy = 0.0
	self.trajectory_samples = 0
	self.fluctuation_sample_size = 0
	self.fluctuation_square_sum = 0.0
	self.last_std_fluc = 0.0
	self.last_mean_fluc = 0.0
	self.last_polymer_step = 0
	self.max_fluc_wt = 6.0
	self.fluctuation_transition_status = False
	
	# Create correlation radius list
	self.create_correlation_radius_list()

        #self.modObject = modulus(self.temp,self.sysVol,4)
        #self.computeNm()

	# Dipole vector
	self.dipole_vector = [0,0,0]
	
	# Global propensity dict
	self.global_propensities = {}
	self.global_propensities['initiate'] = 0.0
	self.global_propensities['fluctuate'] = 0.0
	self.global_propensities['polymerize'] = 0.0
	
	if self.initiation_type=='bimolecular':
	    self.global_propensities['excite'] = 0
	    self.global_propensities['decay'] = 0
	
	self.non_fluctuation_propensity = 0.0

	# Correction factor for fluctuation from failed moves
	self.fluctuation_rate_factor = 1.0

	# Effective diffusion factor for the bimolecular initiation process
	self.diffusion_factor = 1.0
	
	# Entropy of the current 
	self.current_entropy = 0.0
	
	# Counter for the video file images
	self.pixelCount = 0
	# Counter for vtk files
	self.vtk_counter = 0
	
	# Initial effective diffusion factor
	self.D_em = 1.0

	#Change working directory
	#self.oldcwd = os.getcwd()
	#os.chdir("..")
	#tempcwd = os.getcwd()
	#os.system("mkdir "+self.case_directory)
	#os.chdir(tempcwd+"/"+self.case_directory)
	#os.system("mkdir "+self.outName)
	#os.chdir(self.outName)
	#
	self.create_output_files()
	#self.create_ordered_lattice()
	#self.write_vtk_files()
	#self.compute_ordered_correlation()
	#sys.exit()

    def compute_shrinkage_factor(self):
	unaffected_volume = self.MolVol - self.funcNum*40.75
	
	self.shrinkage_factor = 1.0 - (unaffected_volume + 0.9*(self.funcNum*40.75))/self.MolVol
		
    def create_initiators(self):
	# Initial population of initiators
	self.nI = int(0.01*self.nI_pt*self.total_effective_units)
	self.max_monomers = self.total_effective_units - self.nI
	
        # Initial population of excited initiators
        self.nI_excited = 0

	# Initial count of amines
	if self.initiation_type=='bimolecular':
	    self.n_co_I = int(0.01*self.nI_co_pt*(self.lattice_size_half**3))
	    self.lightObject.set_population(self.nI,self.n_co_I,self.sysVol)
	elif self.initiation_type=='unimolecular':
	    self.lightObject.set_population(self.nI,None,None)

    def create_output_files(self):
	# Create output file
	self.DCoutfile = open('Time_DC'+str(self.proc_no)+'.csv','w')
	print >> self.DCoutfile, 'Event No,Time,DC,cl'
	self.DCoutfile.close()
	
	self.Chainoutfile = open('chainstat'+str(self.proc_no)+'.csv','w')
	print >> self.Chainoutfile, 'Total chains,Active chains,Total Nodes,Junction Points'
	self.Chainoutfile.close()
	
	self.polymerfile = open('polydata'+str(self.proc_no)+'.csv','w')
	self.polymerfile.close()
	
	self.event_file = open('eventstat'+str(self.proc_no)+'.csv','w')
	out_string = ''
	for event_type in self.event_counter.keys():
	    out_string += ','+event_type

	print >> self.event_file, 'Time,DoC'+out_string
	self.event_file.close()
	
	self.fluctuation_file = open('fluctuation'+str(self.proc_no)+'.csv','w')
	self.fluctuation_file.close()
	
	self.site_file = open('site_fractions'+str(self.proc_no)+'.csv','w')
	print >> self.site_file, 'Time,DoC,Empty sites,Surface sites,Blocked sites'
	self.site_file.close()
	
	self.g_bar_file = open('g_bars'+str(self.proc_no)+'.csv','w')
	print >> self.g_bar_file, 'Time,DoC'
	self.g_bar_file.close()
	
	self.diffusivity_file = open('diffusivity'+str(self.proc_no)+'.csv','w')
	print >> self.diffusivity_file, 'Time,DoC,Network Size,Molecule Diffusivity,Network Diffusivity'
	self.diffusivity_file.close()
	
	self.propensity_file = open('propensities'+str(self.proc_no)+'.csv','w')
	out_string = ''
	for event_type in self.global_propensities.keys():
	    out_string += ','+event_type

	print >> self.propensity_file, 'Time,DoC'+out_string
	self.propensity_file.close()
	
	self.energy_file = open('energy'+str(self.proc_no)+'.csv','w')
	print >> self.energy_file, 'Time,DoC,Network Size,Average energy'
	self.energy_file.close()
	
	self.bond_length_file = open('bond_lengths'+str(self.proc_no)+'.csv','w')
	print >> self.bond_length_file, 'Time,DoC,2.0,'+str(math.sqrt(5))+','+str(math.sqrt(6))+',3.0,'+str(math.sqrt(10.0))
	self.bond_length_file.close()
	
	self.heat_file = open('heat_info'+str(self.proc_no)+'.csv','w')
	print >> self.heat_file, 'Time,DoC,Heat Production,Temperature'
	self.heat_file.close()
	
	self.modulus_file = open('modulus'+str(self.proc_no)+'.csv','w')
	print >> self.modulus_file, 'Time,DoC,Modulus,Relaxation time'
	self.modulus_file.close()
	
	self.leap_file = open('fluctuation_epr'+str(self.proc_no)+'.csv','w')
	print >> self.leap_file, 'Time,DoC,Leaps'
	self.leap_file.close()
	
	self.gel_file = open('gel_transient'+str(self.proc_no)+'.csv','w')
	print >> self.gel_file, 'Time,DoC,phi_0,phi_c'
	self.gel_file.close()
	
	# Write vtk template
	vtk_writer.set_case_name(self.outName)
	vtk_writer.create_vtk_template(2*int(self.lattice_size_half))

	self.trial_file = open('fluctuation_trials'+str(self.proc_no)+'.csv','w')
	print >> self.trial_file, 'Time,DoC,Success,Failure'
	self.trial_file.close()
	
	self.monomer_file = open('monomer_count'+str(self.proc_no)+'.csv','w')
	print >> self.monomer_file, 'Time,DoC,Remaining Monomers'
	self.monomer_file.close()
	
	self.movable_file = open('movable_nodes'+str(self.proc_no)+'.csv','w')
	print >> self.movable_file, 'Time,DoC,Loose nodes,Stiff nodes,Blocked nodes'
	self.movable_file.close()
	
	self.component_file = open('component_size'+str(self.proc_no)+'.csv','a')
	print >> self.component_file, 'Time,DoC,Component Size1,Component Size2'
	self.component_file.close()
	
	
    def time_loop(self):
	""" This function performs
	"""
	self.start_time = time.time()
	
	self.step_no = 0
	self.relDoC = 1.0
	self.A_total = 1000.0
	leap_counter = 0
	leap_condition = 'no leap'
	step_type = None
	
	self.prop_count = 0
	
	self.global_propensities['fluctuate'] = 0.0

        while self.check_for_termination() == 0:
	#    if self.total_nodes>=1000:
	#	for chain in self.all_chains:
	#	    for boxno in xrange(0,len(self.all_chains[chain])):
	#		print self.get_connections(chain,boxno)
	#	sys.exit()
	    
	    self.create_polymer_propensity_list()
	    self.global_propensities['polymerize'] = self.A_polymer

	    #self.global_propensities['fluctuate'] =
	    self.get_network_diffusivity()
	    
	    if self.initiation_type == 'bimolecular' and self.diffusion_factor>0.0:
		self.global_propensities['excite'] = self.lightObject.get_excitation_propensity(self.current_time)
		self.global_propensities['decay'] = self.lightObject.get_decay_propensity()
		self.lightObject.compute_total_diffusivity(self.d_factor_initiator,self.d_factor_coinitiator)
		self.global_propensities['initiate'] = self.lightObject.get_initiation_propensity()

		if self.diffusion_factor<=0.0:
		    self.global_propensities['excite'] = 0.0
		    self.global_propensities['decay'] = 0.0
		    self.global_propensities['initiate'] = 0.0
	    elif self.initiation_type=='unimolecular':
		self.global_propensities['initiate'] = self.lightObject.get_excitation_propensity(self.current_time)

	    self.get_radical_termination_propensity()
	    
	    if self.step_no>1 and step_type=='fluctuate' and self.fluctuation_transition_status==False:
		this_value = self.fluctuation_samples[-1]
		lidx = self.get_lidx_from_cubidx(self.last_update_box)
		#reverse_prob = self.moving_probs[lidx][self.reverse_move]/sum(self.moving_probs[lidx])
		reverse_prob = self.moving_probs[lidx][self.reverse_move]/self.global_propensities['fluctuate']
		
		try:
		    self.fluctuation_samples[-1] = math.log(this_value/reverse_prob)
		except ZeroDivisionError:
		    print self.moving_probs[lidx], self.reverse_move, lidx
		    new_pos = self.get_new_position(self.last_update_box,IndexUtils.thermalMoves[self.reverse_move])
		    
		    print self.check_super_lattice_thermal(self.last_update_box,new_pos), self.reverse_move
		    sys.exit()
		#self.fluctuation_times.append(self.current_time)
		
		self.fluctuation_sum += self.fluctuation_samples[-1]
		self.trajectory_entropy += self.fluctuation_samples[-1]
		self.trajectory_samples += 1
		#print self.fluctuation_sum, self.fluctuation_samples[-1]
		self.fluctuation_sample_size += 1

	    if self.fluctuation_transition_status==False:
		leap_condition = self.check_for_equilibrium_leap()
		#self.fluctuation_statistics()

	    step_type, step_size, event_prob, leap_condition = self.get_event_type_and_step(leap_condition)
	    self.event_counter[step_type] += 1
	    self.step_no += 1
	    
	    if self.lightObject.lightmode=='Instant' and self.global_propensities['initiate']>0:
		self.step_no -= 1
		step_size = 0
		step_type = 'initiate'

	    if step_type!='fluctuate':
		#self.entropy_list = []
		self.fluctuation_samples = []
		#self.fluctuation_times = []
		self.fluctuation_sum = 0.0
		self.fluctuation_sample_size = 0
		#self.fluctuation_square_sum = 0.0
		#self.max_entropy_rate = 0.0
		self.last_non_fluctuation_step = self.step_no
		self.last_non_fluctuation_time = self.current_time

	    if self.step_no%self.writeInt==0 and self.total_nodes>0:
		self.compute_max_connected_unit()
		self.write_intermediate_output(step_type)

            if step_type == 'excite':
                self.lightObject.update_excitation()
            elif step_type == 'decay':
                self.lightObject.update_decay()
            elif step_type == 'initiate':
                self.lightObject.update_population()
                self.start_chain()
		self.last_polymer_step = self.step_no
            elif step_type == 'polymerize':
                # Select which chain will react with which neighbor
                chain, idx, monomer, react = self.get_chain_neighbor_reaction()
		#self.fluctuation_samples = []
		#self.fluctuation_times = []
		self.last_polymer_step = self.step_no
		#print 'Polymerization type: ', react, chain, self.all_chains[chain][-1]

		if react=='prop':
		    self.prop_count += 1
		    self.all_monomers.update_monomer_population(monomer)

                # Propagate chain if reaction is propagation
                if react=='prop' or react=='cross':
		    self.all_chains[chain].append(idx)
		    self.chain_monomers[chain].append(monomer)
		    #self.min_fluc_propensity = self.approximate_minimum_fluctuation_propensity()
		    self.heat_production += self.reaction_enthalpy
		    self.heat_increment = self.reaction_enthalpy
		    
		    if self.all_chains[chain][-1]==self.all_chains[chain][-2]:
			print 'Node extension ERROR'
			sys.stdout.flush()

		    if react=='prop':
			self.total_bonds += 1
			
		    if IndexUtils.get_distance2(self.all_chains[chain][-1],self.all_chains[chain][-2])==9.0:
			self.bond_populations['low'] += 1
		    else:
			self.bond_populations['high'] += 1

                elif react == 'term':
                    self.all_chains[chain].append(idx)
                    term_chain = self.get_chain_to_terminate(chain)
                    
                    # Now move these two chains from active to closed ones
                    if term_chain != None:
			self.close_chains(chain,term_chain)
			#self.min_fluc_propensity = self.approximate_minimum_fluctuation_propensity()
                    else: # Single chain termination
			print 'zero chain : ', self.totalProbSet[chain]
			sys.stdout.flush()
			sys.exit()
            elif step_type == 'fluctuate':
		self.fluctuation_status = self.take_fluctuation_step(event_prob)
		
		#print 'Fluctuation status: ', self.fluctuation_status
		self.fluctuation_trials[self.fluctuation_status] += 1
	    elif step_type == 'radical termination':
		self.terminate_a_pair()
		
	#    if step_type!=None and step_type!='excite':
	#	self.effected_all = self.influenced_nodes(self.last_update_box)
	#	self.get_thermal_propensity()
		    
	    if step_type in ['polymerize','initiate','fluctuate']:
		#print step_type, self.all_chains[1]
		self.effected_all = self.influenced_nodes(self.last_update_box,step_type)
		self.get_thermal_propensity()
		self.update_diffusivities()
		self.get_total_diffusion_factor()

		self.check_radicals()

		#print 'Step size: ', step_size
		self.current_time += step_size
		self.compute_doc(step_type)
		#self.check_all_bond_lengths()

            #self.get_polarization_density()            
            empty_sites = self.lattice_size_half**3 - self.total_nodes

	    #print 'Nodes: ', sum(sum(sum(self.location_chain))), self.total_nodes
	    #sys.stdout.flush()
	    
	#    if self.total_nodes!=sum(sum(sum(self.location_chain))):
	#	print 'loc chain error', self.total_nodes, sum(sum(sum(self.location_chain)))
	#	sys.exit()
	
	    try:
		self.total_energy_histogram[self.bond_populations['high']] += 1
	    except KeyError:
		self.total_energy_histogram[self.bond_populations['high']] = 1
	
	    if self.verbosity==True or self.step_no%self.writeInt==0:
		self.all_stdouts(step_type,step_size)


	    #if self.current_time>self.curing_time and self.curing_time>(self.current_time-step_size):
		#self.compute_crosslinks_prob('light_off')
		#self.compute_MC_pair_correlation('light_off')

	print 'Bond populations: ', self.bond_populations
	
	self.write_energy_histogram()
	
	self.all_stdouts(step_type,step_size)
	self.write_intermediate_output(step_type)
	#self.compute_free_volume()
	if self.proc_no==0:
	    self.write_final_vtk_files()
	self.write_chains()
		
	self.check_total_effective_units()
	self.compute_mol_wt_dist()
	self.compute_hopping_probs()
	self.compute_phonon_distribution()
	#self.compute_crosslinks_prob('final')
	self.compute_MC_pair_correlation('final')
	self.write_summary()

        #if self.WriteCoords:
        #    self.writeChains()

        os.chdir(self.oldcwd)

    def compute_cg(self):
	coord = [0,0,0]
	
	for idx in self.all_chains[1]:
	    coord[0] += idx[0]
	    coord[1] += idx[1]
	    coord[2] += idx[2]

	coord[0] *= 1.0/float(len(self.all_chains))
	coord[1] *= 1.0/float(len(self.all_chains))
	coord[2] *= 1.0/float(len(self.all_chains))
	
	ofile = open('chain_cg'+str(self.proc_no)+'.csv','w')
	
	print >> ofile, self.current_time,',',coord[0],',',coord[1],',',coord[2]

	ofile.close()	

    def check_for_termination(self):
        check = 0
	self.termination_type = None

	if self.A_total<(1e-12):
	    self.termination_type = 'no propensity'
	    check = 1

        #if self.simulation_time != 'infinity':
	if self.current_time >= self.simulation_time:
	    self.termination_type = 'time over'
	    check = 1

	#if len(self.active_chains)==0 and len(self.all_chains)==self.nI:
	#    self.termination_type = 'nothing active'
	#    check = 1

        if self.transientFuncs == 0:
	    self.termination_type = 'no groups'
	    check = 1
	    
	if self.step_no==self.step_limit:
	    self.termination_type = 'step limit'
	    check = 1
	    
	#if self.non_fluctuation_propensity==0.0 and self.curing_time!=0.0:
	#    status = False #self.check_for_convergence()
	#    if status==True:
	#	self.termination_type = 'Convergence'
	#	check = 1
		
	if check==1:
	    print 'Proper termination'
	    print 'Termination type: ', self.termination_type
	    sys.stdout.flush()

        return check
    
    def check_for_convergence(self):
	exit_status = False
	
	if self.counted_steps<1000:
	    self.fluctuation_values.append(self.global_propensities['fluctuate'])
	    #self.total_fluctuation += self.global_propensities['fluctuate']
	    self.counted_steps += 1
	else:
	    self.current_mean_fluctuation = sum(self.fluctuation_values)/1000.0
	    self.std_fluctuation = np.std(self.fluctuation_values)
	    tolerance = 1.0/math.sqrt(self.total_nodes)
	    self.counted_steps = 0
	    self.fluctuation_values = []
	    
	    if abs(self.std_fluctuation/self.current_mean_fluctuation)<(1e-3):
		exit_status = True
	    else:
		print 'Mean fluctuations: ', abs(self.std_fluctuation/self.current_mean_fluctuation)
	
	return exit_status

    def start_chain(self):
        """ This function picks a random site in the lattice, if it is empty
	and initiates a new macromolecular chain from this site.
	"""
	
	check = 0
	while check == 0:  
            boxIdx1 = rand.randint(0,self.lattice_size-1)
            boxIdx2 = rand.randint(0,self.lattice_size-1)
            boxIdx3 = rand.randint(0,self.lattice_size-1)
	    
            idxs = [boxIdx1,boxIdx2,boxIdx3]
	    
            if self.check_if_empty(idxs) == 'Yes':
                """ Check the current number of polymer chains
		create a new list for this chain
		and add this chain to the active chain dict
		"""
		plen = len(self.all_chains)+1
		self.all_chains[plen] = [idxs]
		self.active_chains.append(plen)
		self.chain_connectivity[plen] = []
		self.chain_monomers[plen] = [-1]
		#self.junctions_in_chains[plen] = 0

                """ Mark the location of the free radical in the lattice
		"""
		self.Nr[idxs[0],idxs[1],idxs[2]] += 1
                self.location_chain[idxs[0],idxs[1],idxs[2]] = plen
		self.location_boxno[idxs[0],idxs[1],idxs[2]] = 0

		""" Label the neighboring sites as neighbors
		and update the surfaces that are not completely taken up by monomers
		"""
		#sites = IndexUtils.nearest_neighbors(idxs)
		sites = IndexUtils.nearest_effected(idxs)
		#self.update_polymerization(idxs)
		initial_surfaces = self.get_regular_surface_count(sites)
                self.update_super_lattice(idxs,'add')
		self.update_neighbor_sites(sites,'add')
		final_surfaces = self.get_regular_surface_count(sites)
		self.surface_site_count += final_surfaces - initial_surfaces
		
                self.energy_level_population[0] += 1
		
		""" Update the simulation time, the total no. of nodes used,
		and the last box whose label was changed
		"""
                self.last_initiation = self.current_time
                self.total_nodes += 1
                self.last_update_box = idxs
                
                check = 1

    def create_polymer_propensity_list(self):
        """ This function computes the propensity of reaction for all the active chain ends present in the system
	this is the sum of the propensity for propagation, cyclization, and termination
	"""
        for chain_no in self.active_chains:
            last_idx_set = self.all_chains[chain_no][-1]
            if self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] == 0:
                print 'Radical location Error: ', last_idx_set, chain_no, self.junction_chain[last_idx_set[0],last_idx_set[1],last_idx_set[2]], sum(sum(sum(self.Nr)))
		sys.stdout.flush()
                sys.exit()

	    # If the radical is in the effected zone then only re-evaluate
            if IndexUtils.get_distance(self.last_update_box,last_idx_set) <= 10.0:
                self.totalProbSet[chain_no], self.neighbor_groups[chain_no], self.non_zero_sites[chain_no], self.neighbor_monomers[chain_no] = self.compute_propensity(last_idx_set, chain_no)
                self.problist[chain_no] = sum(self.totalProbSet[chain_no])

        self.A_polymer = sum(self.problist.values())

    def compute_propensity(self,indx_set, chain_no):
	""" This function computes the propensity of propagation, cyclization, and termination at an active node.
	Args:
	    indx_set (list): Co-ordinates of the node
	    
	Returns:
	    Alist (list), func_values (list): List with non-zero values for next bond sites from a total set of 108
	"""
	# Get the list of next step voxels
	#low_set = IndexUtils.low_sitesfornextmove(indx_set)
	#IndSet = IndexUtils.sitesfornextmove(indx_set)
	center_monomer_type = int(self.location_monomers[indx_set[0],indx_set[1],indx_set[2]])

        # Store the reaction weights for each index site
        #Alist = list(np.zeros(shape=(IndexUtils.maxSites)))
        Alist = []
	site_no = []
	
	# Store the function group values at sites free for propagation
	funcValues = []
	# Store the monomer types at sites free for propagation
	monomer_types = []
        #funcValues = np.zeros(shape=(IndexUtils.maxSites))
	# Maximum propensity
	max_A_polymer = 0.0

        # There are 108 total sites available for the next move
        # We randomly keep on selecting from the all the available sites
        # till we reach the total number of functional groups
        # available for the radical reaches Nm0,
        # which is the maximum available number of functional groups

        # Counter for functional groups
        funcGrpCount = self.Nm0
        
        # Shuffle the list of bonds
	#site_range = rand.sample(xrange(IndexUtils.maxSites),10)
	site_range = []
	this_bond_samples = 0
	
	#factor = self.non_fluctuation_propensity/(self.all_monomers.mean_d_factors*self.fluctuator.mean_lattice_propensity)
	
	#if self.fluctuation_transition_status==False:
	xi_1 = rand.uniform(0.0,1.0)

	if xi_1>self.fluctuator.length_energy_factor:
	    #site_range = IndexUtils.sample_bond_sites()
	    #this_bond_samples = IndexUtils.total_bond_samples
	    site_range = IndexUtils.sample_low_en_sites()
	    this_bond_samples = IndexUtils.low_en_samples
	    
	#    if len(self.all_chains[chain_no])>1:
	#	vec1 = np.array(self.all_chains[chain_no][-2]) - np.array(self.all_chains[chain_no][-1])
	#	
	#	for samp in site_range:
	#	    bond_vec = np.array(IndexUtils.totalList[samp])
	#	    
	#	    cosine = np.dot(vec1,bond_vec)/(np.linalg.norm(vec1,2)*np.linalg.norm(bond_vec,2))
	#	    
	#	    if cosine==-1.0 and xi_2>self.fluctuator.angle_energy_factor:
	#		site_range.remove(samp)
	#		site_range.insert(0,samp)

	if this_bond_samples==0:
	    #print 'Selecting from high energy sites ', xi, self.fluctuator.energy_factor, self.diffusion.phi[0], self.critical_phi_0
	    site_range = IndexUtils.sample_bond_sites()
	    this_bond_samples = IndexUtils.total_bond_samples
	else:
	    site_range = site_range + IndexUtils.sample_high_en_sites()
	    this_bond_samples = IndexUtils.total_bond_samples
	
	IndSet = IndexUtils.nextbondsites(indx_set,site_range)
	
	#print site_range
        #rand.shuffle(IndexUtils.site_range)
	#site_set = rand.sample(IndexUtils.site_range,40)

	# Avoid node
	try:
	    avoid_node = self.all_chains[chain_no][-2]
	except IndexError:
	    avoid_node = None
	
        for k in xrange(0,this_bond_samples):
	    setno = site_range[k]

	#for setno in site_set:
	    if funcGrpCount > 0:
		idxs = IndSet[k]
		if idxs==indx_set:
		    print 'Rad Idx Error: ', idxs, indx_set
		    sys.stdout.flush()
		    sys.exit()
		if idxs!=avoid_node:
		    if self.location_functionalities[idxs[0],idxs[1],idxs[2]]>0:
			this_monomer = int(self.location_monomers[idxs[0],idxs[1],idxs[2]])
			this_propensity = self.all_monomers.k_p[center_monomer_type,this_monomer]*self.location_functionalities[idxs[0],idxs[1],idxs[2]]

			Alist.append(this_propensity)
			funcValues.append(self.location_functionalities[idxs[0],idxs[1],idxs[2]])
			monomer_types.append(this_monomer)
			site_no.append(setno)
			#Alist[setno] += self.r*self.k_p*self.funcGrp2[idxs[0],idxs[1],idxs[2]]
			funcGrpCount += -self.location_functionalities[idxs[0],idxs[1],idxs[2]]
			
			if Alist[-1]>max_A_polymer:
			    max_A_polymer = Alist[-1]
			
			if self.Nr[idxs[0],idxs[1],idxs[2]]>= 1:
			    this_monomer = int(self.location_monomers[idxs[0],idxs[1],idxs[2]])
			    this_propensity = self.all_monomers.k_t[center_monomer_type,this_monomer]*self.Nr[idxs[0],idxs[1],idxs[2]]

			    Alist.append(0.5*this_propensity)
			    funcValues.append(0)
			    monomer_types.append(this_monomer)
			    site_no.append(setno)
			    
			    #Alist[setno] += 0.5*self.k_t*self.Nr[idxs[0],idxs[1],idxs[2]]
			    funcGrpCount += -1
			    
			    if Alist[-1]>max_A_polymer:
				max_A_polymer = Alist[-1]

		    elif self.check_super_lattice(idxs)==8:
			xi_2 = rand.uniform(0.0,1.0)
			if xi_2<self.all_monomers.free_monomer_fraction:
			    this_monomer = int(self.all_monomers.select_a_monomer())
			    this_propensity = self.all_monomers.k_p[center_monomer_type,this_monomer]*self.all_monomers.functionality[this_monomer]
			    Alist.append(this_propensity)
			    site_no.append(setno)
			    monomer_types.append(this_monomer)
			    funcValues.append(self.all_monomers.functionality[this_monomer])
			    funcGrpCount += -self.all_monomers.functionality[this_monomer]
			    
			    if Alist[-1]>max_A_polymer:
				max_A_polymer = Alist[-1]
            else:
                break

	#if funcGrpCount>0:
	#    #xi = xi_1*xi_2/(xi_1+xi_2)
	#    
	#    for k in xrange(0,this_bond_samples):
	#	setno = site_range[k]
	#    #for setno in IndexUtils.site_range:
	#    #for setno in site_set:
	#	if funcGrpCount > 0:
	#	    idxs = IndSet[k]
	#	    
	#	    if idxs==indx_set:
	#		print 'Rad Idx Error: ', idxs, indx_set
	#		sys.stdout.flush()
	#		sys.exit()
	#	    if self.check_super_lattice(idxs)==8:
	#		xi_2 = rand.uniform(0.0,1.0)
	#		if xi_2<self.all_monomers.free_monomer_fraction:
	#		    this_monomer = int(self.all_monomers.select_a_monomer())
	#		    this_propensity = self.all_monomers.k_p[center_monomer_type,this_monomer]*self.all_monomers.functionality[this_monomer]
	#		    Alist.append(this_propensity)
	#		    site_no.append(setno)
	#		    monomer_types.append(this_monomer)
	#		    funcValues.append(self.all_monomers.functionality[this_monomer])
	#		    funcGrpCount += -self.all_monomers.functionality[this_monomer]
	#		    
	#		    if Alist[-1]>max_A_polymer:
	#			max_A_polymer = Alist[-1]
	#	    else:
	#		    break
                 
        return Alist, funcValues, site_no, monomer_types
        
    def get_chain_neighbor_reaction(self):
	""" Uses the polymer reaction propensity data to decide which chain is going to react next
	and the kind of reaction that is expected.
	"""
        # Random number to select which chain is going to react
        xi2 = rand.uniform(1e-10,1.0)
        # Random number to select the type of reaction
        xi3 = rand.uniform(1e-10,1.0)
        # Random number to select the type of reaction
        xi4 = rand.uniform(1e-10,1.0)
        
        temp_sum = 0.0
        
        for chain_no in self.active_chains:
            thisProb = self.problist[chain_no]/self.A_polymer

            # If found select chain no and exit
            if (temp_sum < xi2) and xi2 <= (temp_sum+thisProb):
                selected_chain = chain_no
                break
            else:
                temp_sum += thisProb

        # Get indices of chain ends
        Alist = self.totalProbSet[selected_chain]
        funcValues = self.neighbor_groups[selected_chain]
	site_set = self.non_zero_sites[selected_chain]
	
        last_idx_set = self.all_chains[selected_chain][-1]
	last_monomer = int(self.location_monomers[last_idx_set[0],last_idx_set[1],last_idx_set[2]])
	# Get the second last index set
	try:
	    prev_idx = self.all_chains[selected_chain][-2]
	#    if self.blocked_matrix[prev_idx[0],prev_idx[1],prev_idx[2]]==0:
	#	print 'Error: ', prev_idx, self.all_chains[selected_chain]
	#	sys.stdout.flush()
	#	sys.exit()
	except IndexError:
	    prev_idx = None
	    
        # Get the list of next step voxels
        #IndSet = IndexUtils.sitesfornextmove(last_idx_set)
	IndSet = IndexUtils.nextbondsites(last_idx_set,site_set)
	
	if last_idx_set in IndSet:
	    print 'Index Error: ', selected_chain
	    sys.stdout.flush()
	    sys.exit()

        # Now find which neighboring box is going to react
        sumAlist = sum(Alist)
          
        temp_sum = 0.0

	#diff = -1.0
	#
	#while diff < 0.0:
	#    bno = rand.randint(0,IndexUtils.maxSites-1)
	#    
	#    diff = Alist[bno]/self.max_A[selected_chain] - xi3

	for bno in range(0,len(site_set)):
            thisProb = Alist[bno]/sumAlist

            # If found select chain no and exit
            if (temp_sum < xi3) and xi3 <=(temp_sum+thisProb):
                selected_box = self.non_zero_sites[selected_chain][bno]
		funcValue = self.neighbor_groups[selected_chain][bno]
		monomer = int(self.neighbor_monomers[selected_chain][bno])
		idx = IndexUtils.nextbondsites(last_idx_set,[site_set[bno]])[0]
		break
            else:
                temp_sum += thisProb

        # Now decide whether reaction is propagation or termination
        #idx = IndSet[selected_box]
         
        propWt, cycleWt, termWt = 0.0, 0.0, 0.0

	if funcValue==self.all_monomers.functionality[int(monomer)]:
	    propWt = self.all_monomers.k_p[last_monomer,monomer]*funcValue
	elif funcValue>0:
	    cycleWt = self.all_monomers.k_p[last_monomer,monomer]*funcValue

        termWt = self.all_monomers.k_t[last_monomer,monomer]*self.Nr[idx[0],idx[1],idx[2]]

        # Decide if there is any reaction, and whether the reaction is 
        # propagation and termination
        totalWt = propWt + cycleWt + termWt

        if totalWt == 0.0:
            print 'Something wrong', Alist, bno, Alist[bno], self.problist[selected_chain], funcValue
	    sys.stdout.flush()
            sys.exit()
            react = None    
        else:
            propFac = propWt/totalWt
            cycleFac = (propWt+cycleWt)/totalWt
            
            if xi4 <= propFac:
                react = 'prop'
                # One functional group taken at idx
		self.location_functionalities[idx[0],idx[1],idx[2]] = funcValue - 1
                self.transientFuncs += -1
                self.transientMonomers += -1
                self.total_nodes += 1
  
                self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] += -1
                self.Nr[idx[0],idx[1],idx[2]] += 1
		
                self.location_chain[idx[0],idx[1],idx[2]] = selected_chain
		self.location_monomers[idx[0],idx[1],idx[2]] = monomer
		self.location_boxno[idx[0],idx[1],idx[2]] = len(self.all_chains[selected_chain])

                if IndexUtils.get_distance2(idx,last_idx_set) not in IndexUtils.length_set:
                    print 'prop error: ', idx, last_idx_set
		    sys.stdout.flush()
                    sys.exit()

		#sites = IndexUtils.nearest_neighbors(idx)
		#self.update_polymerization(idx)
		sites = IndexUtils.nearest_effected(idx)
		
		initial_surfaces = self.get_regular_surface_count(sites)
		self.update_super_lattice(idx,'add')
		self.update_neighbor_sites(sites,'add')
		final_surfaces = self.get_regular_surface_count(sites)
		self.surface_site_count += final_surfaces - initial_surfaces
		
                #self.update_neighbors(idx,self.neighbor_set,site_values,'add')
		self.energy_level_population[0] += 1

		self.last_update_box = idx
	    elif xi4 <= cycleFac:
                react = 'cross'
                # Use up the cyclization functional group
		self.location_functionalities[idx[0],idx[1],idx[2]] += -1
                self.transientFuncs += -1
                self.totalJunctionPoints += 1
                
                # Update junction points dict
                self.junction_chain[idx[0],idx[1],idx[2]] = selected_chain
                self.junction_list.append(idx)
		
		if self.location_chain[idx[0],idx[1],idx[2]]!=selected_chain:
		    true_junction = 1
		    self.realJunctionPoints += 1
		else:
		    true_junction = self.check_for_true_junction(selected_chain,idx)
		    self.realJunctionPoints += true_junction
		
		self.junction_boxno[idx[0],idx[1],idx[2]] = len(self.all_chains[selected_chain])
		
		#self.junctions_in_chains[selected_chain] += 1
		#self.junctions_in_chains[self.location_chain[idx[0],idx[1],idx[2]]] += 1

		self.update_chain_connectivity(self.location_chain[idx[0],idx[1],idx[2]],selected_chain)
		
		if IndexUtils.get_distance2(idx,last_idx_set) not in IndexUtils.length_set:
                    print 'prop error: ', idx, last_idx_set
		    sys.stdout.flush()
                    sys.exit()

		#sites = IndexUtils.sites_for_thermal_moves(idx)
		#sites = IndexUtils.nearest_neighbors(idx)
		sites = IndexUtils.nearest_effected(idx)
		
		initial_surfaces = self.get_junction_surface_count(sites)
		#self.update_super_lattice(idx,'add')
		self.update_junction_super_lattice(idx,'add')
		
		if true_junction==1:
		    self.update_junction_neighbors(sites,'add')
		    final_surfaces = self.get_junction_surface_count(sites)
		    self.junction_surface_count += (final_surfaces - initial_surfaces)
		
		#self.junction_interaction_count += self.update_junction_neighbors(idx,'add')
		
                # Now move the radical    
                if self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] == 0:
                    print 'Error: radical move'
		    sys.stdout.flush()
                    sys.exit()  

                self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] += -1
                self.Nr[idx[0],idx[1],idx[2]] += 1
                self.last_update_box = idx
            elif termWt > 0.0:
                react = 'term' 
                # Now move the radical      
                self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] += -1
                self.Nr[idx[0],idx[1],idx[2]] += -1
                
                self.last_update_box = idx
        
        return selected_chain, idx, monomer, react
    
    def update_chain_connectivity(self,chain_1,chain_2):
	if chain_2 not in self.chain_connectivity[chain_1]:
	    self.chain_connectivity[chain_1].append(chain_2)
	
	if chain_1 not in self.chain_connectivity[chain_2]:
	    self.chain_connectivity[chain_2].append(chain_1)

    def check_other_connection(self,chain_1):
	status = 'not connected'
	
	for chain_2 in self.chain_connectivity[chain_1]:
	    if chain_2 != chain_1:
		status = 'connected'
		break
	
	return status

    def check_for_true_junction(self,chain,idx):
	box_no = self.all_chains[chain].index(idx)
	truth_value = 0
	
	for i in xrange(box_no+1,len(self.all_chains[chain])-1):
	    if self.junction_chain[idx[0],idx[1],idx[2]]>0:# and self.junction_chain[idx[0],idx[1],idx[2]]!=chain:
		truth_value = 1
		break
	
	return truth_value

    def terminate_a_pair(self):
	# Pick to chains
	# Delete everything about the second chain
	# Extend the first chain number of times the length of the second
	pair = []
	
	#print 'picking pairs ', self.chain_connectivity, self.active_chains
	#sys.stdout.flush()
	
	while len(pair)<2:
	    n1 = rand.sample(self.active_chains,1)[0]
	    if n1 not in pair:
		if len(self.chain_connectivity[n1])==0:
		    pair.append(n1)
		elif self.check_other_connection(n1)=='not connected':
		    pair.append(n1)
	
	#len_of_1 = len(self.all_chains[pair[0]])
	#len_of_2 = len(self.all_chains[pair[1]])
	idx = self.all_chains[pair[1]][-1]
	self.Nr[idx[0],idx[1],idx[2]] += -1
	
	self.active_chains.remove(pair[1])
	del self.totalProbSet[pair[1]]
        del self.neighbor_monomers[pair[1]]
	del self.neighbor_groups[pair[1]]
        del self.problist[pair[1]]

	#self.delete_a_chain(pair[1])
	#
	#for step_no in xrange(0,len_of_2):
	#    idx, monomer = self.forced_extension(pair[0])
	#    self.all_chains[pair[0]].append(idx)
	    
	idx = self.all_chains[pair[0]][-1]
	self.Nr[idx[0],idx[1],idx[2]] += -1
	
	self.active_chains.remove(pair[0])
	del self.totalProbSet[pair[0]]
        del self.neighbor_monomers[pair[0]]
        del self.neighbor_groups[pair[0]]
        del self.problist[pair[0]]

	self.early_terminations.append(pair[0])
	self.early_terminations.append(pair[1])
	
	#self.all_chains[pair[1]] = self.all_chains[pair[0]][len_of_1:]
	#temp_set = self.all_chains[pair[0]][:len_of_1]
	#del self.all_chains[pair[0]]
	#self.all_chains[pair[0]] = temp_set
	#
	#for idx in self.all_chains[pair[1]]:
	#    self.location_chain[idx[0],idx[1],idx[2]] = pair[1]
	
	#self.update_chain_numbers(pair[1])

    def get_chain_to_terminate(self,thischain):
        allChains = []
        idxSet = self.all_chains[thischain][-1]
        # Loop over all chains to know the ones ending in this box
        for chain_no in self.active_chains:
            if chain_no != thischain:
                # Get indices of chain ends
                last_idx_set = self.all_chains[chain_no][-1]
                
                if idxSet == last_idx_set:
                    allChains.append(chain_no)
        
        if len(allChains)>0:
            termChain = rand.choice(allChains)
        else:
            termChain = None
        
        return termChain
        
    def close_chains(self,chain1,chain2):
        self.active_chains.remove(chain1)
        self.active_chains.remove(chain2)
    
        self.terminated_pairs[chain1] = chain2
        self.terminated_pairs[chain2] = chain1

        del self.totalProbSet[chain1]
        del self.totalProbSet[chain2]
        del self.neighbor_monomers[chain1]
        del self.neighbor_monomers[chain2]
	del self.neighbor_groups[chain1]
        del self.neighbor_groups[chain2]
        del self.problist[chain1]
        del self.problist[chain2]
    
    def compute_chain_stat(self):
        chainLenList = []

        self.mean = 0.0
        self.var = 0.0

        for aChain in range(1,len(self.all_chains)+1):
            chainLenList.append(len(self.all_chains[aChain]))

        if len(chainLenList) != 0:    
            self.mean = sum(chainLenList)/len(chainLenList)
        
            for chainL in chainLenList:
                self.var += (chainL-self.mean)**2
            
            self.var = self.var/len(chainLenList)          

    def compute_doc(self,step_type):
        self.DoC = (self.initialFuncs-self.transientFuncs)/float(self.initialFuncs)
        
#        if self.last_polymer_propensity==0.0 and self.global_propensities['polymerize']>0.0:
#	    self.polymerfile = open('polydata'+str(self.proc_no)+'.csv','a')
#	    print >> self.polymerfile, self.current_time, ',', self.DoC, ',', self.global_propensities['polymerize']
#	    self.polymerfile.close()
#
#	self.last_polymer_propensity=self.global_propensities['polymerize']
#	
#	self.last_node_size = self.total_nodes

    def write_intermediate_output(self,step_type):
	self.DCoutfile = open('Time_DC'+str(self.proc_no)+'.csv','a')
	self.Chainoutfile = open('chainstat'+str(self.proc_no)+'.csv','a')
	self.event_file = open('eventstat'+str(self.proc_no)+'.csv','a')
	self.fluctuation_file = open('fluctuation'+str(self.proc_no)+'.csv','a')
	self.site_file = open('site_fractions'+str(self.proc_no)+'.csv','a')
	self.g_bar_file = open('g_bars'+str(self.proc_no)+'.csv','a')
	self.diffusivity_file = open('diffusivity'+str(self.proc_no)+'.csv','a')
	self.propensity_file = open('propensities'+str(self.proc_no)+'.csv','a')
	self.energy_file = open('energy'+str(self.proc_no)+'.csv','a')
	self.bond_length_file = open('bond_lengths'+str(self.proc_no)+'.csv','a')
	self.heat_file = open('heat_info'+str(self.proc_no)+'.csv','a')
	#self.modulus_file = open('modulus'+str(self.proc_no)+'.csv','a')
	self.trial_file = open('fluctuation_trials'+str(self.proc_no)+'.csv','a')
	self.monomer_file = open('monomer_count'+str(self.proc_no)+'.csv','a')
	self.leap_file = open('fluctuation_epr'+str(self.proc_no)+'.csv','a')
	self.gel_file = open('gel_transient'+str(self.proc_no)+'.csv','a')
	self.movable_file = open('movable_nodes'+str(self.proc_no)+'.csv','a')
	self.component_file = open('component_size'+str(self.proc_no)+'.csv','a')
	
	print >> self.DCoutfile, self.step_no, ',', self.current_time, ',', self.DoC, ',', self.realJunctionPoints/float(self.max_monomers),',', (self.totalJunctionPoints-self.realJunctionPoints)/float(self.max_monomers)
	print >> self.Chainoutfile, len(self.all_chains),',',len(self.active_chains),',',self.total_nodes,',',self.totalJunctionPoints
	print >> self.fluctuation_file, self.current_time,',',self.global_diffusion_propensity,',',self.global_propensities['fluctuate']
	print >> self.site_file, self.current_time,',',self.DoC,',',self.diffusion.pop_out_string
	print >> self.g_bar_file, self.current_time,',',self.DoC,',',self.diffusion.g_bar_string
	print >> self.diffusivity_file, self.current_time,',',self.DoC,',',self.total_nodes,',',self.all_monomers.mean_d_factors,',',self.lattice_propensity
	print >> self.energy_file, self.current_time,',',self.DoC,',',self.total_nodes,',',(self.total_network_energy/float(self.total_nodes))
	print >> self.bond_length_file, self.current_time,',',self.DoC,',',self.bond_pops[4.0],',',self.bond_pops[5.0],',',self.bond_pops[6.0],',',self.bond_pops[9.0],',',self.bond_pops[10.0]
	print >> self.heat_file, self.current_time,',',self.DoC,',',self.heat_production,',',self.current_temperature
	#print >> self.modulus_file, self.current_time,',',self.DoC,',',self.MaterialModel.E_infi,',',self.MaterialModel.relax_tc
	print >> self.trial_file, self.current_time,',',self.DoC,',',self.fluctuation_trials['success'],',',self.fluctuation_trials['failure']
	
	if self.trajectory_samples>0:
	    print >> self.leap_file, self.current_time, ',', self.DoC, ',', self.fluctuation_epr, ',', self.fluctuation_sum, ',', self.global_propensities['fluctuate'],',',(self.trajectory_entropy/float(self.trajectory_samples))
	print >> self.gel_file, self.current_time, ',', self.DoC, ',', self.diffusion.phi[0], ',', self.critical_phi_0
	print >> self.component_file, self.current_time, ',', self.DoC, ',', self.max_connected_unit, ',', self.last_unit
	
	f1 = self.total_highly_movable_nodes/float(self.total_nodes)
	f2 = (self.total_movable_nodes-self.total_highly_movable_nodes)/float(self.total_nodes)
	
	print >> self.movable_file, self.current_time, ',', self.DoC, ',', f1, ',', f2, ',', (1-f1-f2)
	
	out_string_prop = str(self.current_time)+','+str(self.DoC)
	out_string_events = str(self.current_time)+','+str(self.DoC)
	out_string_monomers = str(self.current_time)+','+str(self.DoC)

	for event_type in self.global_propensities.keys():
	    if event_type=='fluctuate':
		out_string_prop += ','+str(self.all_monomers.mean_d_factors*self.global_propensities[event_type])

	    out_string_events += ','+str(self.event_counter[event_type])
	    
	for mo_type in xrange(0,self.all_monomers.no_monomers):
	    out_string_monomers += ','+str(self.all_monomers.added_set[mo_type]/sum(self.all_monomers.added_set))
	    
	print >> self.propensity_file, out_string_prop
	print >> self.event_file, out_string_events
	print >> self.monomer_file, out_string_monomers
	
	self.DCoutfile.close()
	self.Chainoutfile.close()
	self.event_file.close()
	self.fluctuation_file.close()
	self.site_file.close()
	self.g_bar_file.close()
	self.diffusivity_file.close()
	self.propensity_file.close()
	self.energy_file.close()
	self.bond_length_file.close()
	self.heat_file.close()
	#self.modulus_file.close()
	self.trial_file.close()
	self.monomer_file.close()
	self.leap_file.close()
	self.gel_file.close()
	self.movable_file.close()
	self.component_file.close()
	
	if self.write_video=='Yes':
	    self.write_vtk_files()
	    self.write_final_vtk_files()

    def write_chains(self):
        for chain in self.closedChains.keys():
            outfile = open(str(chain)+'chaincoords.txt','w')
            for listno in range(0,len(self.closedChains[chain])):
                print >> outfile, self.closedChains[chain][listno][0], self.closedChains[chain][listno][1], self.closedChains[chain][listno][2]

            outfile.close()
    
    def compute_max_connected_unit(self):
	self.max_connected_unit = 0
	self.last_unit = 0
	total_length = 0
	
	for chain in self.all_chains.keys():
	    local_length = len(self.all_chains[chain])
	    total_length += local_length
	    
	    if len(self.chain_connectivity[chain])>0:
		for i in self.chain_connectivity[chain]:
		    if i!=chain:
			local_length += len(self.all_chains[i])
		    
	    if local_length>self.max_connected_unit:
		self.last_unit = self.max_connected_unit

	    self.max_connected_unit = max(local_length,self.max_connected_unit)
	
	self.last_unit = total_length - self.max_connected_unit

    def write_densities(self):
        outfile = open('DensityFile.txt','w')
        
        totalChains = len(self.closedChains)
        totalJuncs = len(self.closedJuncPoints)
        avgJuncPoints = totalJuncs/(3*totalChains)
        
        print >> outfile , totalChains, totalJuncs, avgJuncPoints
        
        outfile.close()

    def write_pixels(self):
        self.pixelCount += 1
                
        outfile_temp = open('crosspoly'+str(self.pixelCount)+'.txt','w')
        
        dSize = int(self.lattice_size)
        
        #for idx in range(0,dSize):
        #    for idy in range(0,dSize):
        #        for idz in range(0,dSize):
        #            iset = [idx,idy,idz]
        #            if iset in self.activeChains[0]:
        #                if self.neighbor_sites[idx,idy,idz] == 1:
        #                    print >> outfile_temp, '1'
        #                elif self.junctionPoints[idx,idy,idz] == 1:
        #                    print >> outfile_temp, '3'
        #                elif self.blocked_matrix[idx,idy,idz] == 1:
        #                    print >> outfile_temp, '2'
        #            else:
        #                print >> outfile_temp, '0'
        
        outfile_temp.close()

    #def compute_directional_shrinkage(self):
    #    self.maxBonds = (self.lattice_size**2)*(self.lattice_size+1)
    #    
    #    self.shrinkage = np.zeros(shape=(3))
    #    
    #    # Each bond accounts for a reduction from 1.519 A to 1.33 A
    #    for did in range(0,3):
    #        self.shrinkage[did] = (1.519-1.33)*self.directionBonds[did]
    #
    #    # Domain size
    #    domain_size = (0.5*self.initialFuncs*self.MolVol)**(1.0/3.0)
    #
    #    for did in range(0,3):
    #        self.shrinkage[did] *= (1.0/domain_size)
    #        # Average shrinkage in each direction
    #        self.shrinkage[did] *= (1.0/(self.lattice_size**2))
    #
    #    print self.shrinkage, self.directionBonds, domain_size
        
    # This function checks if this box is free for a monomer
    # Or if each of its 8 vertices are not-occupied by the polymer chain
    def check_if_empty(self,idxs):
        check = 'Yes'

	if self.check_super_lattice(idxs)=='not empty':
	    check = 'No'
        
        #checkSum = self.blocked_matrix[idxs[0],idxs[1],idxs[2]] + self.neighbor_sites[idxs[0],idxs[1],idxs[2]]
        #
        #if checkSum > 0:
        #    check = 'No'
        
        return check        

    def make_thermal_move(self):
	status = 'failure'
	
	idxs = self.select_thermal_unit_fast()
	status = self.move_this_node(idxs)

	return status

    def pick_a_radical_node(self):
	chain_no = rand.choice(self.active_chains)
	
	return self.all_chains[chain_no][-1]

    def move_this_node(self,idxs):
	chain_no = self.location_chain[idxs[0],idxs[1],idxs[2]]
	box_no = int(self.location_boxno[idxs[0],idxs[1],idxs[2]])
	#box_no2 = self.all_chains[chain_no].index(idxs)
	#if box_no!=box_no2:
	#    print 'Move this node: ', chain_no, box_no, box_no2
	#    sys.exit()
	
	lidx = self.get_lidx_from_cubidx(idxs)
	trials = 0
	
	# Initialize bond update dict
	self.bond_pop_update[4.0], self.bond_pop_update[5.0], self.bond_pop_update[6.0] = 0, 0, 0
	self.bond_pop_update[9.0], self.bond_pop_update[10.0] = 0, 0
	
	# Default trial status
	status = 'failure'
	check = 'bad move'
	
	#while self.movability[idxs[0],idxs[1],idxs[2]]>0 and check=='bad move' and status=='failure':
	# Default energy allowed move
	#jump_probability, move_dict, move_count, move_prob_dict, max_prob = self.compute_jump_prob(idxs)
	move_set = range(0,6)
	rand.shuffle(move_set)
	move_idx = 0
	
	xi = rand.uniform(0.0,(1.0-1e-3))
	
	while status=='failure' and move_idx<6: #and self.movability[idxs[0],idxs[1],idxs[2]]>0:
	    # Select random number for trial
	    
	    energy_decision = False
	
	    # Select a move
	    move_no = move_set[move_idx]
	    value = self.movable_dict[lidx,move_set[move_idx]]
	    move_idx += 1
	    
	    if value>0:
		trial_move = IndexUtils.thermalMoves[move_no]

		check, second_check = 'not empty', 'not empty'

		trial_move = IndexUtils.thermalMoves[move_no]
		
		if (self.moving_probs[lidx][move_no]/self.max_moving_probs[lidx])>xi:
		    energy_decision = True
		    
		    if self.fluctuation_transition_status==False:
			#self.fluctuation_samples.append(self.moving_probs[lidx][move_no]/sum(self.moving_probs[lidx]))
			self.fluctuation_samples.append(self.moving_probs[lidx][move_no]/self.global_propensities['fluctuate'])
			
			self.reverse_move = 3 + (2-move_no)
			
			self.bond_populations['high'] += self.bond_motion_update[lidx,move_no]
			self.bond_populations['low'] += -self.bond_motion_update[lidx,move_no]


		new_pos = self.get_new_position(idxs,trial_move)
		
		if self.junction_chain[idxs[0],idxs[1],idxs[2]] > 0:
		    if self.location_chain[idxs[0],idxs[1],idxs[2]] == chain_no:
			second_chain = self.junction_chain[idxs[0],idxs[1],idxs[2]]
			second_box = int(self.junction_boxno[idxs[0],idxs[1],idxs[2]])
		    else:
			second_chain = self.location_chain[idxs[0],idxs[1],idxs[2]]
			second_box = int(self.location_boxno[idxs[0],idxs[1],idxs[2]])

		    if energy_decision==True:
			self.update_chain_position(chain_no,box_no,new_pos)
			self.update_second_chain_position(second_chain,second_box,new_pos)

			self.location_chain[new_pos[0],new_pos[1],new_pos[2]] = chain_no
			self.location_chain[idxs[0],idxs[1],idxs[2]] = 0
			
			self.location_boxno[new_pos[0],new_pos[1],new_pos[2]] = self.location_boxno[idxs[0],idxs[1],idxs[2]]
			self.location_boxno[idxs[0],idxs[1],idxs[2]] = -1

			self.junction_chain[new_pos[0],new_pos[1],new_pos[2]] = second_chain
			self.junction_chain[idxs[0],idxs[1],idxs[2]] = 0
			
			self.junction_boxno[new_pos[0],new_pos[1],new_pos[2]] = self.junction_boxno[idxs[0],idxs[1],idxs[2]]
			self.junction_boxno[idxs[0],idxs[1],idxs[2]] = -1
			
			self.junction_list.remove(idxs)
			self.junction_list.append(new_pos)
			
			if self.thermal_propensities[lidx]>=1.0:
			    self.highly_movable_nodes.remove(idxs)
			    self.total_highly_movable_nodes += -1
			    self.highly_movable_propensities += -self.thermal_propensities[lidx]
			    self.high_propensities.remove(self.thermal_propensities[lidx])

			self.movable_nodes.remove(idxs)
			
			self.global_propensities['fluctuate'] += -self.thermal_propensities[lidx]
			self.fluctuation_quanta[int(self.thermal_propensities[lidx])] += -1
			self.thermal_propensities[lidx] = 0
			
			for k in xrange(0,6):
			    self.movable_dict[lidx,k] = 0
			    self.moving_probs[lidx,k] = 0
			    
			self.max_moving_probs[lidx] = 0
			self.total_movable_nodes += -1
			
			status = 'success'
			
			self.last_update_box = new_pos
			self.last_move = trial_move
		else:
		    if energy_decision==True:
			self.update_chain_position(chain_no,box_no,new_pos)
			
			self.location_chain[new_pos[0],new_pos[1],new_pos[2]] = chain_no
			self.location_chain[idxs[0],idxs[1],idxs[2]] = 0
			
			self.location_boxno[new_pos[0],new_pos[1],new_pos[2]] = self.location_boxno[idxs[0],idxs[1],idxs[2]]
			self.location_boxno[idxs[0],idxs[1],idxs[2]] = -1
			
			if self.thermal_propensities[lidx]>=1.0:
			    self.highly_movable_nodes.remove(idxs)
			    self.total_highly_movable_nodes += -1
			    self.highly_movable_propensities += -self.thermal_propensities[lidx]
			    self.high_propensities.remove(self.thermal_propensities[lidx])

			self.movable_nodes.remove(idxs)
			
			self.global_propensities['fluctuate'] += -self.thermal_propensities[lidx]
			self.fluctuation_quanta[int(self.thermal_propensities[lidx])] += -1
			self.thermal_propensities[lidx] = 0
			
			for k in xrange(0,6):
			    self.movable_dict[lidx,k] = 0
			    self.moving_probs[lidx,k] = 0
			self.max_moving_probs[lidx] = 0
			
			self.total_movable_nodes += -1
			
			status = 'success'
			
			self.last_update_box = new_pos
			self.last_move = trial_move

	    self.fluctuation_trials[status] += 1

	return status

#    def decide_fluctuation(self,xi,idx,energy_change):
#	node_type = self.get_node_type(idx)
#	associated_energy, node_type = self.get_associated_energy(idx)
#	
#	transition_odds = math.exp(-energy_change/self.fluctuator.beta)
#	
#	if xi<=transition_odds:
#	#if energy_change<0.0:
#	    decision = True
#	    self.total_network_energy += energy_change
#	    self.update_bond_population()
#	    
#	    if self.fluctuation_transition_status==False:
#		epr_value = min(1.0,transition_odds)/min(1.0,1.0/transition_odds)
#		
#		self.fluctuation_samples.append(epr_value)
#		self.fluctuation_times.append(self.current_time)
#	else:
#	    decision = False
#	
#	if node_type=='isolated':
#	    decision = True
#
#	return decision

    def get_new_position(self,thisbox,trial_move):
    	new_pos = [0,0,0]

        #thisbox = self.all_chains[chain_no][box_no]
        for i in range(0,3):
            new_pos[i] = thisbox[i] + trial_move[i]

            if new_pos[i] >= self.lattice_size:
                new_pos[i] = new_pos[i] - self.lattice_size
            elif new_pos[i] < 0:
                new_pos[i] = self.lattice_size + new_pos[i]

	#check = self.check_super_lattice_thermal(thisbox,new_pos)
	
	#if check!='empty':
	#    print check, sum(sum(sum(self.super_lattice)))/float(self.total_nodes)
	#    sys.exit()
	    
	return new_pos

    def update_bond_population(self):
	for key in self.bond_pops.keys():
	    self.bond_pops[key] += self.bond_pop_update[key]

    def check_for_motion(self,chain_no,box_no,trial_move):
	"""This function computes the change in bond length and angle energy for each move
	    ----boxlow_2----boxlow_1----center----boxup_1----boxup_2----
	"""
	new_pos = [0,0,0]

        U_length_1, U_length_2 = 0.0, 0.0
	U_angle_low, U_angle_up, U_angle_c = 0.0, 0.0, 0.0

        thisbox = self.all_chains[chain_no][box_no]
        for i in range(0,3):
            new_pos[i] = thisbox[i] + trial_move[i]

            if new_pos[i] >= self.lattice_size:
                new_pos[i] = new_pos[i] - self.lattice_size
            elif new_pos[i] < 0:
                new_pos[i] = self.lattice_size + new_pos[i]
                    
        check = self.check_super_lattice_thermal(thisbox,new_pos)

        boxlow_1, boxlow_2, boxup_1, boxup_2 = None, None, None, None
        
        if check=='empty':
            if box_no == 0:
                boxlow_1, boxlow_2 = None, None
            else:
		try:
		    boxlow_1 = self.all_chains[chain_no][box_no-1]
		except IndexError:
		    boxlow_1 = None
		    
		if box_no-1==0:
		    boxlow_2 = None
		else:
		    try:
			boxlow_2 = self.all_chains[chain_no][box_no-2]
		    except IndexError:
			boxlow_2 = None

            if box_no == len(self.all_chains[chain_no])-1:
                if chain_no not in self.active_chains and chain_no not in self.early_terminations:
                    other_chain = self.terminated_pairs[chain_no]
                    
		    try:
			boxup_1 = self.all_chains[other_chain][-2]
		    except IndexError:
			boxup_1 = None

		    try:
			boxup_2 = self.all_chains[other_chain][-3]
		    except IndexError:
			boxup_2 = None
                else:
                    boxup_1, boxup_2 = None, None
            else:
                boxup_1 = self.all_chains[chain_no][box_no+1]
		if box_no+1 == len(self.all_chains[chain_no])-1:
		    if chain_no not in self.active_chains and chain_no not in self.early_terminations:
			other_chain = self.terminated_pairs[chain_no]
			
			try:
			    boxup_2 = self.all_chains[other_chain][-2]
			except IndexError:
			    boxup_2 = None
		else:
		    try:
			boxup_2 = self.all_chains[chain_no][box_no+2]
		    except IndexError:
			boxup_2 = None
		
            # Check if this move is compliant with bond distances
            check, dist1, dist2 = IndexUtils.check_motion(new_pos,boxlow_1,boxup_1)

            if check=='good move' and self.fluctuator.bond_energy>0.0:
		#if boxlow_1!=None and boxup_1!=None:
		#    U_angle_c = self.fluctuator.get_angle_energy_change(boxlow_1,boxup_1,thisbox,new_pos)
		    
                if boxlow_1!=None:
                    U_length_1, r_new, r_old = self.fluctuator.get_bond_energy_change(thisbox,new_pos,boxlow_1)

		#    if boxlow_2!=None:
		#	angle_low_i = self.fluctuator.get_angle(boxlow_2,boxlow_1,thisbox)
		#	angle_low_f = self.fluctuator.get_angle(boxlow_2,boxlow_1,new_pos)
		#	U_angle_low = self.fluctuator.get_angle_energy(angle_low_f)-self.fluctuator.get_angle_energy(angle_low_i)

                if boxup_1!=None:
                    U_length_2, r_new, r_old = self.fluctuator.get_bond_energy_change(thisbox,new_pos,boxup_1)

		#    if boxup_2!=None:
		#	angle_up_i = self.fluctuator.get_angle(boxup_2,boxup_1,thisbox)
		#	angle_up_f = self.fluctuator.get_angle(boxup_2,boxup_1,new_pos)
		#	U_angle_up = self.fluctuator.get_angle_energy(angle_up_f)-self.fluctuator.get_angle_energy(angle_up_i)

        #return check, U_length_1, U_length_2, U_angle_c, U_angle_low, U_angle_up, boxlow_1, boxup_1
	return check, U_length_1, U_length_2, 0.0, 0.0, 0.0, boxlow_1, boxup_1


    def update_chain_position(self,chainNo,boxno,new_pos):	
	oldIdx = self.all_chains[chainNo][boxno]
        self.all_chains[chainNo][boxno] = new_pos
	
	self.location_functionalities[new_pos[0],new_pos[1],new_pos[2]] = self.location_functionalities[oldIdx[0],oldIdx[1],oldIdx[2]]
	self.location_functionalities[oldIdx[0],oldIdx[1],oldIdx[2]] = 0
	
	self.location_monomers[new_pos[0],new_pos[1],new_pos[2]] = self.location_monomers[oldIdx[0],oldIdx[1],oldIdx[2]]
	self.location_monomers[oldIdx[0],oldIdx[1],oldIdx[2]] = 0
	
	try:
	    if self.all_chains[chainNo][boxno]==self.all_chains[chainNo][boxno-1] and len(self.all_chains[chainNo])>1:
		print 'Bad chain move: ', self.all_chains[chainNo][boxno], self.all_chains[chainNo][boxno-1], boxno, (boxno-1), len(self.all_chains[chainNo])
		sys.stdout.flush()
		sys.exit()
	except IndexError:
	    pass
	
	try:
	    if self.all_chains[chainNo][boxno]==self.all_chains[chainNo][boxno+1] and len(self.all_chains[chainNo])>1:
		print 'Bad chain move: ', self.all_chains[chainNo][boxno], self.all_chains[chainNo][boxno+1], boxno, (boxno-1), len(self.all_chains[chainNo])
		sys.stdout.flush()
		sys.exit()
	except IndexError:
	    pass
	
	if boxno==(len(self.all_chains[chainNo])-1) and self.Nr[oldIdx[0],oldIdx[1],oldIdx[2]]==0 and chainNo not in self.early_terminations:
            print self.terminated_pairs
	    other_chain = self.terminated_pairs[chainNo]
            self.all_chains[other_chain][-1] = new_pos

	if self.Nr[oldIdx[0],oldIdx[1],oldIdx[2]]>0:
	    #print ''
	    #print 'Radical moved', oldIdx, self.new_pos, chainNo
	    #print ''
	    self.Nr[new_pos[0],new_pos[1],new_pos[2]] = self.Nr[oldIdx[0],oldIdx[1],oldIdx[2]]
	    self.Nr[oldIdx[0],oldIdx[1],oldIdx[2]] = 0

	#neighbor_set_old = IndexUtils.nearest_neighbors(oldIdx)
	neighbor_set_old = IndexUtils.nearest_effected(oldIdx)
	#initial_surfaces_old = self.get_regular_surface_count(neighbor_set_old)
	
	#neighbor_set_new = IndexUtils.nearest_neighbors(new_pos)
	neighbor_set_new = IndexUtils.nearest_effected(new_pos)
	#initial_surfaces_new = self.get_regular_surface_count(neighbor_set_new)
	
	total_sites = self.get_total_sites(neighbor_set_old,neighbor_set_new)
	initial_surfaces = self.get_regular_surface_count(total_sites)

        self.update_super_lattice(oldIdx,'remove')
	self.update_neighbor_sites(neighbor_set_old,'remove')
	
	self.update_super_lattice(new_pos,'add')
	self.update_neighbor_sites(neighbor_set_new,'add')
        #self.update_neighbors(self.new_pos,self.neighbor_set_new,site_values,'add')

	final_surfaces = self.get_regular_surface_count(total_sites)
	#final_surfaces_old = self.get_regular_surface_count(neighbor_set_old)
	#final_surfaces_new = self.get_regular_surface_count(neighbor_set_new)
	
	self.surface_site_count += final_surfaces - initial_surfaces
	
	#self.surface_site_count += final_surfaces_new + final_surfaces_old - initial_surfaces_new - initial_surfaces_old
	
    def update_second_chain_position(self,chainNo,boxno,new_pos):
	oldIdx = self.all_chains[chainNo][boxno]
        self.all_chains[chainNo][boxno] = new_pos
    
	try:
	    if self.all_chains[chainNo][boxno]==self.all_chains[chainNo][boxno-1] and len(self.all_chains[chainNo])>1:
		print 'Bad chain move: ', self.all_chains[chainNo][boxno], self.all_chains[chainNo][boxno-1], boxno, (boxno-1), len(self.all_chains[chainNo])
		sys.stdout.flush()
		sys.exit()
	except IndexError:
	    pass
	
	try:
	    if self.all_chains[chainNo][boxno]==self.all_chains[chainNo][boxno+1] and len(self.all_chains[chainNo])>1:
		print 'Bad chain move: ', self.all_chains[chainNo][boxno], self.all_chains[chainNo][boxno+1], boxno, (boxno-1), len(self.all_chains[chainNo])
		sys.stdout.flush()
		sys.exit()
	except IndexError:
	    pass
	
	jn_neighbors_old = IndexUtils.nearest_effected(oldIdx)
	#jn_initial_surfaces_old = self.get_junction_surface_count(jn_neighbors_old)
	
	jn_neighbors_new = IndexUtils.nearest_effected(new_pos)
	#jn_initial_surfaces_new = self.get_junction_surface_count(jn_neighbors_new)
	
	total_sites = self.get_total_sites(jn_neighbors_old,jn_neighbors_new)
	initial_surfaces = self.get_junction_surface_count(total_sites)

	self.update_junction_super_lattice(oldIdx,'remove')
	self.update_junction_neighbors(jn_neighbors_old,'remove')
	
	#lost_surfaces = self.update_junction_neighbors(oldIdx,'remove')
	
	self.update_junction_super_lattice(new_pos,'add')
	self.update_junction_neighbors(jn_neighbors_new,'add')
	
	final_surfaces = self.get_junction_surface_count(total_sites)
	
	#jn_final_surfaces_old = self.get_junction_surface_count(jn_neighbors_old)
	#jn_final_surfaces_new = self.get_junction_surface_count(jn_neighbors_new)
	
	self.junction_surface_count += final_surfaces - initial_surfaces

	if boxno==(len(self.all_chains[chainNo])-1):
	    if chainNo not in self.active_chains and chainNo not in self.early_terminations:
		other_chain = self.terminated_pairs[chainNo]
		self.all_chains[other_chain][-1] = new_pos

    def update_super_lattice(self,idx,update_type):
	if update_type == 'add':
	    #self.blocked_matrix[idx[0],idx[1],idx[2]] = 1
	    self.blocked_site_count += 1
	elif update_type == 'remove':
	    #self.blocked_matrix[idx[0],idx[1],idx[2]] = 0
	    self.blocked_site_count -= 1

	IdxSet = IndexUtils.get_super_lattice_pts_periodic(idx)

	for idx in IdxSet:
	    if update_type == 'add':
		self.super_lattice[idx[0],idx[1],idx[2]] = 1

            if update_type == 'remove':
                self.super_lattice[idx[0],idx[1],idx[2]] = 0

    def check_super_lattice(self,idx):        
        IdxSet = IndexUtils.get_super_lattice_pts(idx)

        count = 0
        
        for idxLoop in IdxSet:
            if self.super_lattice[idxLoop[0],idxLoop[1],idxLoop[2]] == 0:
                count += 1
                
        if count == 8:
            check = 'empty'
        else:
            check = 'not empty'
            
        return count

    def check_super_lattice_thermal(self,oldidx, newidx): 
        oldIdxSet = IndexUtils.get_super_lattice_pts_periodic(oldidx)
        newIdxSet = IndexUtils.get_super_lattice_pts_periodic(newidx)

        count = 0

        for idxLoop in newIdxSet:
            if idxLoop not in oldIdxSet:
                if self.super_lattice[idxLoop[0],idxLoop[1],idxLoop[2]] == 0:
                    count += 1

        if count == 4:
            check = 'empty'
        else:
            check = 'not empty'

	#print count, check

        return check
    
    def get_associated_energy(self,idx):
	total_energy = 0.0
	if total_energy>(1.001*self.fluctuator.maximum_energies[3]):
	    print 'Error: ', total_energy
	    sys.exit()

	chain_no = self.location_chain[idx[0],idx[1],idx[2]]
	box_no = self.all_chains[chain_no].index(idx)
	    
	# Find the location of the previous two nodes
	try:
	    prev_1 = self.all_chains[chain_no][box_no-1]
	except IndexError:
	    prev_1 = None

	try:
	    prev_2 = self.all_chains[chain_no][box_no-2]
	except IndexError:
	    prev_2 = None
	    
	if prev_1!=None:
	    total_energy += 0.5*self.fluctuator.get_bond_energy(idx,prev_1)
	if prev_2!=None:
	    angle = self.fluctuator.get_angle(prev_2,prev_1,idx)
	    total_energy += self.fluctuator.get_angle_energy(angle)
	    
	# Find the location of the next two nodes    
	try:
	    next_1 = self.all_chains[chain_no][box_no+1]
	except IndexError:
	    next_1 = None
	
        try:
	    next_2 = self.all_chains[chain_no][box_no+2]
	except IndexError:
	    next_2 = None
	    
	# Decide node type
	if prev_1==idx:
	    node_type = 'isolated'
	else:
	    node_type = 'regular'
	
    	if next_1!=None:
	    total_energy += 0.5*self.fluctuator.get_bond_energy(idx,next_1)
	if next_2!=None:
	    angle = self.fluctuator.get_angle(next_2,next_1,idx)
	    total_energy += self.fluctuator.get_angle_energy(angle)
	    
	if self.junction_chain[idx[0],idx[1],idx[2]]>0:
	    jn_chain = self.junction_chain[idx[0],idx[1],idx[2]]
	    
	    if jn_chain != chain_no:
		second_box = self.all_chains[jn_chain].index(idx)
	    else:
		for k, try_idx in enumerate(self.all_chains[jn_chain]):
		    if try_idx==idx and k!=box_no:
			second_box = k

	    # Find the location of the previous two nodes
	    try:
		prev_1 = self.all_chains[jn_chain][second_box-1]
	    except IndexError:
		prev_1 = None
	
	    try:
		prev_2 = self.all_chains[jn_chain][second_box-2]
	    except IndexError:
		prev_2 = None

	    if prev_1!=None:
		total_energy += 0.5*self.fluctuator.get_bond_energy(idx,prev_1)
	    if prev_2!=None:
		angle = self.fluctuator.get_angle(prev_2,prev_1,idx)
		total_energy += self.fluctuator.get_angle_energy(angle)

	    # Find the location of the next two nodes    
	    try:
		next_1 = self.all_chains[jn_chain][second_box+1]
	    except IndexError:
		next_1 = None
	    
	    try:
		next_2 = self.all_chains[jn_chain][second_box+2]
	    except IndexError:
		next_2 = None

	    if next_1!=None:
		total_energy += 0.5*self.fluctuator.get_bond_energy(idx,next_1)
	    if next_2!=None:
		angle = self.fluctuator.get_angle(next_2,next_1,idx)
		total_energy += self.fluctuator.get_angle_energy(angle)

	return total_energy, node_type
		
    def get_thermal_propensity(self):
	total_chains = len(self.all_chains)
        usedMonomers = self.TotalMonomers - self.transientMonomers

        try:
            total_propensity = self.global_propensities['fluctuate']
        except KeyError:
            total_propensity = 0.0
	    
	#print 'Effected nodes: ', len(self.effected_all), self.total_nodes

	if total_chains!=0 and self.total_nodes!=0:
	#    for chain in self.all_chains:
	#	for idx in self.all_chains[chain]:
	    for idx in self.effected_all:
		    if self.location_chain[idx[0],idx[1],idx[2]]==0:
			print 'Wrong'
			sys.exit()

		    i = self.get_lidx_from_cubidx(idx)
		    
		    mo_type = self.location_monomers[idx[0],idx[1],idx[2]]
	
		    #self.associated_energy[idx[0],idx[1],idx[2]], node_type = self.get_associated_energy(idx)
		    
		    node_type = -1
		    if self.junction_chain[idx[0],idx[1],idx[2]]==0 and self.Nr[idx[0],idx[1],idx[2]]>=1:
			# max_displacement_energy_index
			tr_idx = 0
			node_type = 0
		    elif self.Nr[idx[0],idx[1],idx[2]]>=1:
			# max_displacement_energy_index
			tr_idx = 0
			node_type = 2
		    elif self.junction_chain[idx[0],idx[1],idx[2]]>0:
			# max_displacement_energy_index
			tr_idx = 3
			node_type = 3
		    else:
			# max_displacement_energy_index
			tr_idx = 1
			node_type = 1

		    if node_type!=-1:
			#energy_index = int(self.associated_energy[idx[0],idx[1],idx[2]]/self.fluctuator.energy_quanta)
			#jump_probability = self.fluctuator.escape_probabilities[node_type][energy_index]
			jump_probability, move_dict, move_count, move_prob_dict, max_prob, bond_update = self.compute_jump_prob(idx)
			self.fluctuation_quanta[int(jump_probability)] += 1
		    else:
			jump_probability = 1.0
			
		    if self.thermal_propensities[i]>0.0:
			total_propensity += -self.thermal_propensities[i]
			
			if self.thermal_propensities[i]>=1.0:
			    self.highly_movable_nodes.remove(idx)
			    self.total_highly_movable_nodes += -1
			    self.highly_movable_propensities += -self.thermal_propensities[i]
			    self.high_propensities.remove(self.thermal_propensities[i])
			
			self.movable_nodes.remove(idx)
			self.thermal_propensities[i] = 0

			for k in xrange(0,6):
			    self.movable_dict[i,k] = 0
			    self.moving_probs[i,k] = 0
			    self.bond_motion_update[i,k] = 0
			    
			self.max_moving_probs[i] = 0
			self.total_movable_nodes += -1

		    if move_count>0:
			self.thermal_propensities[i] = jump_probability
			
			for k in xrange(0,6):
			    self.movable_dict[i,k] = move_dict[k]
			    self.moving_probs[i,k] = move_prob_dict[k]
			    self.bond_motion_update[i,k] = bond_update[k]
			    
			self.max_moving_probs[i] = max_prob
			self.movable_nodes.append(idx)
			self.total_movable_nodes += 1
			
			if jump_probability>=1.0:
			    self.highly_movable_nodes.append(idx)
			    self.total_highly_movable_nodes += 1
			    self.highly_movable_propensities += self.thermal_propensities[i]
			    self.high_propensities.append(self.thermal_propensities[i])
			
			total_propensity += self.thermal_propensities[i]
		
	self.global_propensities['fluctuate'] = total_propensity

	return total_propensity

    def fluctuation_statistics(self):
	self.fluctuation_samples.append(self.global_propensities['fluctuate'])
	self.fluctuation_times.append(self.current_time)
	
	#self.fluctuation_sum += self.global_propensities['fluctuate']
	#self.fluctuation_square_sum += self.global_propensities['fluctuate']**2

    def get_radical_termination_propensity(self):
	if len(self.active_chains)<=1:
	    self.global_propensities['radical termination'] = 0.0
	else:
	    self.global_propensities['radical termination'] = 0.0
	    avg_chain_length = 0.0
	    chain_count = 0
	    
	    for chain_no in self.active_chains:
		if len(self.chain_connectivity[chain_no])==0:
		    avg_chain_length += len(self.all_chains[chain_no])
		    chain_count += 1
		elif self.check_other_connection(chain_no)=='not connected':
		    avg_chain_length += len(self.all_chains[chain_no])
		    chain_count += 1
		    
	    if chain_count>1:
		avg_chain_length *= (1.0/float(chain_count))
		
		self.oligomer_eps_self = avg_chain_length*self.all_monomers.mean_eps_self
		
		self.d_factor_oligomer = self.diffusion.compute_effective_diffusivity(self.oligomer_eps_self)
		
		if math.isnan(self.d_factor_oligomer)==True:
		    propensity = 0.0
		else:
		    combined_diffusivity = 2.0*self.fluctuator.mean_diff*self.d_factor_oligomer/float(avg_chain_length)
		    sigma12 = 4.0*self.lattice_h*(avg_chain_length**0.6)
		    
		    reduced_mass = 0.5*self.molWt*avg_chain_length*(1E-26)/6.022
		    vel12 = (1E10)*np.sqrt((8.0*self.kB*self.temp)/(math.pi*reduced_mass))
		    
		    propensity = 4*math.pi*combined_diffusivity*vel12*(sigma12**2)
		    propensity *= 1.0/(self.sysVol*(4.0*combined_diffusivity + sigma12*vel12))
		    propensity *= 0.5*len(self.active_chains)*(len(self.active_chains)-1)

		self.global_propensities['radical termination'] = propensity

    def get_network_diffusivity(self):
	self.update_diffusivities()
	if self.total_bonds>0:
	    #self.node_propensity = 26.0/max(self.surface_site_count,1.0)
	    self.node_propensity = 1.0 #- float(self.total_bonds)/float(self.total_sites)
	else:
	    self.node_propensity = 0.0
	
    def update_diffusivities(self):
	#surface_count = sum(sum(sum(self.junc_neighbors)))
	#
	#if surface_count!=self.junction_surface_count:
	#    print 'Surface count error: ', surface_count, self.junction_surface_count
	#    sys.exit()
	
	#self.diffusion.compute_interface_probabilities(self.totalJunctionPoints,self.junction_surface_count,self.total_nodes,self.surface_site_count)
	
	self.diffusion.compute_interface_probabilities(self.realJunctionPoints,self.junction_surface_count,self.total_nodes,self.surface_site_count)
	
	#self.d_factor_monomer = self.diffusion.compute_effective_diffusivity(self.eps_self)
	
	for mo_type in xrange(0,self.all_monomers.no_monomers):
	    self.diffusion.set_initial_guess(self.all_monomers.d_factor_monomers[mo_type])
	    self.all_monomers.d_factor_monomers[mo_type] = self.diffusion.compute_effective_diffusivity(self.all_monomers.eps_self[mo_type])

	self.all_monomers.construct_prob_dist()

	self.diffusion.set_initial_guess(self.mean_diffusion_factor)
	self.mean_diffusion_factor = self.diffusion.compute_effective_diffusivity(self.all_monomers.mean_eps_self)
	
	if self.total_nodes>0:
	    self.lattice_propensity = self.all_monomers.mean_d_factors*self.global_propensities['fluctuate']/(6.0*float(self.total_nodes))
	else:
	    self.lattice_propensity = 1.0

	if self.initiation_type=='bimolecular':
	    try:
		self.diffusion.set_initial_guess(self.d_factor_initiator)
	    except AttributeError:
		self.diffusion.set_initial_guess(1.0)
	    self.d_factor_initiator = self.diffusion.compute_effective_diffusivity(self.eps_I)
	    
	    try:
		self.diffusion.set_initial_guess(self.d_factor_coinitiator)
	    except AttributeError:
		self.diffusion.set_initial_guess(1.0)

	    self.d_factor_coinitiator = self.diffusion.compute_effective_diffusivity(self.eps_co_I)
	    #print 'Diffusivities: ', self.d_factor_initiator, self.d_factor_coinitiator

    def get_lattice_diffusivity(self):
	# Computes the scaling factor for the diffusivity for any lattice site, network or not
	# Free node diffusivities
	node_diffusivity = 6.0*self.fluctuator.mean_propensity*self.diffusion_factor
	free_sites = self.total_sites - self.total_nodes
	# Free site contributions
	free_sites_total = free_sites*node_diffusivity
	
	# Total contribution for all network nodes
	network_sites_total = self.network_propensity*self.global_propensities['fluctuate']
	
	# Effective lattice propensity
	self.global_lattice_propensity = (free_sites_total+network_sites_total)/(6.0*self.fluctuator.mean_propensity*self.total_sites)
	
    def get_event_type_and_step(self,leap_status):
	xi = rand.uniform(0.0,(1.0-1e-3))

	self.A_total = 0.0
	self.non_fluctuation_propensity = 0.0
	self.fluctuation_propensity_multiplier = 0.0
	max_prop = 0.0
	
	#if self.total_nodes<25000:
	#    leap_status = 'leap'
	    
	#if len(self.active_chains)==1:# and len(self.all_chains[1])<1000:
	#    self.global_propensities['initiate'] = 0
	
	if self.max_length!=None and bool(self.all_chains)==True:
	    self.global_propensities['initiate'] = 0.0
	    if len(self.all_chains[1])<self.max_length:
		flag = -1
	    elif len(self.all_chains[1])==self.max_length:
		self.global_propensities['polymerize'] = 0.0
		flag = 1
	    else:
		flag = 1
	else:
	    flag = 1
	    
	if leap_status=='leap':
	    if int(sum(self.global_propensities.values()))==int(self.global_propensities['fluctuate']):
	        leap_status = 'no leap'
	    
	for event_type in self.global_propensities.keys():
	    if math.isnan(self.global_propensities[event_type])==True:
		self.global_propensities[event_type] = 0.0

	    if event_type=='fluctuate':
		if leap_status=='no leap' and flag==1:
		    self.fluctuation_propensity_multiplier = (self.all_monomers.mean_d_factors*self.fluctuator.mean_lattice_propensity)
		    self.A_total += self.fluctuation_propensity_multiplier*self.global_propensities[event_type]
		    #print event_type,': ',(self.network_propensity*self.global_propensities[event_type])
		    max_prop = max(max_prop,self.fluctuation_propensity_multiplier*self.global_propensities[event_type])
	    else:
		self.A_total += self.global_propensities[event_type]
		self.non_fluctuation_propensity += self.global_propensities[event_type]
		max_prop = max(max_prop,self.global_propensities[event_type])
		#print event_type,': ',(self.global_propensities[event_type])

        step_type, step_size = None, None

	#if self.fluctuation_transition_status==False and self.total_nodes>100:
	#    check1 = self.non_fluctuation_propensity/(self.all_monomers.mean_d_factors*self.fluctuator.mean_lattice_propensity)
	#    check2 = 1.0
	#    
	#    check = check1 - check2
	#    
	#    if check>0:
	#	print 'Fluctuation status changed ', self.DoC, self.non_fluctuation_propensity, self.all_monomers.mean_d_factors, self.fluctuator.mean_lattice_propensity
	#	self.fluctuation_transition_status = True
		
	#if self.diffusion.phi[0]<1.0/3.0:
	#    self.equilibrium_tolerance *= 0.1
	    #self.fluctuation_transition_status = True
	    
	if self.A_total > 0.0:
            for key in self.global_propensities.keys():
		if key=='fluctuate':
		    if leap_status=='no leap':
			this_prob = self.fluctuation_propensity_multiplier*self.global_propensities[key]/max_prop
		else:
		    this_prob = self.global_propensities[key]/max_prop
		    
		if this_prob>=xi:
		    step_type = key
		    break

	    step_size = -np.log(rand.uniform(1E-16,1.0))/self.A_total

#        prob_sum = 0.0
#
#        if self.A_total > 0.0:
#            for key in self.global_propensities.keys():
#		if key=='fluctuate':
#		    if leap_status=='no leap':
#			this_prob = self.fluctuation_propensity_multiplier*self.global_propensities[key]/self.A_total
#		else:
#		    this_prob = self.global_propensities[key]/self.A_total
#		
#		if prob_sum<xi and xi<=(prob_sum+this_prob):
#		    step_type = key
#
#	        prob_sum += this_prob
#
#            step_size = -np.log(rand.uniform(1E-16,1.0))/self.A_total

	#if leap_status=='leap':
	#    if step_type!='polymerize':
	#	return_value = 'leap'
		
	if step_type==None:
	    print sum(self.global_propensities.values()), self.global_propensities['fluctuate']
	    print 'Propensities wrong!\n', leap_status, self.A_total, xi, max_prop, flag, self.step_no
	    print self.global_propensities, max_prop, self.fluctuation_propensity_multiplier, self.diffusion.phi
	    print leap_status
	    
	if leap_status=='leap':
	    leap_status = 'no leap'
	    
	if leap_status=='equilibriate':
	    step_type = 'fluctuate'
	    step_size = 0.0
	
	#if len(self.active_chains)==1 and len(self.all_chains[1])<1000:
	#    step_type = 'polymerize'
	
	#if len(self.all_chains)==1 and len(self.all_chains[1])==1000:
	#    self.global_propensities['initiate'] = 0.0
	#    self.global_propensities['polymerize'] = 0.0
	#    step_type = 'fluctuate'

        return step_type, step_size, xi, leap_status

    def select_thermal_unit_fast(self):
	diff = -1.0
	
	real_wt = 0.0
	
	xi_1 = rand.uniform(0.0,1.0)
	
	check = self.highly_movable_propensities/self.global_propensities['fluctuate'] - xi_1
	
	attempts = 0
	
	for i in range(6,-1,-1):
	    if self.fluctuation_quanta[i]>0:
		self.max_fluc_wt = float(i)
		break
	
	if self.total_highly_movable_nodes>0 and check>0.0:
	    while diff<0.0:
		attempts += 1
		xi = rand.uniform(0.0,(1.0-1e-3))
		# Uses rejection sampling to draw a node to fluctuate
		
		selected_node = self.highly_movable_nodes[rand.randint(0,self.total_highly_movable_nodes-1)]
		#selected_node = self.movable_nodes[node_set[attempts%self.total_movable_nodes]]
		
	    #    if self.movability[selected_node[0],selected_node[1],selected_node[2]]==0:
	    #	print 'Error', self.movable_dict[self.get_lidx_from_cubidx(selected_node)],self.movability[selected_node[0],selected_node[1],selected_node[2]]
	    #	sys.exit()
	    #    else:		
		i = self.get_lidx_from_cubidx(selected_node)

		real_wt = self.thermal_propensities[i]
		
		diff = real_wt/self.max_fluc_wt - xi
	else:
	    while diff<0.0:
		attempts += 1
		
		xi = rand.uniform(0.0,(1.0-1e-3))
		
		selected_node = self.movable_nodes[rand.randint(0,self.total_movable_nodes-1)]
		i = self.get_lidx_from_cubidx(selected_node)

		#real_wt = self.all_thermal_propensities[i]
		real_wt = self.thermal_propensities[i]
		
		if real_wt<1.0:
		    diff = real_wt/self.fluctuator.max_transition_factor - xi
		else:
		    diff = -1.0

	#print 'Succesful after ', attempts, 'attempts'
	
	self.max_fluc_trials = max(self.max_fluc_trials,attempts)
	
	return selected_node	    

    def get_lidx_from_cubidx(self,box_idx):
        l_idx = box_idx[0]+box_idx[1]*self.lattice_size+box_idx[2]*(self.lattice_size**2)
        
        return l_idx          
    
    def influenced_nodes(self,updated_idx,step_type):
	if step_type=='fluctuate':
	    all_idxs = IndexUtils.get_influenced_monomers_thermal(updated_idx,self.last_move,self.location_chain)
        else:
	    all_idxs = IndexUtils.get_influenced_monomers(updated_idx,self.location_chain)
	#all_idxs = IndexUtils.get_influenced_monomers(updated_idx)
	
        effected_idxs = list(all_idxs)
	#movable_idxs = []
	
	self.effected_mov_count = 0

#            if self.movability[idx[0],idx[1],idx[2]] > 0:
#		self.effected_mov_count += 1
#                movable_idxs.append(idx)

	#print 'Nodes to update: ', self.effected_mov_count
	
	loc_chain = int(self.location_chain[updated_idx[0],updated_idx[1],updated_idx[2]])
	box_no = int(self.location_boxno[updated_idx[0],updated_idx[1],updated_idx[2]])
	idxs = self.get_connections(loc_chain,box_no)
	
	for t_idx in idxs:
	    if t_idx!=None:
		if self.location_chain[t_idx[0],t_idx[1],t_idx[2]]>0:
		    effected_idxs.append(t_idx)
		
	if self.junction_chain[updated_idx[0],updated_idx[1],updated_idx[2]]>0:
	    loc_chain = int(self.junction_chain[updated_idx[0],updated_idx[1],updated_idx[2]])
	    box_no = int(self.junction_boxno[updated_idx[0],updated_idx[1],updated_idx[2]])
	    idxs = self.get_connections(loc_chain,box_no)
	    
	    for t_idx in idxs:
		if t_idx!=None:
		    if self.location_chain[t_idx[0],t_idx[1],t_idx[2]]>0:
			effected_idxs.append(t_idx)

        return effected_idxs #, movable_idxs
    
    def get_connections(self,chain_no,box_no):
	boxlow_1, boxup_1 = None, None
	
	if box_no == 0:
	    boxlow_1 = None
	else:
	    try:
		boxlow_1 = self.all_chains[chain_no][box_no-1]
	    except IndexError:
		boxlow_1 = None

	if box_no == len(self.all_chains[chain_no])-1:
	    if chain_no not in self.active_chains and chain_no not in self.early_terminations:
		other_chain = self.terminated_pairs[chain_no]
		
		try:
		    boxup_1 = self.all_chains[other_chain][-2]
		except IndexError:
		    boxup_1 = None
	else:
	    boxup_1 = self.all_chains[chain_no][box_no+1]

	return boxlow_1, boxup_1
    
    def check_neighbors(self,neighbor_set):
	filled_sites = []
	
	for neigh_idx in neighbor_set:
	    if self.check_super_lattice(neigh_idx)<8:
		# A value of 1 for the occupied site
		filled_sites.append(1)
	    else:
		# A value of 0 for the unoccupied site
		filled_sites.append(0)

	return filled_sites

    def update_neighbor_count(self,site_values,update_type,site_type):
	update = 26 - sum(site_values)
	if update_type=='add':
	    self.neighbor_site_count += update
	    
	    if site_type=='junction':
		self.junction_neighbor_count += update
	elif update_type=='remove':
	    self.neighbor_site_count -= update
	    
	    if site_type=='regular':
		self.junction_neighbor_count -= update

    def local_surface_count(self,neighbor_set):
	surface_site_update = 0
	
	for neigh_idx in neighbor_set:
	    site_value = self.check_super_lattice(neigh_idx)
	    #if site_value!=0 and site_value!=8:
	    if site_value>=4 and site_value<8:
		surface_site_update += 1
		
	return surface_site_update
    
    def compute_pair_correlation(self,tag):
	# Fractional occupancy of the lattice by the nodes
	self.rho_rs = self.total_nodes/self.total_sites
	
	# Radius limit for computing the pair correlation
	self.R_lim = 30
	
	self.node_correlations = np.zeros(shape=(self.R_lim+1))
	self.jn_correlations = np.zeros(shape=(self.R_lim+1))
	
	# Site density
	rho_rs = float(self.total_nodes)/float(self.total_sites)

	# Junction density
	rho_jn = float(self.totalJunctionPoints)/float(self.total_sites)
	id_max = int(self.lattice_size)
	
	for x_idx in range(0,id_max):
	    for y_idx in range(0,id_max):
		for z_idx in range(0,id_max):
		    if self.location_chain[x_idx,y_idx,z_idx]!=0:
			this_g_r, this_j_r = self.correlation_fn([x_idx,y_idx,z_idx])
			for r_val in range(0,self.R_lim):
			    self.node_correlations[r_val] += this_g_r[r_val]/rho_rs
			    self.jn_correlations[r_val] += this_j_r[r_val]/rho_jn
	
	pair_file = open('pc_'+tag+'_'+str(self.proc_no)+'.csv','w')
	
	print >> pair_file, 'Radial distance,Node correlation,Junction correlation'
	
	for r_val in range(0,self.R_lim):
	    self.node_correlations[r_val] *= (1.0/self.total_nodes)
	    self.jn_correlations[r_val] *= (1.0/self.totalJunctionPoints)
	    
	    print >> pair_file, r_val, ',', self.node_correlations[r_val], ',', self.jn_correlations[r_val]
	    
	pair_file.close()

    def correlation_fn(self,idx):
	id_max = int(self.lattice_size)
	idx_lim = self.R_lim
	idx2 = [0,0,0]
	g_r = np.zeros(shape=(self.R_lim+1))
	j_r = np.zeros(shape=(self.R_lim+1))
	
	N_ls_pop = np.zeros(shape=(self.R_lim+1))
	N_rs_pop = np.zeros(shape=(self.R_lim+1))
	N_js_pop = np.zeros(shape=(self.R_lim+1))
	
	for x_idx in range(idx[0]-idx_lim,idx[0]+idx_lim):
	    if x_idx>=self.lattice_size:
                idx2[0] = x_idx - id_max
	    elif x_idx<0:
		idx2[0] = id_max + x_idx
	    else:
		idx2[0] = x_idx

	    for y_idx in range(idx[1]-idx_lim,idx[1]+idx_lim):
		if y_idx>=self.lattice_size:
		    idx2[1] = y_idx - id_max
		elif y_idx<0:
		    idx2[1] = id_max + y_idx
		else:
		    idx2[1] = y_idx
		    
		for z_idx in range(idx[2]-idx_lim,idx[2]+idx_lim):
		    if z_idx>=self.lattice_size:
		        idx2[2] = z_idx - id_max
	            elif z_idx<0:
		        idx2[2] = id_max + z_idx
		    else:
		        idx2[2] = z_idx

		    dist = self.correlation_distance(idx,[x_idx,y_idx,z_idx])
		    r_idx = int(dist)

		    if r_idx<(idx_lim+1):
			N_ls_pop[r_idx] += 1
			
			if self.location_chain[idx2[0],idx2[1],idx2[2]]>0:
			    N_rs_pop[r_idx] += 1

			if self.junction_chain[idx2[0],idx2[1],idx2[2]]>0:
			    N_js_pop[r_idx] += 1
			    
	g_r = np.zeros(shape=(self.R_lim+1))
	j_r = np.zeros(shape=(self.R_lim+1))
	
	for r_val in range(0,self.R_lim+1):
	    if N_ls_pop[r_val] != 0:
		g_r[r_val] = float(N_rs_pop[r_val])/float(N_ls_pop[r_val])
		j_r[r_val] = float(N_js_pop[r_val])/float(N_ls_pop[r_val])
	    else:
		g_r[r_val] = 1.0
		j_r[r_val] = 1.0

	return g_r, j_r

    def correlation_distance(self,idx1,idx2):
	total = 0.0

	for i in range(0,3):
	    # Shortest distance consider periodicity
	    short1 = abs(idx1[i] - idx2[i])
	    
	    total += short1**2

	dist = total

	return dist

    def compute_MC_pair_correlation(self,tag):
	id_max = int(self.lattice_size)
	
	# Monte Carlo sampling size
	#mc_sample_size = self.total_sites*int(math.log10(self.total_sites))
	mc_sample_size = (len(self.rad_list))*1000
	#mc_sample_size = 100000
	
	sample_count = 0
	n_site_count, j_site_count = 0,0
	#node_pair_count = 0
	#junction_pair_count = 0
	
	node_pair_count, junction_pair_count = mc_sample_size, mc_sample_size
	
	N_ls_pop_n = np.zeros(shape=(len(self.rad_list)))
	N_ls_pop_j = np.zeros(shape=(len(self.rad_list)))
	
	N_rs_pop = np.zeros(shape=(len(self.rad_list)))
	N_js_pop = np.zeros(shape=(len(self.rad_list)))
	
	# Pair correlation functions
	g_r_n = np.zeros(shape=(len(self.rad_list)))
	g_r_j = np.zeros(shape=(len(self.rad_list)))
	
	g_r_n_var = np.zeros(shape=(len(self.rad_list)))
	g_r_j_var = np.zeros(shape=(len(self.rad_list)))
	
	# Crosslink data structures
	junctions_sampled = 0
	sampled_coordinates = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	
	for chain_no in self.all_chains.keys():
	    for idset in self.all_chains[chain_no]:
		if self.junction_chain[idset[0],idset[1],idset[2]]<=chain_no and sampled_coordinates[idset[0],idset[1],idset[2]]==0:
		    sampled_coordinates[idset[0],idset[1],idset[2]] = -1
		    
		    if self.junction_chain[idset[0],idset[1],idset[2]]==0:
			ls_pop_local, rs_local = self.compute_sample_rdfs(idset)
		    else:
			junctions_sampled += 1
			ls_pop_local, rs_local, ls_pop_j_local, rs_j_local = self.compute_sample_rdfs(idset)

		    for j in xrange(0,len(self.rad_list)):
			N_ls_pop_n[j] += ls_pop_local[j]
			N_rs_pop[j] += rs_local[j]/float(ls_pop_local[j])
			
			if ls_pop_local[j]>0:
			    local_rho = rs_local[j]/float(ls_pop_local[j])
			else:
			    local_rho = 0
			
			g_r_n_var[j] += local_rho**2
			
			if self.junction_chain[idset[0],idset[1],idset[2]]>0:
			    N_ls_pop_j[j] += ls_pop_j_local[j]
			    N_js_pop[j] += rs_j_local[j]/float(ls_pop_j_local[j])
			    
			    if ls_pop_j_local[j]>0:
				local_j_rho = rs_j_local[j]/float(ls_pop_j_local[j])
			    else:
				local_j_rho = 0
				
			    g_r_j_var[j] += local_j_rho**2
	
	#for i in xrange(0,self.total_nodes):
	#    chain_no = rand.randint(1,len(self.all_chains))
	#    node_no = rand.randint(0,len(self.all_chains[chain_no])-1)
	#    idset = self.all_chains[chain_no][node_no]
	#    
	#    ls_pop_local, rs_local = self.compute_sample_rdfs(idset)
	#    
	#    for j in xrange(0,len(self.rad_list)):
	#	N_ls_pop_n[j] += ls_pop_local[j]
	#	N_rs_pop[j] += rs_local[j]
	
	#while node_pair_count>=0 or junction_pair_count>=0:
	#    print node_pair_count, junction_pair_count
	#    # Pick the first random node from the network
	#    chain_no = rand.randint(1,len(self.all_chains))
	#    node_no = rand.randint(0,len(self.all_chains[chain_no])-1)
	#    idset1 = self.all_chains[chain_no][node_no]
	#    sample_count += 1
	#    
	#    idset2, dist, c_idx1 = self.select_random_pair_loc(idset1)
	#    
	#    # Pick junction node
	#    j_id = self.junction_list[rand.randint(0,self.totalJunctionPoints-1)]
	#    j_id2, dist2, c_idx2 = self.select_random_pair_loc(j_id)
	#    
	#    if dist>=4.0 and dist<(self.R_lim+1)**2:
	#        r_idx = self.rad_list.index(dist)
	#	
	#	if N_ls_pop_n[r_idx]<100:
	#	    N_ls_pop_n[r_idx] += 1
	#	    n_site_count += 1
	#
	#	    if self.location_chain[idset1[0],idset1[1],idset1[2]]>0 and self.location_chain[idset2[0],idset2[1],idset2[2]]>0:
	#		N_rs_pop[r_idx] += 1
	#		node_pair_count += -1
	#
	#    if dist>=4.0 and dist2<(self.R_lim+1)**2:
	#        r_idx = self.rad_list.index(dist)
	#	
	#	if N_ls_pop_j[r_idx]<100:
	#	    N_ls_pop_j[r_idx] += 1
	#	    j_site_count += -1
	#	     
	#	    if self.junction_chain[j_id[0],j_id[1],j_id[2]]>0 and self.junction_chain[j_id2[0],j_id2[1],j_id2[2]]>0:
	#		N_js_pop[r_idx] += 1
	#		junction_pair_count += -1

	# Site density
	rho_rs_true = float(self.total_nodes)/float(self.total_sites)
	#rho_rs = float(node_pair_count)/float(n_site_count)
	# Junction density
	rho_jn_true = float(self.totalJunctionPoints)/float(self.total_sites)
	#rho_jn = float(junction_pair_count)/float(j_site_count)
	
	#print 'rhos: ', rho_rs, rho_jn, sample_count, n_site_count, j_site_count, node_pair_count, junction_pair_count
	print 'true rhos : ', rho_rs_true, rho_jn_true
	print 'junctions: ', self.totalJunctionPoints, junctions_sampled
	sys.stdout.flush()

	pair_file = open('pc_'+tag+'_'+str(self.proc_no)+'.csv','w')
	pair_var_file = open('pcvar_'+tag+'_'+str(self.proc_no)+'.csv','w')
	#eriod_file = open('periodicity'+str(self.proc_no)+'.csv','w')
		
	print >> pair_file, 'Radial distance,Node correlation,Junction correlation'
	
	#print >> period_file, 'Lattice distance,Node X,Node Y, Node Z, Junction X, Junction Y, Junction Z'
		
	node_correlations = np.zeros(shape=(len(self.rad_list)))
	jn_correlations = np.zeros(shape=(len(self.rad_list)))
	
	for r_idx in range(0,len(self.rad_list)):
	    if N_ls_pop_n[r_idx]>0:
		node_correlations[r_idx] = float(N_rs_pop[r_idx])/(float(self.total_nodes))
	    else:
		node_correlations[r_idx] = 1.0

	    if N_ls_pop_j[r_idx]>0:
		jn_correlations[r_idx] = float(N_js_pop[r_idx])/(float(self.totalJunctionPoints))
	    else:
		jn_correlations[r_idx] = 1.0
		
	    g_r_n_var[r_idx] *= 1.0/float(self.total_nodes)
	    g_r_n_var[r_idx] += -node_correlations[r_idx]**2
	    g_r_n_var[r_idx] = g_r_n_var[r_idx]/(node_correlations[r_idx]**2)
	    node_correlations[r_idx] *= 1.0/rho_rs_true
	    
	    g_r_j_var[r_idx] *= 1.0/float(junctions_sampled)
	    g_r_j_var[r_idx] += -jn_correlations[r_idx]**2
	    g_r_j_var[r_idx] = g_r_j_var[r_idx]/(jn_correlations[r_idx]**2)
	    jn_correlations[r_idx] *= 1.0/rho_jn_true

	    print >> pair_file, math.sqrt(self.rad_list[r_idx]), ',', node_correlations[r_idx], ',', jn_correlations[r_idx]
	    print >> pair_var_file, math.sqrt(self.rad_list[r_idx]), ',', g_r_n_var[r_idx], ',', g_r_j_var[r_idx]
	    
	    #print >> period_file, math.sqrt(self.rad_list[r_idx]), ',', x_r_p, ',', y_r_p, ',', z_r_p, ',', x_j_p, ',', y_j_p, ',', z_j_p
	    
	pair_file.close()
	pair_var_file.close()
	#period_file.close()
	
	#self.compute_oroganov_order_parameters(node_correlations,jn_correlations)
	
    def compute_sample_rdfs(self,idx):
	N_ls_pop_n = np.zeros(shape=(len(self.rad_list)))
	
	N_rs_pop = np.zeros(shape=(len(self.rad_list)))
	
	if self.junction_chain[idx[0],idx[1],idx[2]]>0:
	    N_ls_pop_j = np.zeros(shape=(len(self.rad_list)))
	    N_js_pop = np.zeros(shape=(len(self.rad_list)))
	
	for x in range(-self.R_lim,self.R_lim+1):
	    for y in range(-self.R_lim,self.R_lim+1):
		for z in range(-self.R_lim,self.R_lim+1):
		    trial_x = [idx[0]+x,idx[1]+y,idx[2]+z]
		    
		    for ix in xrange(0,3):
			if trial_x[ix] >= self.lattice_size:
			    trial_x[ix] = trial_x[ix] - self.lattice_size
			elif trial_x[ix] < 0:
			    trial_x[ix] = self.lattice_size + trial_x[ix]
			    
		    dist = self.correlation_distance(idx,trial_x)
		    
		    if dist>=4.0 and dist<(self.R_lim+1)**2:
			r_idx = self.rad_list.index(dist)
			N_ls_pop_n[r_idx] += 1
			
			if self.location_chain[trial_x[0],trial_x[1],trial_x[2]]>0:
			    N_rs_pop[r_idx] += 1
			    
			if self.junction_chain[idx[0],idx[1],idx[2]]>0:
			    N_ls_pop_j[r_idx] += 1
			    if self.junction_chain[trial_x[0],trial_x[1],trial_x[2]]>0:
				N_js_pop[r_idx] += 1

	if self.junction_chain[idx[0],idx[1],idx[2]]>0:
	    return N_ls_pop_n, N_rs_pop, N_ls_pop_j, N_js_pop
	else:
	    return N_ls_pop_n, N_rs_pop
	
    def create_correlation_radius_list(self):
	# Create a dict with the squared distance as keys and the number of sites
	self.rad_list = []
	#self.rad_list = np.linspace(2,self.R_lim,self.R_lim-1)

	#self.total_site_list = np.zeros(shape=(len(self.rad_list)))
			    
	for x in range(0,self.R_lim+1):
	    for y in range(0,self.R_lim+1):
		for z in range(0,self.R_lim+1):
		    dist2 = self.correlation_distance([0,0,0],[x,y,z])
		    if dist2>=4.0 and dist2<(self.R_lim+1)**2:
			if dist2 not in self.rad_list:
			    self.rad_list.append(dist2)
	
	self.rad_list.sort()
		
	#total_ring_sites = sum(self.total_site_list)
	#
	#for r_idx in range(0,len(self.rad_list)):
	#    self.total_site_list[r_idx] = float(self.total_site_list[r_idx])/float(total_ring_sites)
	    
    def select_random_pair_loc(self,idx,rad_value):
	dist = 0.0
	
	ul, ll = [], []
	
	for ix in range(0,3):
	    ul.append(min(int(self.lattice_size)-1,idx[ix]+self.R_lim+1))
	    ll.append(max(0,idx[ix]-self.R_lim-1))
	
	trial_x = [0,0,0]
	
	dist_val = [0,0,0]
	
	while dist<4.0 or dist>=(self.R_lim+1)**2:
	    for ix in range(0,3):
		trial_x[ix] = rand.randint(idx[ix]-self.R_lim-1,idx[ix]+self.R_lim+1)
		dist_val[ix] = abs(trial_x[ix] - idx[ix])
		
	    dist = self.correlation_distance(idx,trial_x)
	    
	    for ix in range(0,3):
		if dist>=4.0 and dist<(self.R_lim+1)**2:
		    if trial_x[ix] >= self.lattice_size:
			trial_x[ix] = trial_x[ix] - self.lattice_size
		    elif trial_x[ix] < 0:
			trial_x[ix] = self.lattice_size + trial_x[ix]

	return trial_x, dist, dist_val

    def compute_oroganov_order_parameters(self,node_corr,junc_corr):
	self.node_order = 0.0
	self.junction_order = 0.0
	
	for idx in range(0,len(node_corr)):
	    self.node_order += node_corr[idx]**2
	    self.junction_order += junc_corr[idx]**2
	    
	# Normalize values
	self.node_order = self.node_order/len(node_corr)
	self.junction_order = self.junction_order/len(junc_corr)

    def rate_function(self,func_arg):
	func_val = sp.gammainc(1,func_arg) - self.event_val
	
	return func_val
	
    def test_scale_separation(self,rate):
	# Compute the maximum number of fluctuation trials possible
	
	if self.initiation_type=='bimolecular' and self.diffusion_factor==0.0:
	    a_R = self.global_propensities['polymerize'] + self.global_propensities['initiate']
	else:
	    a_R = self.A_total - self.global_propensities['fluctuate']
	
	if a_R>0.0:
	    max_fluctuation_trials = max(1,int(rate/a_R))
	else:
	    max_fluctuation_trials = 1

	return max_fluctuation_trials
    
    def write_summary(self):
	summary_file = open('summary'+str(self.proc_no)+'.csv','w')
	
	print >> summary_file, 'Total nodes',',','Total junctions',',','Total sites',',','Final density'
	print >> summary_file, self.total_nodes,',', self.totalJunctionPoints,',',self.total_sites,',',self.final_density
	
	summary_file.close()
	
    def write_energy_histogram(self):
	histogram_file = open('energy_histogram'+str(self.proc_no)+'.csv','w')
	
	for key in self.total_energy_histogram.keys():
	    print >> histogram_file, str(key)+','+str(self.total_energy_histogram[key])
	    
	histogram_file.close()

    def take_fluctuation_step(self,event_prob):
	#self.check_for_dark_polymerization()
	max_trials = 1
	
	status = self.make_thermal_move()
	
	return status

    def decide_radical_motion(self):
	""" This function is used to compare the total propensity of fluctuation propensity of a radical
	"""
	total_chemical_propensity = sum(self.global_propensities.values()) - self.global_propensities['fluctuate']
	
	if total_chemical_propensity>self.slow_rate:
	    radical_motion = False
	    total_trials = None
	else:
	    radical_motion = True
	    total_trials = self.get_scale_separation(total_chemical_propensity,self.slow_rate)
	    # Also decide the maximum number of radical motion trials possible
	    
	return radical_motion
    
    def get_scale_separation(self,rate_1,rate_2):
	""" This function considers firing rate of two events
	    and computes the number of events of type 2 can happen
	    during the average time period of event 1
	"""
	
	max_trials = max(1,int(rate_1/rate_2))

	return max_trials
	
    def correct_step_size(self,total_trials,initial_guess,event_prob,rate):
	print 'Initial guess: ',initial_guess
	sys.stdout.flush()
	# Use Bisection routine to Newton-Raphson method to estimate time for total_trials>1
	self.event_val = event_prob
	self.event_num = total_trials
	
	initial_rate = sp.gammaincinv(total_trials,event_prob)
	ll = 0.0
	ul = 2.0*initial_rate
	
	#this_rate = bisect(self.rate_function,ll,ul)
	
	real_step_size = initial_guess*sp.gammaincinv(total_trials,event_prob)/sp.gammaincinv(1,event_prob)
	
	if real_step_size<initial_guess:
	    print 'error: ', sp.gammaincinv(1,event_prob), sp.gammaincinv(total_trials,event_prob)
	    sys.exit()
	else:
	    self.fluctuation_rate_factor = initial_guess/real_step_size
	
	print initial_guess, real_step_size
	
	return real_step_size

    def select_a_radical(self):
	max_wt = self.max_radical_fluctuation_propensity/self.total_radical_propensity

	diff = -1.0
	
	while diff<0.0:
	    xi_sample = rand.uniform(0.0,1.0)
	
	    print self.non_zero_radicals
	    sys.stdout.flush()
	    
	    chain_no = rand.choice(self.non_zero_radicals)
	    
	    this_density = self.radical_fluctuation_propensities[chain_no]/self.total_radical_propensity
	    
	    diff = (this_density/max_wt) - xi_sample

	#    print 'radical propensity: ', this_density, max_wt, xi_sample, diff
	#
	#print 'Radical selected'
	#sys.stdout.flush()
	
	if diff>0.0:
	    active_node = self.all_chains[chain_no][-1]
	    return active_node, chain_no
	else:
	    return None, None
	
    def approximate_minimum_fluctuation_propensity(self):
	# No. of radicals
	n_r = len(self.active_chains)
	# No. of regular nodes
	n_nodes = self.total_nodes
	# No. of junction points
	n_jn = self.totalJunctionPoints
	
	# Total transition factor
	total_transition = n_r*self.fluctuator.minimum_escape_probs[0] + (n_nodes - n_jn - n_r)*self.fluctuator.minimum_escape_probs[1] + n_jn*self.fluctuator.minimum_escape_probs[3]
	
	propensity = 6.0*self.fluctuator.mean_propensity*total_transition/float(n_nodes)
	
	return propensity

    def compute_equilibrium_entropy(self):
	"""This routine computes the equilibrium entropy of the system
	    given the current population of nodes, junctions, and radical ends
	"""
	# population fractions
	n_r = float(2.0*len(self.active_chains))/float(self.total_nodes)
	n_n = float(self.total_nodes - 2*len(self.active_chains) - self.totalJunctionPoints)/float(self.total_nodes)
	n_j = float(self.totalJunctionPoints)/float(self.total_nodes)
	
	eq_probs = np.zeros(shape=(11))
	energy_density = np.zeros(shape=(11))
	
	self.eq_entropy = 0.0
	self.old_entropy = self.current_entropy
	self.current_entropy = 0.0
	
	for en_id in range(0,7):
	    eq_probs[en_id] = n_r*self.fluctuator.nr_MB_probs[en_id] + n_j*self.fluctuator.jn_MB_probs[en_id] + n_n*self.fluctuator.nn_MB_probs[en_id]
	    if eq_probs[en_id]>0.0:
		self.eq_entropy += -eq_probs[en_id]*math.log(eq_probs[en_id])
		
	    if self.energy_level_population[en_id]>0:
		energy_density[en_id] = self.energy_level_population[en_id]/self.total_nodes
		self.current_entropy += -energy_density[en_id]*math.log(energy_density[en_id])
		
	print 'Equilibrium entropy: ', self.eq_entropy, sum(self.energy_level_population), self.total_nodes
	print 'Current entropy: ', self.current_entropy, self.old_entropy
	
	return self.eq_entropy, self.current_entropy
    
    def fit_func(self,x,a,b):
	return a*(x-x[0]) + b

    def check_for_equilibrium_leap(self):
	"""This routine compares the current fluctuation propensity with the minimum possible value
	    If the difference is less than 1% then we assume the chain configuration is at equilibrium
	"""
	eq_leap = 'no leap'
	
	if self.DoC<0.01:
	    fluc_limit = self.equilibrium_tolerance
	    window_size = self.minimal_sample_size
	else:
	    fluc_limit = self.equilibrium_tolerance
	    window_size = self.minimal_sample_size

	#window_size = 10#max(int(self.global_propensities['fluctuate']),100)
	#window_size = max(self.total_nodes,2)
	#sample_size = len(self.fluctuation_samples)
	
	if self.fluctuation_sample_size>window_size:# and self.non_fluctuation_propensity>1.0:
	#    for k in xrange(-sample_size,-window_size):
	#	self.fluctuation_sum += -self.fluctuation_samples[k]
		#self.fluctuation_square_sum += -self.fluctuation_samples[-(window_size+1)]**2

	    #rsd = (self.fluctuation_square_sum/self.fluctuation_sum - self.fluctuation_sum/window_size)
	    factor1 = self.fluctuation_sum/float(self.fluctuation_sample_size)#(self.stepurrent_time - self.last_non_fluctuation_step)
	    factor2 = (self.current_time - self.last_non_fluctuation_time)/float(self.fluctuation_sample_size)
	    
	    #factor1 = self.fluctuation_sum/(self.current_time)
	    #factor2 = (self.current_time)/self.step_no
	    
	    self.fluctuation_epr = factor1/factor2
	    
	    #print self.fluctuation_epr
	    
	    #if self.fluctuation_epr<(self.equilibrium_tolerance*self.fluctuator.epr_quanta):
	    if abs(factor1)<fluc_limit:
		eq_leap = 'leap'

	#else:
	#    self.fluctuation_sum += self.fluctuation_samples[-1]
	#    self.fluctuation_square_sum += self.fluctuation_samples[-1]**2
	
	if eq_leap=='no leap':
	    if self.global_propensities['polymerize']==0.0 and self.global_propensities['initiate']>0.0:
		eq_leap = 'leap'

	return eq_leap

    def get_total_diffusion_factor(self):
	"""This function computes the multiplicative factor for observing any fluctuation event in the lattice
	This includes the propensity of seeing a fluctuating network node
	and the propensity of seeing a fluctuating monomer
	"""
	# Propensity factor
	total_propensity = 0.0
	
	for mo_type in xrange(0,self.all_monomers.no_monomers):
	    total_propensity += self.all_monomers.monomer_population[mo_type]*self.fluctuator.lattice_propensity[mo_type]
	
	total_propensity *= 1.0/float(self.total_nodes)
	
	self.global_diffusion_propensity = self.global_propensities['fluctuate']/total_propensity

    def write_chain_coords(self):
	ofile = open('coords.txt','w')
	print >> ofile, 'x','y','z'
	for node_no in range(0,len(self.all_chains[1])):
	    this_node = self.all_chains[1][node_no]
	    print >> ofile, this_node[0], this_node[1], this_node[2]
	    
	ofile.close()

    def get_node_type(self,idx):
	if self.junction_chain[idx[0],idx[1],idx[2]]==0 and self.Nr[idx[0],idx[1],idx[2]]>=1:
	    # max_displacement_energy_index
	    node_type = 0
	elif self.Nr[idx[0],idx[1],idx[2]]>=1:
	    # max_displacement_energy_index
	    node_type = 2
	elif self.junction_chain[idx[0],idx[1],idx[2]]>0:
	    # max_displacement_energy_index
	    node_type = 3
	else:
	    # max_displacement_energy_index
	    node_type = 1
	    
	return node_type

    def check_common_super_lattice_pts(self,idx_1,idx_2):
	set_1 = IndexUtils.get_super_lattice_pts(idx_1)
	set_2 = IndexUtils.get_super_lattice_pts(idx_2)
	
	for pts in set_2:
	    if pts not in set_1:
		print pts, self.super_lattice[pts[0],pts[1],pts[2]]
		
    def all_stdouts(self,step_type,step_size):
	print 'Event propensities:'
	
	for event_type in self.global_propensities.keys():
	    if event_type=='fluctuate':
		print event_type, ': ', (self.all_monomers.mean_d_factors*self.fluctuator.mean_lattice_propensity)*self.global_propensities[event_type], self.fluctuator.mean_lattice_propensity
	    else:
		print event_type, ': ', self.global_propensities[event_type]
	
	print 'Step type : ', step_type, (time.time()-self.start_time), self.max_fluc_trials
	#print 'Step size : ', step_size, self.A_total
	self.start_time = time.time()
	
	self.max_fluc_trials = 0
	
	print 'Step No : ', self.step_no, 'Time: ', self.current_time, 'DoC', self.DoC
	print 'Movable nodes : ', self.total_highly_movable_nodes, self.total_movable_nodes, self.total_nodes, self.max_fluc_wt
	print 'Thermal propensities : ', self.highly_movable_propensities, self.global_propensities['fluctuate']
	print '\nActive chains : ', len(self.active_chains), len(self.all_chains)
	
	#print 'Molecule diffusivity: ', self.all_monomers.mean_eps_self, self.all_monomers.mean_d_factors, math.exp(-0.5*self.all_monomers.eps_self[0]/(0.0259*300.0/self.temp))
	#print 'Oligomer diffusivity: ', self.oligomer_eps_self, self.diffusion.compute_effective_diffusivity(self.oligomer_eps_self)
	print 'Node diffusivity: ', self.mean_diffusion_factor
	print 'Lattice propensity: ', self.lattice_propensity
	print 'EPR: ', self.fluctuation_epr, self.fluctuation_sample_size
	print 'Total bonds: ', self.bond_populations['high']
	
	#print 'Total propensities: ', self.A_total
	#print 'Surface sites: ', self.surface_site_count
	#print 'Neighbor sites: ', self.neighbor_site_count
	sys.stdout.flush()

    def write_vtk_files(self):
	self.vtk_counter += 1
	filename = vtk_writer.case_and_step_file(str(self.vtk_counter))

	ofile = open(filename,'a')
	box_size = int(self.lattice_size)

	for x in range(0,box_size):
	    for y in range(0,box_size):
		for z in range(0,box_size):
		    neighset = IndexUtils.nearest_neighbors([x,y,z])
		    thermal_set = IndexUtils.sites_for_thermal_moves([x,y,z])
		    neighset.append([x,y,z])
		    thermal_set.append([x,y,z])
		    value = self.check_node_juncs(thermal_set,[x,y,z])
		    
		    #super_lat_value = self.check_super_lattice([x,y,z])
		    #pixel_wt = 8 -super_lat_value
		    
		    print >> ofile, str(value)
		##    if self.super_lattice[x,y,z]>0:
		##	print >> ofile, '1'
		#    if self.Nr[x,y,z]>0:
		#	print >> ofile, '3'
		#    elif self.junction_chain[x,y,z]>0:
		#	print >> ofile, '2'
		#    elif self.location_chain[x,y,z]>0:
		#	print >> ofile, '1'
		##    elif self.check_super_lattice([x,y,z])==0:
		##	print >> ofile, '2'
		##    elif self.check_super_lattice([x,y,z])<8:
		##	print >> ofile, '1'
		#    else:
		#	print >> ofile, '0'

	ofile.close()
	
    def check_node_juncs(self,indset,loc_id):
	value = 0
	
	for idx in indset:
	    x, y, z = idx[0], idx[1], idx[2]
	    
	#    if self.Nr[x,y,z]>0 and value<3:
	#	value = 2
	    if self.junction_chain[x,y,z]>0 and value<3:
		value = 1
	    elif self.location_chain[x,y,z]>0 and value<2:
		value = 1
	
	x, y, z = loc_id[0], loc_id[1], loc_id[2]
	
	if self.location_chain[x,y,z]>0:
	    value = 2
	if self.junction_chain[x,y,z]>0:
	    value = 3

	return value
	
    def create_ordered_lattice(self):
	box_size = int(self.lattice_size)
	self.total_nodes = 0
	self.all_chains[1] = []
	
	for x in range(0,box_size):
	    for y in range(0,box_size):
		for z in range(0,box_size):
		    if x%3==0 and y%3==0 and z%3==0:
			self.total_nodes += 1
			self.junction_chain[x,y,z] = 1
			self.all_chains[1].append([x,y,z])
			self.update_super_lattice([x,y,z],'add')

    def compute_ordered_correlation(self):
	id_max = int(self.lattice_size)
	
	# Monte Carlo sampling size
	#mc_sample_size = self.total_sites*int(math.log10(self.total_sites))
	mc_sample_size = self.total_sites*int(math.log(self.total_sites))
	#mc_sample_size = 100000
	print mc_sample_size
	sys.stdout.flush()
	
	sample_count = 0
	n_site_count, j_site_count = 0,0
	node_pair_count = 0
	junction_pair_count = 0
	
	N_ls_pop_j = np.zeros(shape=(len(self.rad_list)))
	
	N_js_pop = np.zeros(shape=(len(self.rad_list)))
	
	# Pair correlation functions
	N_ls_pop_x_n = np.zeros(shape=(len(self.rad_list)+1))
	N_ls_pop_y_n = np.zeros(shape=(len(self.rad_list)+1))
	N_ls_pop_z_n = np.zeros(shape=(len(self.rad_list)+1))
	
	N_js_x_n = np.zeros(shape=(len(self.rad_list)+1))
	N_js_y_n = np.zeros(shape=(len(self.rad_list)+1))
	N_js_z_n = np.zeros(shape=(len(self.rad_list)+1))
	
	# Pair correlation functions
	g_r_n = np.zeros(shape=(len(self.rad_list)))
	g_r_j = np.zeros(shape=(len(self.rad_list)))
	
	while junction_pair_count<mc_sample_size:
	    # Pick the first random node from the network
	    chain_no = rand.randint(1,len(self.all_chains))
	    node_no = rand.randint(0,len(self.all_chains[chain_no])-1)
	    j_id = self.all_chains[chain_no][node_no]
	    sample_count += 1
	    
	    j_id2, dist, c_idx2 = self.select_random_pair_loc(j_id)

	    if int(dist)>1.0 and int(dist)<(self.R_lim+1):
	        r_idx = int(dist)-2
		N_ls_pop_j[r_idx] += 1
		j_site_count += 1
		
		if c_idx2[0]>0:
		    N_ls_pop_x_n[c_idx2[0]-1] += 1
		if c_idx2[1]>0:
		    N_ls_pop_y_n[c_idx2[1]-1] += 1
		if c_idx2[2]>0:
		    N_ls_pop_z_n[c_idx2[2]-1] += 1 
		 
		if self.junction_chain[j_id[0],j_id[1],j_id[2]]>0 and self.junction_chain[j_id2[0],j_id2[1],j_id2[2]]>0:
		    N_js_pop[r_idx] += 1

		    if c_idx2[0]>0:
			N_js_x_n[c_idx2[0]-1] += 1
		    if c_idx2[1]>0:
			N_js_y_n[c_idx2[1]-1] += 1
		    if c_idx2[2]>0:
			N_js_z_n[c_idx2[2]-1] += 1

		    junction_pair_count += 1
		    
	# Junction density
	rho_jn_true = float(self.total_nodes)/float(self.total_sites)
	rho_jn = float(junction_pair_count)/float(j_site_count)
	   
	pair_file = open('pair_correlations'+str(self.proc_no)+'.csv','w')
	period_file = open('periodicity'+str(self.proc_no)+'.csv','w')
	
	print >> pair_file, 'Radial distance,Junction correlation'
	
	print >> period_file, 'Lattice distance,Node X,Node Y, Node Z, Junction X, Junction Y, Junction Z'
	
	ntotal_j = sum(N_ls_pop_j)
	
	jn_correlations = np.zeros(shape=(len(self.rad_list)))
	
	for r_idx in range(0,len(self.rad_list)):
	    jn_correlations[r_idx] = float(N_js_pop[r_idx])/(float(N_ls_pop_j[r_idx])*rho_jn)
	    
	    x_j_p = float(N_js_x_n[r_idx])/(float(N_ls_pop_x_n[r_idx])*rho_jn)
	    y_j_p = float(N_js_y_n[r_idx])/(float(N_ls_pop_y_n[r_idx])*rho_jn)
	    z_j_p = float(N_js_z_n[r_idx])/(float(N_ls_pop_z_n[r_idx])*rho_jn)
	    
	    print >> pair_file, self.rad_list[r_idx], ',', jn_correlations[r_idx]
	    
	    print >> period_file, self.rad_list[r_idx], ',', x_j_p, ',', x_j_p, ',', x_j_p

	pair_file.close()
	period_file.close()

    def check_for_dark_polymerization(self):
	total_chem_prop = self.global_propensities['initiate']
	
	try:
	    total_chem_prop += self.global_propensities['excite']
	except KeyError:
	    pass
	
	if total_chem_prop==0.0:
	    return True
	else:
	    return False
	
    def compute_free_volume(self):
	self.free_vol = 0
	
	for ix in range(0,int(self.lattice_size)):
	    for iy in range(0,int(self.lattice_size)):
		for iz in range(0,int(self.lattice_size)):
		    if self.location_chain[ix,iy,iz]==0 and self.check_super_lattice([ix,iy,iz])=='empty':
			all_sites = IndexUtils.sites_for_thermal_moves([ix,iy,iz])
			count = 0
			
			for site in all_sites:
			    print self.check_super_lattice_thermal(site), ',', count
			    if self.check_super_lattice_thermal(site)=='empty':
				count += 1
			
			if count==0:
			    self.free_vol += 1
			    
	print self.free_vol

    def get_directional_stiffness(self,moves,chain_no,node_no,second_chain,node_no_2):
	K_values = [0,0,0,0,0,0]
	total_energy_change = [0,0,0,0,0,0]
	
	for move_no in xrange(0,len(moves)):
	    energy_change = 0.0
	    
	    if moves[move_no]!=0:
		self.move = IndexUtils.thermalMoves[move_no]
	
		check, U1_L1, U1_L2, U1_ac, U_a1, U_a2, b1, b2 = self.check_for_motion(chain_no,node_no)

		if check=='good move':
		    energy_change = U1_L1+U1_L2+U1_ac+U_a1+U_a2

		    if second_chain!=None:
			second_check, U2_L1, U2_L2, U2_ac, U_sa1, U_sa2, sb1, sb2 = self.check_for_motion(second_chain,node_no_2)
			
			if second_check=='good move':
			    energy_change += (U2_L1+U2_L2+U2_ac+U_sa1+U_sa2)
		else:
		    energy_change = 0.0

	    for x_idx in xrange(0,3):
		if moves[move_no]==1:
		    total_energy_change[x_idx] = energy_change
		elif moves[move_no]==-1:
		    total_energy_change[x_idx+3] = energy_change
	
	for idx in xrange(0,6):
	    K_values[idx] = self.get_stiffness(total_energy_change[idx])

	return K_values

    def get_stiffness(self,energy_change):
	K = abs(energy_change/(self.lattice_h_final**2))
	
	return K

    def get_omega(self,local_stiffness,sample_k,lattice_h):
	omega = []
	
	for k_idx in xrange(0,3):
	    factor_1 = math.sqrt(local_stiffness[k_idx]*self.k_quanta/self.mol_wt_Kg)
	    omega_value = 2.0*factor_1*abs(math.sin(0.5*self.lattice_h_final*sample_k[k_idx]))
	    omega.append(omega_value)
	
	return omega

    def get_second_chain_and_box_no(self,chain_no,box_no):
	idxs = self.all_chains[chain_no][box_no]
	
	if self.location_chain[idxs[0],idxs[1],idxs[2]] == chain_no:
	    second_chain = self.junction_chain[idxs[0],idxs[1],idxs[2]]
	else:
	    second_chain = self.location_chain[idxs[0],idxs[1],idxs[2]]
	    
	if second_chain != chain_no:
	    second_box = self.all_chains[second_chain].index(idxs)    
	else:
	    for k, try_idx in enumerate(self.all_chains[second_chain]):
		if try_idx==idxs and k!=box_no:
		    second_box = k
		    
	return second_chain, second_box

    def compute_phonon_distribution(self):
	"""This function calculates the Vibrational Density of States
	for the lattice polymer network.
	"""
	self.lattice_h_final = self.lattice_h#*((self.DoC*self.shrinkage_factor)**(1.0/3.0))
	
	# Frequency scaling factor
	self.freq_scale = math.sqrt((1.602E-19)/(((1E-10)**2)*(1.66054E-27)))
	
	# Stiffness quanta (in N/m)
	self.k_quanta = 1.0#2.0*(self.eps_self/(self.lattice_h_final**2))
	
	# Moleculart weight
	self.mol_wt_Kg = self.molWt
	
	# Monte Carlo sampling size
	global_mass = 0
	self.occupied_holes = (self.total_effective_units-self.total_nodes)/float(self.total_sites)
	
	for chain_no in self.all_chains.keys():
	    chain_length = len(self.all_chains[chain_no])
	    for node_no in xrange(0,len(self.all_chains[chain_no])):
		idx = self.all_chains[chain_no][node_no]
		if self.junction_chain[idx[0],idx[1],idx[2]]<=chain_no:
		    this_neighbors = IndexUtils.get_mass_neighbors(idx)
		    
		    local_mass = 0
		    
		    for loc_idx in this_neighbors:
			local_mass += min(1,self.location_chain[loc_idx[0],loc_idx[1],loc_idx[2]])

		    global_mass += local_mass

		    total_stiffness = [0.0,0.0,0.0,0.0,0.0,0.0]

		    if self.junction_chain[idx[0],idx[1],idx[2]]>0:
			second_chain, node_no_2 = self.get_second_chain_and_box_no(chain_no,node_no)
		    else:
			second_chain, node_no_2 = None, None
	
		    lidx = self.get_lidx_from_cubidx(idx)
		    
		    neighbor_stiffness = self.compute_interaction_stiffness(this_neighbors)
		    
		    if chain_length>1:
			bond_stiffness = self.get_bond_stiffness(chain_no,node_no,second_chain,node_no_2)
	
		    for k in xrange(0,6):
			total_stiffness[k] = neighbor_stiffness[k] + bond_stiffness[k]
			
		    self.disorder_calculator.add_stiffness_and_mass_to_bin(total_stiffness,local_mass)

	self.disorder_calculator.compute_histograms(self.total_effective_units,global_mass)
	self.disorder_calculator.write_histograms(self.proc_no)
	#self.disorder_calculator.compute_mean_stiffness()
	#self.disorder_calculator.initialize_structs()
	#self.disorder_calculator.compute_Gamma_omega()
	#self.disorder_calculator.compute_DoS_MFP()
	#self.disorder_calculator.write_DoS_MFP(self.proc_no)
	#self.disorder_calculator.write_gamma_file(self.proc_no)
	#self.disorder_calculator.compute_thermal_properties(self.proc_no)
	
    def compute_interaction_stiffness(self,neighbors):
	neighbor_stiffness = [0,0,0,0,0,0]
	
	for neighbor_no in IndexUtils.mass_neigbors_set:
	    loc_stiff = None
	    idx = neighbors[neighbor_no]
	    
	    if self.location_chain[idx[0],idx[1],idx[2]]>0:
		loc_stiff = IndexUtils.stiff_vecs[neighbor_no]
	    elif self.check_super_lattice(idx)==8:
		xi = rand.uniform(0.0,1.0)
		
		if xi<self.occupied_holes:
		    loc_stiff = IndexUtils.stiff_vecs[neighbor_no]

	    if loc_stiff!=None:
		for k in [0,1,2]:
		    if loc_stiff[k]>=0:
			neighbor_stiffness[k] += loc_stiff[k]
		    else:
			neighbor_stiffness[k+3] += abs(loc_stiff[k])

	return list(neighbor_stiffness)

    def get_bond_stiffness(self,chain_no,node_no,second_chain,node_no_2):
	bond_stiffness = [0,0,0,0,0,0]
	#factor = 10.0*self.eps_self/(self.lattice_h_final**2)
	
	center = self.all_chains[chain_no][node_no]
	
	# Find the location of the surrounding two nodes
	try:
	    prev_1 = self.all_chains[chain_no][node_no-1]
	    if center==prev_1:
		print center, prev_1, len(self.all_chains[chain_no])
		sys.stdout.flush()
	    else:
		unit_vec, norm = IndexUtils.get_unit_vector(center,prev_1)
		
		for x_idx in xrange(0,3):
		    if unit_vec[x_idx]>0:
			bond_stiffness[x_idx] += unit_vec[x_idx]#/(norm**2)
		    elif unit_vec[x_idx]<0:
			bond_stiffness[x_idx+3] += abs(unit_vec[x_idx])#/(norm**2)
	except IndexError:
	    pass

	try:
	    next_1 = self.all_chains[chain_no][node_no+1]
	    if center==next_1:
		print center, next_1, len(self.all_chains[chain_no])
		sys.stdout.flush()
	    else:
		unit_vec, norm = IndexUtils.get_unit_vector(center,next_1)
		
		for x_idx in xrange(0,3):
		    if unit_vec[x_idx]>0:
			bond_stiffness[x_idx] += unit_vec[x_idx]#/(norm**1)
		    elif unit_vec[x_idx]<0:
			bond_stiffness[x_idx+3] += abs(unit_vec[x_idx])#/(norm**1)
	except IndexError:
	    pass
	
	if second_chain!=None:
	    # Find the location of the surrounding two nodes
	    try:
		prev_1 = self.all_chains[second_chain][node_no_2-1]
		if center==prev_1:
		    print center, prev_1, len(self.all_chains[second_chain])
		    sys.stdout.flush()
		else:
		    unit_vec, norm = IndexUtils.get_unit_vector(center,prev_1)
		    
		    for x_idx in xrange(0,3):
			if unit_vec[x_idx]>0:
			    bond_stiffness[x_idx] += unit_vec[x_idx]#/(norm**1)
			elif unit_vec[x_idx]<0:
			    bond_stiffness[x_idx+3] += abs(unit_vec[x_idx])#/(norm**1)
	    except IndexError:
		pass

	    try:
		next_1 = self.all_chains[second_chain][node_no_2+1]
		if center==next_1:
		    print center, next_1, len(self.all_chains[second_chain])
		    sys.stdout.flush()
		else:
		    unit_vec, norm = IndexUtils.get_unit_vector(center,next_1)
		    
		    for x_idx in xrange(0,3):
			if unit_vec[x_idx]>0:
			    bond_stiffness[x_idx] += unit_vec[x_idx]#/(norm**2)
			elif unit_vec[x_idx]<0:
			    bond_stiffness[x_idx+3] += abs(unit_vec[x_idx])#/(norm**2)
	    except IndexError:
		pass

	return list(bond_stiffness)

    def get_local_density(self,neighbor_sites):
	occupied_sites = 1
	max_mass = 27.0
	
	for site_idx in neighbor_sites:
	    vertex_count = self.check_super_lattice(site_idx)
	    
	    if vertex_count==8 or self.location_chain[site_idx[0],site_idx[1],site_idx[2]]>0:
		occupied_sites += 1

	local_density = float(occupied_sites)/max_mass
	
	return local_density

#    def compute_phonon_distributions(self):
#	for chain_no in self.all_chains.keys():
#	    for node_num in xrange(0,len(self.all_chains[chain_no])-1):
#		total_stiffness = [0.0,0.0,0.0,0.0,0.0,0.0]
#		dist = IndexUtils.get_distance(self.all_chains[chain_no][node_num],self.all_chains[chain_no][node_num+1])
#		idx = self.all_chains[chain_no][node_num]
#		if self.junction_chain[idx[0],idx[1],idx[2]]>0:
#		    second_chain, node_no_2 = self.get_second_chain_and_box_no(chain_no,node_num)
#		else:
#		    second_chain, node_no_2 = None, None
#		    
#		neighbor_stiffness = self.compute_interaction_stiffness(idx)
#
#		neighbor_sites = IndexUtils.get_mass_location_sites(idx)
#		
#		bond_stiffness = self.get_bond_stiffness(chain_no,node_num,second_chain,node_no_2)
#		
#		local_density = self.get_local_density(neighbor_sites)
#		
#		for k in xrange(0,6):
#		    total_stiffness[k] = neighbor_stiffness[k] + bond_stiffness[k]
#
#		self.disorder_calculator.add_stiffness_and_density_to_bin(total_stiffness,local_density)
#		    
#		if dist<2.0:
#		    print 'Bond error'
#		    print 'Chain no: ', chain_no
#		    print 'Chain lenght: ', len(self.all_chains[chain_no])
#		    print 'Wrong node: ', node_num
#		    sys.stdout.flush()
#		    sys.exit()
#	
#	#if self.proc_no==1:
#	#self.disorder_calculator.setup_scaling_factors(self.eps_self,self.lattice_h,self.molWt)
#
#	self.disorder_calculator.compute_histograms()
#	self.disorder_calculator.write_histograms(self.proc_no)
#	self.disorder_calculator.compute_mean_stiffness()
#	self.disorder_calculator.initialize_structs()
#	self.disorder_calculator.compute_Gamma_omega()
#	self.disorder_calculator.compute_DoS_MFP()
#	self.disorder_calculator.write_DoS_MFP(self.proc_no)
#	self.disorder_calculator.write_gamma_file(self.proc_no)
#	#self.disorder_calculator.compute_thermal_properties(self.proc_no)
#		    
#	print 'Checking for consecutive nodes Successful'
#	sys.stdout.flush()
	
    def check_total_effective_units(self):
	occupied_sites = 0
	
	for x_idx in xrange(0,int(self.lattice_size)):
	    for y_idx in xrange(0,int(self.lattice_size)):
		for z_idx in xrange(0,int(self.lattice_size)):
		    vertex_count = self.check_super_lattice([x_idx,y_idx,z_idx])

		    if vertex_count==8 or self.location_chain[x_idx,y_idx,z_idx]>0:
			occupied_sites += 1

	print 'Total effective units in lattice: ', occupied_sites, self.prop_count
	print 'Out of a total of: ', (self.lattice_size_half**3)
	
	self.final_density = occupied_sites*self.molWt*(1.66E-24)/(self.sysVol*(1E-24))
	
    def forced_extension(self,selected_chain):
	""" Uses the polymer reaction propensity data to decide which chain is going to react next
	and the kind of reaction that is expected.
	"""
        # Random number to select the type of reaction
        xi3 = rand.uniform(1e-10,1.0)
        # Random number to select the type of reaction
        xi4 = rand.uniform(1e-10,1.0)
        temp_sum = 0.0
    
        # Get indices of chain ends
        Alist = self.totalProbSet[selected_chain]
        funcValues = self.neighbor_groups[selected_chain]
	site_set = self.non_zero_sites[selected_chain]
	
        last_idx_set = self.all_chains[selected_chain][-1]
	last_monomer = self.location_monomers[last_idx_set[0],last_idx_set[1],last_idx_set[2]]
	# Get the second last index set
	try:
	    prev_idx = self.all_chains[selected_chain][-2]
	#    if self.blocked_matrix[prev_idx[0],prev_idx[1],prev_idx[2]]==0:
	#	print 'Error: ', prev_idx, self.all_chains[selected_chain]
	#	sys.stdout.flush()
	#	sys.exit()
	except IndexError:
	    prev_idx = None
	    
        # Get the list of next step voxels
        IndSet = IndexUtils.sitesfornextmove(last_idx_set)
	
	if last_idx_set in IndSet:
	    print 'Index Error: ', selected_chain
	    sys.stdout.flush()
	    sys.exit()

        # Now find which neighboring box is going to react
        sumAlist = sum(Alist)
          
        temp_sum = 0.0

	propWt, cycleWt, termWt = 0.0, 1.0, 1.0

	while termWt>0.0 and propWt==0.0:
	    for bno in range(0,len(site_set)):
		thisProb = Alist[bno]/sumAlist

		# If found select chain no and exit
		if (temp_sum < xi3) and xi3 <=(temp_sum+thisProb):
		    selected_box = self.non_zero_sites[selected_chain][bno]
		    funcValue = self.neighbor_groups[selected_chain][bno]
		    monomer = self.neighbor_monomers[selected_chain][bno]
		    break
		else:
		    temp_sum += thisProb

	    # Now decide whether reaction is propagation or termination
	    idx = IndSet[selected_box]

	    if funcValue==self.all_monomers.functionality[int(monomer)]:
		propWt = self.all_monomers.k_p[last_monomer,monomer]*funcValue

	    termWt = self.all_monomers.k_t[last_monomer,monomer]*self.Nr[idx[0],idx[1],idx[2]]

        # Decide if there is any reaction, and whether the reaction is 
        # propagation and termination
        totalWt = propWt + cycleWt + termWt

        if totalWt == 0.0:
            print 'Something wrong', Alist, bno, Alist[bno], self.problist[selected_chain]
	    sys.stdout.flush()
            sys.exit()
            react = None    

        react = 'prop'
	# One functional group taken at idx
	self.location_functionalities[idx[0],idx[1],idx[2]] = funcValue - 1
	self.transientFuncs += -1
	self.transientMonomers += -1
	self.total_nodes += 1
		
	self.Nr[last_idx_set[0],last_idx_set[1],last_idx_set[2]] += -1
	self.Nr[idx[0],idx[1],idx[2]] += 1
	
	self.location_chain[idx[0],idx[1],idx[2]] = selected_chain
	self.location_monomers[idx[0],idx[1],idx[2]] = monomer
	
	if IndexUtils.get_distance(idx,last_idx_set) < 2.0:
	    print 'prop error: ', idx, last_idx_set
	    sys.stdout.flush()
	    sys.exit()

	self.neighbor_set = IndexUtils.nearest_neighbors(idx)

	initial_surfaces = self.local_surface_count(self.neighbor_set)
	site_values = self.check_neighbors(self.neighbor_set)
	self.update_super_lattice(idx,'add')
	final_surfaces = self.local_surface_count(self.neighbor_set)
	self.surface_site_count += final_surfaces - initial_surfaces
	
	#self.update_neighbors(idx,self.neighbor_set,site_values,'add')
	self.energy_level_population[0] += 1
	
	self.last_update_box = idx

        return idx, monomer
	
    def update_chain_numbers(self,chain_no):
	for chain_key in self.all_chains.keys():
	    if chain_key>chain_no:
		self.all_chains[chain_key-1] = self.all_chains[chain_key]
		del self.all_chains[chain_key]
	
	for chain_key in self.terminated_pairs.keys():
	    if chain_key>chain_no:
		self.terminated_pairs[chain_key-1] = self.terminated_pairs[chain_key]
		del self.terminated_pairs[chain_key]

	for chain_key in self.terminated_pairs.keys():
	    if self.terminated_pairs[chain_key]>chain_no:
		self.terminated_pairs[chain_key] += -1
		
    def check_radicals(self):
	for chain_no in self.active_chains:
	    idx = self.all_chains[chain_no][-1]

	    #print 'Chain no: ', chain_no, ' is ', self.Nr[idx[0],idx[1],idx[2]], idx
	    
	    if self.Nr[idx[0],idx[1],idx[2]]==0:
		#self.Nr[idx[0],idx[1],idx[2]] = 1.0
		sys.exit()
		
    def compute_crosslinks_prob(self,tag):
	total_samples = 0
	total_prob_values = 0
	
	for idx_set in self.junction_list:
	    #next_sites = IndexUtils.sitesfornextmove(idx_set)
	    next_sites, next_idx = IndexUtils.get_blocking_monomers(idx_set,self.location_chain)
	    local_junctions = 0
	    local_sites = 0
	    total_samples += 1
	    
	    for site in next_sites:
		if self.location_chain[site[0],site[1],site[2]]>0:
		    local_sites += 1
		if self.junction_chain[site[0],site[1],site[2]]>0:
		    local_junctions += 1
	    
	    if local_sites!=0:
		total_prob_values += float(local_junctions)/float(local_sites)
	    else:
		total_prob_values += 0

	self.conditional_prob = total_prob_values/float(total_samples)
	self.regular_density = float(total_samples)/float(self.TotalMonomers)
	
	print 'Conditional: ', self.regular_density, ' and ', self.conditional_prob
	
	summary_file = open('crosslink_summary_'+tag+'_'+str(self.proc_no)+'.csv','w')
	
	print >> summary_file, 'Crosslinks density',',','Conditional Probability'
	print >> summary_file, self.regular_density,',', self.conditional_prob
	
	summary_file.close()
	
    def compute_jump_prob(self,idxs):
	chain_no = self.location_chain[idxs[0],idxs[1],idxs[2]]
	box_no = int(self.location_boxno[idxs[0],idxs[1],idxs[2]])
	#box_no2 = self.all_chains[chain_no].index(idxs)
	#if box_no!=box_no2:
	#    print 'Computing jump prob: ', chain_no, box_no, box_no2
	#    print idxs, self.all_chains[chain_no]
	#    sys.exit()
	
	lidx = self.get_lidx_from_cubidx(idxs)
	trials = 0
	move_dict = {}
	move_count = 0
	move_prob_dict = {}
	max_move_prob = 0.0
	bond_update = {}
	
	# Initialize bond update dict
	self.bond_pop_update[4.0], self.bond_pop_update[5.0], self.bond_pop_update[6.0] = 0, 0, 0
	self.bond_pop_update[9.0], self.bond_pop_update[10.0] = 0, 0
	
	prob_total = 0.0
	
	#while self.movability[idxs[0],idxs[1],idxs[2]]>0 and check=='bad move' and status=='failure':
	# Default energy allowed move
	energy_decision = False
	
	for move_no in xrange(0,6):
	    #check = 'bad move'
	    move_dict[move_no] = 0.0
	    move_prob_dict[move_no] = 0.0
	    bond_update[move_no] = 0
	    this_prob = 0.0

	    trial_move = IndexUtils.thermalMoves[move_no]
	    
	    check, U1_L1, U1_L2, U1_ac, U_a1, U_a2, b1, b2 = self.check_for_motion(chain_no,box_no,trial_move)
	    
	    if check=='good move':
		solvent_prob = self.check_motion_interaction(chain_no,box_no,trial_move)

	    total_length_energy_change = U1_L1+U1_L2
	    total_angle_energy_change = U1_ac+U_a1+U_a2
	    
	    if check == 'good move':
		if self.junction_chain[idxs[0],idxs[1],idxs[2]] > 0:
		    check = 'bad move'
		
		    if self.location_chain[idxs[0],idxs[1],idxs[2]] == chain_no:
			second_chain = self.junction_chain[idxs[0],idxs[1],idxs[2]]
			second_box = int(self.junction_boxno[idxs[0],idxs[1],idxs[2]])
		    else:
			second_chain = self.location_chain[idxs[0],idxs[1],idxs[2]]
			second_box = int(self.location_boxno[idxs[0],idxs[1],idxs[2]])

		    check, U2_L1, U2_L2, U2_ac, U_sa1, U_sa2, sb1, sb2 = self.check_for_motion(second_chain,second_box,trial_move)

		    if check == 'good move':
			total_length_energy_change += U2_L1+U2_L2
			total_angle_energy_change += U2_ac+U_sa1+U_sa2
			
			length_index = int(total_length_energy_change/self.fluctuator.bond_energy_quanta)
			angle_index = 0
			#angle_index = int(total_angle_energy_change/self.fluctuator.bond_angle_energy_quanta)
			this_key = (length_index,angle_index)
			
			this_prob = self.fluctuator.transition_odds[this_key]*solvent_prob
		else:
		    length_index = int(total_length_energy_change/self.fluctuator.bond_energy_quanta)
		    angle_index = 0
		    #angle_index = int(total_angle_energy_change/self.fluctuator.bond_angle_energy_quanta)
		    this_key = (length_index,angle_index)
		    
		    this_prob = self.fluctuator.transition_odds[this_key]*solvent_prob
		
	    if this_prob>0.0:
		move_dict[move_no] = 1.0
		move_prob_dict[move_no] = this_prob
		bond_update[move_no] = int(total_length_energy_change/self.fluctuator.energy_eps)
		
		max_move_prob = max(max_move_prob,this_prob)
		move_count += 1
		prob_total += this_prob
	
	return prob_total, move_dict, move_count, move_prob_dict, max_move_prob, bond_update

    def check_motion_interaction(self,chain_no,box_no,trial_move):
	self.check_for_motion(chain_no,box_no,trial_move)
	idx = self.all_chains[chain_no][box_no]
	this_node = self.location_monomers[idx[0],idx[1],idx[2]]
	
	initial_energy, new_energy = 0.0, 0.0
	
	neighbors = IndexUtils.get_short_mass_neighbors(idx)
	
	for id2 in neighbors:
	    loc_monomer = self.location_monomers[id2[0],id2[1],id2[2]]
	    if loc_monomer>-1:
		initial_energy += self.all_monomers.intermat[this_node,loc_monomer]
	    else:
		if self.check_super_lattice(id2)=='empty':
		    sample_monomer = self.all_monomers.sample_a_monomer()
		    initial_energy += self.all_monomers.intermat[this_node,sample_monomer]
		    
	new_pos = [0,0,0]
	
	for i in range(0,3):
            new_pos[i] = idx[i] + trial_move[i]

            if new_pos[i] >= self.lattice_size:
                new_pos[i] = new_pos[i] - self.lattice_size
            elif new_pos[i] < 0:
                new_pos[i] = self.lattice_size + new_pos[i]

	neighbors = IndexUtils.get_short_mass_neighbors(new_pos)
	
	for id2 in neighbors:
	    loc_monomer = self.location_monomers[id2[0],id2[1],id2[2]]
	    if loc_monomer>-1:
		new_energy += self.all_monomers.intermat[this_node,loc_monomer]
	    else:
		if self.check_super_lattice(id2)=='empty':
		    xi = rand.uniform(0.0,1.0)
		    
		    if xi<self.all_monomers.occupancy:
			sample_monomer = self.all_monomers.sample_a_monomer()
			new_energy += self.all_monomers.intermat[this_node,sample_monomer]
		    
	delta_energy = (new_energy-initial_energy)
		    
	this_key = int(delta_energy/self.fluctuator.beta)
	
	try:
	    this_prob = self.fluctuator.solvent_prob_keys[this_key]
	except:
	    self.fluctuator.solvent_prob_keys[this_key] = min(1.0,math.exp(-delta_energy/self.fluctuator.beta))

	return self.fluctuator.solvent_prob_keys[this_key]

    def get_prob(self,energy_change):
	transition_odds = math.exp(-energy_change/self.fluctuator.beta)
	
	this_prob = min(1.0,transition_odds)

	return this_prob

    def get_regular_surface_count(self,sites):
	surface_site_update = 0
	
	for idx in sites:
	    if self.neighbor_sites[idx[0],idx[1],idx[2]]>0:
		surface_site_update += 1
		
	return surface_site_update

    def get_junction_surface_count(self,sites):
	surface_site_update = 0
	
	for idx in sites:
	    if self.junction_neighbor_sites[idx[0],idx[1],idx[2]]>0:
		surface_site_update += 1
		
	return surface_site_update
	
    def check_junction_super_lattice(self,idx):
	IdxSet = IndexUtils.get_super_lattice_pts(idx)
	count = 0
	
	for idx in IdxSet:
	    if self.junction_super_lattice[idx[0],idx[1],idx[2]] == 1:
		count += 1

	return count

    def update_junction_super_lattice(self,idx,update_type):
	IdxSet = IndexUtils.get_super_lattice_pts(idx)

	for idx in IdxSet:
	    if update_type == 'add':
		self.junction_super_lattice[idx[0],idx[1],idx[2]] = 1

            if update_type == 'remove':
		self.junction_super_lattice[idx[0],idx[1],idx[2]] = 0
	
    def update_neighbor_sites(self,sites,update_type):
	for idx in sites:
	    if update_type=='add':
		self.neighbor_sites[idx[0],idx[1],idx[2]] += 1
	    elif update_type=='remove':
	    	self.neighbor_sites[idx[0],idx[1],idx[2]] += -1

    def update_junction_neighbors(self,primary_sites,update_type):
	for idx in primary_sites:
	    if update_type=='add':
		self.junction_neighbor_sites[idx[0],idx[1],idx[2]] += 1

	    if update_type=='remove':
	    	self.junction_neighbor_sites[idx[0],idx[1],idx[2]] += -1

    def get_total_sites(self,list_1,list_2):
	temp = list(list_1)
	
	for idx in list_2:
	    if idx not in list_1:
		temp.append(idx)
	
	return temp

    def update_polymerization(self,idx):
	if self.neighbor_sites[idx[0],idx[1],idx[2]]>0:
	    self.neighbor_sites[idx[0],idx[1],idx[2]]=0
	    self.surface_site_count += -1
	if self.junction_neighbor_sites[idx[0],idx[1],idx[2]]>0:
	    self.junction_neighbor_sites[idx[0],idx[1],idx[2]]=0
	    self.junction_surface_count += -1
	    
    def check_all_bond_lengths(self):
	for chain in self.all_chains:
	    length = len(self.all_chains[chain])
	    
	    for idx in self.all_chains[chain]:
		if self.location_chain[idx[0],idx[1],idx[2]]==0:
		    print 'Wrong network value'
		    sys.exit()

	#    if length>1:
	#	for node_no in xrange(0,length-1):
	#	    if IndexUtils.get_distance2(self.all_chains[chain][node_no],self.all_chains[chain][node_no+1]) not in IndexUtils.length_set:
	#		print 'Wrong segment length'
	#		sys.exit()
	
    def compute_radius_of_gyration(self):
	mean_c = [0,0,0]
	
	radius = 0.0
	
	for idx in self.all_chains[1]:
	    for j in [0,1,2]:
		radius += idx[j]**2
		
		mean_c[j] += idx[j]
		
	for j in [0,1,2]:
	    mean_c[j] *= 1.0/len(self.all_chains[1])
	    
	print mean_c
	
	
	radius2 = radius/len(self.all_chains[1]) - (mean_c[0]**2 + mean_c[1]**2 + mean_c[2]**2)
	
	of = open('radius.csv','a')
	
	print >> of, self.current_time, ',', math.sqrt(radius2)
	
	of.close()
	
    def compute_hopping_probs(self):
	self.jump_prob_dict = {}
	
	sampled_coordinates = np.zeros(shape=(self.lattice_size,self.lattice_size,self.lattice_size))
	
	self.jump_prob_dict[str(1.0)] = 0
	
	# Initialize counter
	for k in xrange(0,len(self.fluctuator.max_displacement_energies)):
	    factor = math.exp(-self.fluctuator.max_displacement_energies[k]/self.fluctuator.beta)
	    #factor = -self.fluctuator.max_displacement_energies[k]/self.fluctuator.beta
	    self.jump_prob_dict[str(factor)] = 0
	    
	self.jump_prob_dict[str(0.0)] = 0
	
	for chain_no in self.all_chains.keys():
	    for idset in self.all_chains[chain_no]:
		if self.junction_chain[idset[0],idset[1],idset[2]]<=chain_no and sampled_coordinates[idset[0],idset[1],idset[2]]==0:
		    sampled_coordinates[idset[0],idset[1],idset[2]] = -1

		i = self.get_lidx_from_cubidx(idset)
		    
		for k in xrange(0,6):
		    for key in self.jump_prob_dict.keys():
			if self.moving_probs[i,k]<=float(key):
			    self.jump_prob_dict[key] += 1
	
	# Add contribution due to total monomers
	for key in self.jump_prob_dict.keys():
	    if float(key)<=1.0:
		self.jump_prob_dict[key] += 6*self.all_monomers.total_monomers
		
	max_value = max(self.jump_prob_dict.values())

	for key in self.jump_prob_dict.keys():
	    self.jump_prob_dict[key] *= 1.0/float(max_value)
	    
	ofile = open('hopping_probs'+str(self.proc_no)+'.csv','w')
	
	for key in self.jump_prob_dict.keys():
	    print >> ofile, key+','+str(self.jump_prob_dict[key])
	    
	ofile.close()
	
    def compute_mol_wt_dist(self):
	samples = []
	
	for chain in self.all_chains.keys():
	    #samples.append(len(self.all_chains[chain]))
	    if chain not in self.terminated_pairs.keys():
		samples.append(len(self.all_chains[chain]))
	    else:
		if chain<self.terminated_pairs[chain]:
		    samples.append(len(self.all_chains[chain])+len(self.all_chains[self.terminated_pairs[chain]]))

	m_hist, bin_edges = np.histogram(samples,50,density=True)
	
	ofile = open('molwt_dist.csv','w')
	
	for i in xrange(0,50):
	    print >> ofile, 0.5*(bin_edges[i]+bin_edges[i+1]),',',m_hist[i]

	ofile.close()
	
    def write_final_vtk_files(self):
	filename1 = vtk_writer.case_and_step_file('nodetypes')
	filename2 = vtk_writer.case_and_step_file('composition')
	filename3 = vtk_writer.create_vtk_bond_file(int(self.lattice_size),str(self.vtk_counter))

	ofile1 = open(filename1,'a')
	ofile2 = open(filename2,'a')
	ofile3 = open(filename3,'a')
	box_size = int(self.lattice_size)

	for x in range(0,box_size):
	    for y in range(0,box_size):
		for z in range(0,box_size):
		#    if self.super_lattice[x,y,z]>0:
		#	print >> ofile, '1'
		    if self.Nr[x,y,z]>0:
			print >> ofile1, '5'
		    elif self.junction_chain[x,y,z]>0:
			print >> ofile1, '4'
		    elif self.location_chain[x,y,z]>0:
			print >> ofile1, '3'
		    elif self.check_super_lattice([x,y,z])==0:
			print >> ofile1, '2'
		    elif self.check_super_lattice([x,y,z])<8:
			print >> ofile1, '1'
		    else:
			print >> ofile1, '0'

		    print >> ofile2, str(int(self.location_monomers[x,y,z]+1))

	ofile1.close()
	ofile2.close()
	
	total_line_segments = 0
	
	for chain in self.all_chains:
	    total_line_segments += len(self.all_chains[chain])-1
	
	print >> ofile3, '\n'
	
	#print >> ofile3, 'LINES '+str(total_line_segments)+' '+str(3*total_line_segments)
	
	count = 0
	
	line_sets = []
	
	for chain in self.all_chains:
	    for i in xrange(0,len(self.all_chains[chain])-1):
		idx1, idx2 = self.all_chains[chain][i], self.all_chains[chain][i+1]
		
		if IndexUtils.get_nonperiodic_distance2(idx1,idx2)<=10.0:
		    lidx1, lidx2 = self.get_lidx_from_cubidx(idx1), self.get_lidx_from_cubidx(idx2)
		    
		    outstring = str(int(2))+' '+str(int(lidx1))+' '+str(int(lidx2))
		    
		    line_sets.append(outstring)
		    
	
	total_line_segments = len(line_sets)
	
	print >> ofile3, 'LINES '+str(total_line_segments)+' '+str(3*total_line_segments)
	
	for k in xrange(0,total_line_segments):
	    print >> ofile3 , line_sets[k]
	
	ofile3.close()
	
    def write_chains(self):
	os.system("mkdir "+str(self.proc_no)+"chains")
	os.chdir(str(self.proc_no)+"chains")

	for chain in self.all_chains:
	    ofile = open('chain'+str(chain)+'.csv','w')
	    
	    for node in self.all_chains[chain]:
		print >> ofile, str(node[0]),',',str(node[1]),',',str(node[2])
		
	    ofile.close()
	    
	ofile = open('location_monomers.csv','w')
	
	for x in xrange(0,int(self.lattice_size)):
	    for y in xrange(0,int(self.lattice_size)):
		for z in xrange(0,int(self.lattice_size)):
		    if self.location_monomers[x,y,z]>=0:
			print >> ofile, str(x),',',str(y),',',str(z),',',self.location_monomers[x,y,z]
			
	ofile.close()
	
	ofile = open('lattice.txt','w')
	
	print >> ofile, 'Lattice size = '+str(int(self.lattice_size))
	print >> ofile, 'Site occupancy fraction = '+str(self.all_monomers.occupancy)
	
	ofile.close()
	    
	os.chdir("..")
