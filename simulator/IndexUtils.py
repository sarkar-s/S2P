# Module to permute the list
import itertools
import random as rand
import numpy as np
import math
import sys

maxIdx = [None]
all_neighbors = []
secondary_neighbors = []

def make_neighbor_set():
    iSet = [1,0,-1]

    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                all_neighbors.append([idX,idY,idZ])
    
    all_neighbors.remove([0,0,0])
    
    iSet = range(-2,3)
    
    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                secondary_neighbors.append([idX,idY,idZ])
    
    for site in all_neighbors:
        secondary_neighbors.remove(site)

make_neighbor_set()

# Permutation function to create the set of bond step directions
def permuteList(idx_set):
    all_list = list(itertools.permutations(idx_set))
        
    idx_list = []
        
    for this_list in all_list:
        if list(this_list) not in idx_list:
            idx_list.append(list(this_list))
            # idx_list.append(list(tuple(-1*x for x in this_list)))
    
    return idx_list

# Create the set of bond step directions
totalList = permuteList([2,0,0]) + permuteList([-2,0,0])
totalList = totalList + permuteList([2,1,0]) + permuteList([-2,1,0])
totalList = totalList + permuteList([-2,-1,0]) + permuteList([2,-1,0])
totalList = totalList + permuteList([2,1,1]) + permuteList([-2,1,1])
totalList = totalList + permuteList([-2,-1,1]) + permuteList([2,-1,1])
totalList = totalList + permuteList([2,-1,-1]) + permuteList([-2,-1,-1]) 
totalList = totalList + permuteList([2,2,1]) + permuteList([-2,2,1])
totalList = totalList + permuteList([-2,2,-1]) + permuteList([2,2,-1])
totalList = totalList + permuteList([-2,-2,1]) + permuteList([-2,-2,-1])
totalList = totalList + permuteList([3,0,0]) + permuteList([-3,0,0])
totalList = totalList + permuteList([3,1,0]) + permuteList([-3,1,0])
totalList = totalList + permuteList([-3,-1,0]) + permuteList([3,-1,0])

length_set = [4.0,5.0,6.0,9.0,10.0]

maxSites = len(totalList)
site_range = range(0,maxSites)

cardinal_set = permuteList([1,0,0]) + permuteList([-1,0,0])
# Dict for grouping bonds along the direction of orientation
arranged_bond_set = {}
# Dict with 6 bonds one towards each cardinal direction
bond_groups = {}

for i in xrange(0,6):
    arranged_bond_set[i] = []

def unique_max(cos_set):
    max_cos = max(cos_set)
    c = cos_set.count(max_cos)
    
    if c>1:
        idx = 0
        max_v = -1
        max_idx = 0
        
        for k in xrange(0,6):
            if cos_set[k]!=max_cos and cos_set[k]>max_v:
                max_v = cos_set[k]
                max_idx = k
    else:
        max_idx = cos_set.index(max_cos)
        
    return max_idx

def classify_bonds():
    for i in xrange(0,maxSites):
        max_cos = 0.0
        max_cos_idx = 0
        
        cos_set = []
    
        for j in xrange(0,6):
            cos_set.append(np.inner(totalList[i],cardinal_set[j]))
            
        max_cos_idx = unique_max(cos_set)
            
        arranged_bond_set[max_cos_idx].append(i)
    
    size = len(arranged_bond_set[0])
    
    for idx in xrange(0,size):
        bond_groups[idx] = []
        
        for k in xrange(0,6):
            bond_groups[idx].append(arranged_bond_set[k][idx])

classify_bonds()
group_indices = range(0,len(bond_groups.keys()))

directional_samples = maxSites/6
total_bond_samples = directional_samples*6

def sample_bond_sites():
    samples = []
    rand.shuffle(group_indices)
    
    for i in group_indices:
        #this_set = rand.sample(arranged_bond_set[i],directional_samples)
        rand.shuffle(bond_groups[i])
        
        samples = samples + bond_groups[i]
        
        #for idx in this_set:
        #    samples.append(idx)
    
    return samples
    
low_en_list = permuteList([3,0,0]) + permuteList([-3,0,0])

low_en_grp = []

for idx in low_en_list:
    low_en_grp.append(totalList.index(idx))
    
low_en_samples = len(low_en_grp)
    
def sample_low_en_sites():
    samples =[]
    
    rand.shuffle(low_en_grp)
    
    samples = list(low_en_grp)
    
    return samples

high_en_groups = {}

for j in group_indices:
    high_en_groups[j] = list(bond_groups[j])
    
    for idx in low_en_list:
        if idx in high_en_groups[j]:
            high_en_groups[j].remove(idx)
            
low_en_samples = len(low_en_grp) 

def sample_high_en_sites():
    samples = []
    rand.shuffle(group_indices)
    
    for i in group_indices:
        #this_set = rand.sample(arranged_bond_set[i],directional_samples)
        rand.shuffle(high_en_groups[i])
        
        samples = samples + high_en_groups[i]
        
        #for idx in this_set:
        #    samples.append(idx)
    
    return samples

#low_en_range = range(0,len(low_en_list))

# Compute and store mean bond length
bsum = 0.0
for bset in totalList:
    bsum += (bset[0]**2 + bset[1]**2 + bset[2]**2)**(0.5)
    
meanb = bsum/maxSites
    
lmin = 2.0
lmax = (10.0)**(0.5)

# Thermal moves
thermalMoves = permuteList([1,0,0]) + permuteList([0,0,-1])

# Get thermal move indices in the neighbor list
thermal_idx = []

for move in thermalMoves:
    thermal_idx.append(all_neighbors.index(move))

maxmoves = len(thermalMoves)  
move_nos = range(0,maxmoves) 

# Location for blocking monomers
blockers = []

def make_blocker_set():
    iSet = [1,0,-1]
    temp_set = []
    
    for idX in iSet:
        for idY in iSet:
            temp_set.append([idX,idY])
            
    for k in [0,1,2]:   
        for idset in temp_set:
            this_2 = [idset[0],idset[1]]
            this_2.insert(k,-2)
            blockers.append(this_2)
            
    for k in [2,1,0]:
        for idset in temp_set:
            this_1 = [idset[0],idset[1]]
            this_1.insert(k,2)
            blockers.append(this_1)

make_blocker_set()

blockers_size = len(blockers)
blocker_range = range(0,blockers_size)

# List of influenced monomers
influenced_sites = []
influenced_2D_sites = []

# Spring vectors for the influenced sites
stiff_vecs = []
stiff_norms = []

lumped_mass_neighbors = []
short_mass_neighbors = []

def get_stiffness_components(blk_idx):
    norm = math.sqrt(blk_idx[0]**2 + blk_idx[1]**2 + blk_idx[2]**2)
    #k_ij = compute_stiffness(norm)

    stiff_vec = []
    
    for idx in blk_idx:
        if norm!=0.0:
            stiffness_factor = 4.0/(norm**2)
            value = stiffness_factor*float(idx)/norm
            
            stiff_vec.append(stiffness_factor*float(idx)/norm)
        else:
            stiff_vec.append(0.0)
    
    return stiff_vec

def make_lumped_mass_neighbors():
    iSet = range(-3,4)
    
    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                lumped_mass_neighbors.append([idX,idY,idZ])
            
    lumped_mass_neighbors.remove([0,0,0])
    
    for idx in all_neighbors:
        lumped_mass_neighbors.remove(idx)
                
    for idset in lumped_mass_neighbors:
        this_stiff = get_stiffness_components(idset)
        
        stiff_vecs.append(this_stiff)
        #stiff_norms.append(this_norm)
        
def make_short_mass_neighbors():
    iSet = range(-2,3)
    
    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                lumped_mass_neighbors.append([idX,idY,idZ])
            
    lumped_mass_neighbors.remove([0,0,0])
    
    for idx in all_neighbors:
        lumped_mass_neighbors.remove(idx)
                
make_lumped_mass_neighbors()
make_short_mass_neighbors()

mass_neigbors_set = range(0,len(lumped_mass_neighbors))

# Compute stiffness for this site
def compute_stiffness(r):
    if r!=0.0:
        # phi_prime_prime contribution
        k_ij = (1.0/(r**2))*(156*(2.0/r)**12 - 42*(2.0/r)**6)
        k_ij = (1.0/(r**2))*(-12*(2.0/r)**12 + 6*(2.0/r)**6)
        
        k_ij *= (4.0/3.0)
    else:
        k_ij = 0.0
    
    return k_ij

def make_infuenced_set():
    #iSet = [-4,-3,-2,-1,0,1,2,3,4]
    iSet = range(-2,3)
    
    for idX in iSet:
        for idY in iSet:
            influenced_2D_sites.append([idX,idY])
            for idZ in iSet:
                influenced_sites.append([idX,idY,idZ])

    #for idset in all_neighbors:
    #    influenced_sites.remove(idset)
        
    #influenced_sites.remove([0,0,0])
    
    for idset in influenced_sites:
        this_stiff = get_stiffness_components(idset)
        
        stiff_vecs.append(this_stiff)
        #stiff_norms.append(this_norm)
        
    for idx in all_neighbors:
        influenced_sites.remove(idx)

make_infuenced_set()
full_sites_size = len(influenced_sites)

interaction_set = []

def make_interaction_set(inter_range):
    iSet = range(-inter_range,inter_range+1)
    
    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                interaction_set.append([idX,idY,idZ])
                
    interaction_set.remove([0,0,0])

# These are the list of lattice sites around the center node [0,0,0]
# That would be empty if all the bonds are 
density_count_sites = []

def make_density_count_set():
    iSet = [-2,-1,0,1,2]
    
    for idX in iSet:
        for idY in iSet:
            for idZ in iSet:
                density_count_sites.append([idX,idY,idZ])

    for idset in all_neighbors:
        density_count_sites.remove(idset)
        
    density_count_sites.remove([0,0,0])

make_density_count_set()

full_sites_size = len(density_count_sites)

# Create the set of bond angles
all_angles = []

def create_bond_angles():    
    for vector_1 in totalList:
        for vector_2 in totalList:
            d_prod = np.inner(vector_1,vector_2)
            lv1 = np.linalg.norm(vector_1)
            lv2 = np.linalg.norm(vector_2)
            
            cos_theta = d_prod/(lv1*lv2)

            if abs(cos_theta)>1.0:
                cos_theta += -np.sign(cos_theta)*np.finfo(float).eps 
            
            theta = math.acos(cos_theta)
            
            if theta not in all_angles:
                all_angles.append(theta)
                
create_bond_angles()

def get_max_index(_maxIdx):
    maxIdx[0] = _maxIdx
    
# Get indices for a given box number
def nearest_neighbors(boxidx):
    allidxs = []

    for idxs in all_neighbors:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId == maxIdx[0]:
                thisId = 0
            elif thisId < 0:
                thisId = maxIdx[0] - 1 
        
            this_set.append(thisId)

        allidxs.append(this_set)

    return allidxs

# Get indices for a given box number
def nearest_effected(boxidx):
    allidxs = []
    allidxs.append(boxidx)

    for idxs in interaction_set: 
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId
        
            this_set.append(thisId)

        allidxs.append(this_set)
    
    #for idxs in blockers:
    #    this_set = []
    #    for k in range(0,3):
    #        thisId = boxidx[k]+idxs[k]
    #        # If index number is the last domain
    #        if thisId >= maxIdx[0]:
    #            thisId = thisId - maxIdx[0]
    #        elif thisId < 0:
    #            thisId = maxIdx[0] + thisId
    #    
    #        this_set.append(thisId)
    #
    #    allidxs.append(this_set)

    return allidxs

# Get indices for a given box number
def second_nearest_neighbors(boxidx):
    allidxs = []

    for idxs in secondary_neighbors:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId == maxIdx[0]:
                thisId = 0
            elif thisId < 0:
                thisId = maxIdx[0] - 1 
        
            this_set.append(thisId)

        allidxs.append(this_set)

    return allidxs

# Get the next move index numbers in the lattice
def sitesfornextmove(point_idx):
    global_list = []
        
    for this_idx in totalList:
        this_list = [0,0,0]
        for i in range(0,3):
            this_list[i] = point_idx[i] + this_idx[i]

            if this_list[i] >= maxIdx[0]:
                this_list[i] = this_list[i] - maxIdx[0]
            elif this_list[i] < 0:
                this_list[i] = maxIdx[0] + this_list[i]
        
        global_list.append(this_list)
            
    return global_list

# Get the next move index numbers in the lattice
def nextbondsites(point_idx,bonds):
    global_list = []
        
    for bond in bonds:
        this_idx = totalList[bond]
        this_list = [0,0,0]
        for i in range(0,3):
            this_list[i] = point_idx[i] + this_idx[i]

            if this_list[i] >= maxIdx[0]:
                this_list[i] = this_list[i] - maxIdx[0]
            elif this_list[i] < 0:
                this_list[i] = maxIdx[0] + this_list[i]
        
        global_list.append(this_list)
            
    return global_list

# Get the next move index numbers in the lattice
def sites_for_thermal_moves(point_idx):
    global_list = []
        
    for this_idx in thermalMoves:
        this_list = [0,0,0]
        for i in range(0,3):
            this_list[i] = point_idx[i] + this_idx[i]

            if this_list[i] >= maxIdx[0]:
                this_list[i] = this_list[i] - maxIdx[0]
            elif this_list[i] < 0:
                this_list[i] = maxIdx[0] + this_list[i]
        
        global_list.append(this_list)
            
    return global_list

# Get low energy sites for next move
def low_sitesfornextmove(point_idx):
    global_list = []
        
    for this_idx in low_en_list:
        this_list = [0,0,0]
        for i in range(0,3):
            this_list[i] = point_idx[i] + this_idx[i]

            if this_list[i] >= maxIdx[0]:
                this_list[i] = this_list[i] - maxIdx[0]
            elif this_list[i] < 0:
                this_list[i] = maxIdx[0] + this_list[i]
        
        global_list.append(this_list)
            
    return global_list

# Select a move for random motion
def getthermalmove():
    moveno = rand.randint(0,5)

    return thermalMoves[moveno]

# Compute distance between two boxes
def get_distance(idx1,idx2):
    total = 0.0
    
    for i in range(0,3):
        # Shortest distance consider periodicity
        short1 = abs(idx1[i] - idx2[i])
        short2 = abs(abs(idx1[i] - idx2[i]) - maxIdx[0])
        
        if short1 < short2:
            total += short1**2
        else:
            total += short2**2
        
    dist = total**(0.5)
    
    return dist

# Compute distance between two boxes
def get_unit_vector(idx1,idx2):
    distance = get_distance(idx1,idx2)
    
    unit_vec = []
    
    for k in range(0,3):
        thisId = idx2[k] - idx1[k]
        # If index number is the last domain
        if thisId > 3:
            thisId = thisId - maxIdx[0]
        elif thisId < -3:
            thisId = maxIdx[0] + thisId 
        
        unit_vec.append(float(thisId)/distance)
    
    return unit_vec, distance

# Compute distance between two boxes
def get_distance2(idx1,idx2):
    total = 0.0
    
    for i in range(0,3):
        # Shortest distance consider periodicity
        short1 = abs(idx1[i] - idx2[i])
        short2 = abs(abs(idx1[i] - idx2[i]) - maxIdx[0])
        
        if short1 < short2:
            total += short1**2
        else:
            total += short2**2
    
    return total

def get_nonperiodic_distance2(idx1,idx2):
    total = 0.0
    
    for i in range(0,3):
        # Shortest distance consider periodicity
        short1 = abs(idx1[i] - idx2[i])

        total += short1**2

    return total

# Check bond motion
def check_motion(idx0,idx1,idx2):
    if idx1 == None:
        dist1 = 4.0
    else:
        dist1 = get_distance2(idx0,idx1)

    if idx2 == None:
        dist2 = 4.0
    else:
        dist2 = get_distance2(idx0,idx2)
        
    check = 'good move'

    if dist1 not in length_set:
        check = 'bad move'

    if dist2 not in length_set:
        check = 'bad move'

    return check, dist1, dist2

# Get the set of indices in superlattice surrounding a monomer
def get_super_lattice_pts(idx):
    allIdx = []
    
    for ix in [0,1]:
        x = idx[0]+ix
        
        if x>=maxIdx[0]:
            x = 0
            
        for iy in [0,1]:
            y = idx[1]+iy
            
            if y>=maxIdx[0]:
                y = 0
                
            for iz in [0,1]:
                z = idx[2]+iz
                
                if z>=maxIdx[0]:
                    z = 0
                
                allIdx.append([x,y,z])
                
    return allIdx
    
# Get super lattice point sets for movability
def get_blocking_monomers(boxidx,location_chain):
    allidxs = []
    node_idx = []
    count = 0
    
    for idxs in blockers:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
            
        if location_chain[this_set[0],this_set[1],this_set[2]]>0:
            allidxs.append(this_set)
            node_idx.append(count)
            
        count += 1

    return allidxs, node_idx

# Get super lattice point sets for movability
def get_blocking_sites(boxidx):
    allidxs = []
    
    for idxs in blockers:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
            
        allidxs.append(this_set)

    return allidxs

# Get number of monomers blocking each side
def side_blockers_count(boxidx,location_chain):
    allidxs = []
    node_idx = []
    count = 0
    
    side_count = [0,0,0,0,0,0]
    
    neighbor_stiffness = [0,0,0]
    
    for idxs in blockers:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
            
        stiff_vec, norm = get_stiffness_components(idxs)
            
        if location_chain[this_set[0],this_set[1],this_set[2]]>0:
            factor = 1.0

            for i in range(0,3):
                neighbor_stiffness[i] += abs(stiff_vec[i])/(norm**2)
            
            #side_count[(count/9)] += 1
            
        count += 1
    
    return neighbor_stiffness#side_count

# Get number of monomers blocking each side
def local_stiffness_matrix(boxidx,location_chain):
    allidxs = []
    node_idx = []
    count = 0
    
    side_count = [0,0,0,0,0,0]
    
    neighbor_stiffness = [0,0,0,0,0,0]
    
    neighbor_sites = []
    
    for id_no in xrange(0,full_sites_size):
        idxs = influenced_sites[id_no]
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
 
        if location_chain[this_set[0],this_set[1],this_set[2]]>0:
            this_stiff = stiff_vecs[id_no]
            #this_norm = stiff_norms[id_no]
            
            for loc_id in xrange(0,3):
                if this_stiff[loc_id]>0:
                    neighbor_stiffness[loc_id] += this_stiff[loc_id]
                elif this_stiff[loc_id]<0:
                    neighbor_stiffness[loc_id+3] += abs(this_stiff[loc_id])

        count += 1
    
    return neighbor_stiffness

def get_mass_location_sites(boxidx):
    allidxs = []
    node_idx = []
    count = 0
    
    side_count = [0,0,0,0,0,0]
    
    neighbor_stiffness = [0,0,0,0,0,0]
    
    neighbor_sites = []
    
    for idxs in density_count_sites:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
        
        neighbor_sites.append(this_set)

    return neighbor_sites

# Get neighboring nodes influenced by an update
def get_influenced_monomers(boxidx,location_chain):
    allidxs = []

    for idxs in influenced_sites:
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)
        
        if location_chain[this_set[0],this_set[1],this_set[2]]>0:
            allidxs.append(this_set)

    return allidxs

def get_influenced_monomers_thermal(boxidx,move,location_chain):
    allidxs = get_influenced_monomers(boxidx,location_chain)
    
    for k in [0,1,2]:
        if move[k]!=0:
            set_idx = k
            set_sign = -move[k]
            
    for k in xrange(0,len(influenced_2D_sites)):
        this_idx = list(influenced_2D_sites[k])
        this_idx.insert(set_idx,set_sign*3)
    
        this_set = []
        for k in range(0,3):
            thisId = boxidx[k]+this_idx[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)

        if location_chain[this_set[0],this_set[1],this_set[2]]>0:
            allidxs.append(this_set)

    return allidxs

def get_lattice_set_from_super_lattice(super_idx):
    lattice_idxs = []
    
    for x in [super_idx[0],super_idx[0]-1]:
        # If index number is the last domain
        if x >= maxIdx[0]:
            x = 0
        elif x < 0:
            x = maxIdx[0] - 1
            
        for y in [super_idx[1],super_idx[1]-1]:
            # If index number is the last domain
            if y >= maxIdx[0]:
                y = 0
            elif y < 0:
                y = maxIdx[0] - 1
                
            for z in [super_idx[2],super_idx[2]-1]:
                # If index number is the last domain
                if z >= maxIdx[0]:
                    z = 0
                elif y < 0:
                    z = maxIdx[0] - 1

                lattice_idxs.append([x,y,z])
                
    return lattice_idxs

def get_move_type(moving_idx,fix_idx):
    # Find out which type of thermal move is blocked for
    # the movable site w.r.t. a trial site
    #print moving_idx, fix_idx, maxIdx[0]
    
    move = [0,0]
    
    for k in [0,1,2]:
        diff = fix_idx[k] - moving_idx[k]
        
        if diff>2:
            diff = diff - maxIdx[0]
        elif diff<-2:
            diff = diff + maxIdx[0]
            
        if abs(diff)>1:
            move.insert(k,diff/2)
    #print move
    move_no = thermalMoves.index(move)
    
    #moved_idx = []
    #
    #for k in range(0,3):
    #    thisId = moving_idx[k] + move[k]
    #    # If index number is the last domain
    #    if thisId >= maxIdx[0]:
    #        thisId = thisId - maxIdx[0]
    #    elif thisId < 0:
    #        thisId = maxIdx[0] + thisId
    #    
    #    moved_idx.append(thisId)
    
    return move_no#, moved_idx

def get_moved_location(moving_idx,move_no):
    move = thermalMoves[move_no]
    
    moved_idx = []
    
    for k in range(0,3):
        thisId = moving_idx[k] + move[k]
        # If index number is the last domain
        if thisId >= maxIdx[0]:
            thisId = thisId - maxIdx[0]
        elif thisId < 0:
            thisId = maxIdx[0] + thisId
        
        moved_idx.append(thisId)
        
    return moved_idx

def get_super_lattice_pts_periodic(idx):
    allIdx = []
    
    for ix in [0,1]:
        x = idx[0] + ix
        
        if x >= maxIdx[0]:
            x = 0
            
        for iy in [0,1]:
            y = idx[1] + iy
            
            if y >= maxIdx[0]:
                y = 0
                
            for iz in [0,1]:
                z = idx[2] + iz
                
                if z >= maxIdx[0]:
                    z = 0
                allIdx.append([x,y,z])
                
    return allIdx

bl_move_type = []

for blocker in blockers:
    bl_move_type.append(get_move_type(blocker,[0,0,0]))
    
def get_mass_neighbors(idx):
    allidxs = []
    
    for idxs in lumped_mass_neighbors:
        this_set = []
        for k in range(0,3):
            thisId = idx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)

        allidxs.append(this_set)

    return allidxs

def get_short_mass_neighbors(idx):
    allidxs = []
    
    for idxs in short_mass_neighbors:
        this_set = []
        for k in range(0,3):
            thisId = idx[k]+idxs[k]
            # If index number is the last domain
            if thisId >= maxIdx[0]:
                thisId = thisId - maxIdx[0]
            elif thisId < 0:
                thisId = maxIdx[0] + thisId 
        
            this_set.append(thisId)

        allidxs.append(this_set)

    return allidxs