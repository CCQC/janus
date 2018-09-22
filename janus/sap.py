from .aqmmm import AQMMM
from .system import System
import itertools as it
from copy import deepcopy
import numpy as np

class SAP(AQMMM):

    def __init__(self, param, qm_wrapper, mm_wrapper):
        """
        Initializes the SAP class object
    
        Parameters
        ----------
        See parameters for AQMMM class 

        Returns
        -------
        A SAP class object

        Examples
        --------
        sap = SAP(param, psi4_wrapper, openmm_wrapper)
        """
        
        super().__init__(param, qm_wrapper, mm_wrapper, 'SAP')
        self.modified_variant = param['modified_variant']

    def partition(self, qm_center=None, info=None): 
        """
        Finds the partitions as required by the SAP method 
        and saves each partition as a system object.
        Saves all systems in the dictionary self.systems

        Parameters
        ----------
        qm_center: list of atoms that define the qm center, 
                   default is None

        Returns
        -------
        None

        Examples
        --------
        partition([0])
        """
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            self.partitions = self.get_combos(list(self.buffer_groups))

            for i, part in enumerate(self.partitions):
                sys = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID=i)
                for group in part:
                    sys.qm_residues.append(group)
                    for idx in self.buffer_groups[group].atoms:
                        sys.qm_atoms.append(idx)
                
                # each partition has a copy of its buffer groups - 
                # don't know if this is actually needed
                sys.buffer_groups = {k: self.buffer_groups[k] for k in part}
                self.systems[self.run_ID][sys.partition_ID] = sys

    def run_aqmmm(self):
        """
        Interpolates the energy and gradients from each partition
        according to the SAP method
        """
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_forces

        else:

            self.get_switching_functions()

            # getting first term of ap energy and forces (w/o gradient of switching function)
            energy = deepcopy(self.systems[self.run_ID]['qm'].qmmm_energy)
            self.systems[self.run_ID]['qm'].aqmmm_forces = deepcopy(self.systems[self.run_ID]['qm'].qmmm_forces)
            forces = self.systems[self.run_ID]['qm'].aqmmm_forces

            for i, buf in self.buffer_groups.items():
                energy *= (1 - buf.phi_i)
                forces.update((x, y*(1 - buf.phi_i)) for x,y in forces.items())
            self.systems[self.run_ID]['qm'].aqmmm_energy = deepcopy(energy)

            qmmm_forces = deepcopy(forces)

            # getting rest of the terms of sap energy and forces (w/o gradient of switching function)
            for i, part in enumerate(self.partitions):
                part_energy = deepcopy(self.systems[self.run_ID][i].qmmm_energy)
                self.systems[self.run_ID][i].aqmmm_forces = deepcopy(self.systems[self.run_ID][i].qmmm_forces)
                forces = self.systems[self.run_ID][i].aqmmm_forces
                for j, buf in self.buffer_groups.items():
                    if (b_j in part and buf.order == i):
                        part_energy *= buf.phi_i
                        forces.update((x, y*buf.phi_i) for x,y in forces.items())
                    elif b_j not in part:
                        part_energy *= (1 - buf.phi_i)
                        forces.update((x, y*(1 - buf.phi_i)) for x,y in forces.items())

                self.systems[self.run_ID][i].aqmmm_energy = deepcopy(part_energy)
                energy += part_energy

            self.systems[self.run_ID]['qmmm_energy'] = energy

            if self.modified_variant is False:
                # computing forces due to gradient of switching function for SAP
                forces_sf = self.compute_sf_gradient()

                # adding forces together
                for i, force in forces_sf.items():
                    if i in qmmm_forces:
                        qmmm_forces[i] += force
                    else:
                        qmmm_forces[i] = force

            # combining all forces
            for i, part in enumerate(self.partitions):
                forces = self.systems[self.run_ID][i].aqmmm_forces
                for i, force in forces.items():
                    if i in qmmm_forces:
                        qmmm_forces[i] += force
                    else:
                        qmmm_forces[i] = force

            self.systems[self.run_ID]['qmmm_forces'] = forces
            
    def compute_sf_gradient(self):
        """
        Computes forces due to the gradient of the switching function
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Forces due to gradient of switching function as a dictionary

        Examples
        --------
        forces = compute_sf_gradient
        """

        # computing forces due to gradient of switching function for SAP
        forces_sf = {self.qm_center[0]: np.zeros((3))}

        for i, b_i in self.buffer_groups.items():

            # get scaler terms
            b_i.energy_scaler = -1 * self.systems[self.run_ID]['qm'].aqmmm_energy / (1 - b_i.s_i)

            for p, part in enumerate(self.partitions):
                aqmmm_energy = self.systems[self.run_ID][p].aqmmm_energy
                if (i in part and b_i.order == p):
                    b_i.energy_scaler += aqmmm_energy / buf.s_i
                elif i not in part:
                    b_i.energy_scaler -= aqmmm_energy / (1 - buf.s_i)
            
            # get d_phi
            for j, b_j in self.buffer_groups.items():

                force_j = b_i.energy_scaler * b_j.d_s_i * b_j.COM_coord * b_i.d_phi_i[j] * b_i.d_phi_i_scaler

                forces[self.qm_center[0]] -= force_j

                for idx, ratio in b_j.weight_ratio.items():

                    if idx not in forces_sf:
                        forces_sf[idx] = ratio * force_j 
                    else:
                        forces_sf[idx] += ratio * force_j 

        return forces_sf


    def get_combos(self, items=None, buffer_distance=None):
        """
        A given list of indices is sorted by distance and 
        the combinations is found based on distance order.
        See example

        Parameters
        ----------
        items: list of indices to get combinations for
        buffer_distance: dictionary with indices from items 
                         and their distances from the qm center
    
        Returns
        -------
        List of all possible combinations 

        Examples    
        --------
        combos = get_combos([1,2,3])

        If 1 was closer to the center than 2, and 2 closer to the center than 3,
        functions returns following:
        [(1), (1,2), (1,2,3)]
        """
        if buffer_distance is None:
            buffer_distance = self.buffer_distance

        groups = sorted(buffer_distance, key=buffer_distance.get)
        self.sap_order = groups

        all_combo = []
        combo = []

        for g in groups:
            combo.append(g)
            all_combo.append(deepcopy(combo))

        return all_combo


    def get_switching_functions(self):
        """
        Computes switching function for SAP computations
        and saves to each buffer group object
        
        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        get_switching_function()
        """

        sf = self.buffer_groups
    
        for i, b_i in enumerate(self.sap_order):

            sf[b_i].order = i
            sf[b_i].d_phi_i = {}

            chi = (1 - sf[b_i].s_i)/sf[b_i].s_i
            sf[b_i].d_phi_i[b_i] = -1 / sf[b_i].s_i**2

            for j, b_j in enumerate(self.sap_order):
                if j < i:
                    chi += (1 - sf[b_j].s_i)/(sf[b_j].s_i - sf[b_i].s_i)
                    sf[b_i].d_phi_i[b_i] += (1 - sf[b_j].s_i) / (sf[b_j].s_i - sf[b_i].s_i)**2
                    sf[b_i].d_phi_i[b_j]  = (sf[b_i].s_i - 1) / (sf[b_j].s_i - sf[b_i].s_i)**2
                    
                elif j > i:
                    chi += ((1 - sf[b_i].s_i)/(sf[b_i].s_i - sf[b_j].s_i)) * sf[b_j].s_i
                    sf[b_i].d_phi_i[b_i] += (sf[b_j].s_i * (sf[b_j].s_i) - 1) / (sf[b_i].s_i - sf[b_j].s_i)**2
                    sf[b_i].d_phi_i[b_j]  = (sf[b_i].s_i * (1 - sf[b_i].s_i)) / (sf[b_i].s_i - sf[b_j].s_i)**2

            sf[b_i].chi_i = chi
            sf[b_i].phi_i = 1/((1 + chi)**3)
            sf[b_i].d_phi_i_scaler = -3/((1 + chi)**4)

                
