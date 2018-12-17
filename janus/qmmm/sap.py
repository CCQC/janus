from janus.qmmm import AQMMM
from janus.system import System
import itertools as it
from copy import deepcopy
import numpy as np

class SAP(AQMMM):
    """
    Class for the SAP adaptive QM/MM method.
    Inherits from AQMMM class

    Parameters
    ----------
        hl_wrapper : :class:`~janus.mm_wrapper.MMWrapper` subclass or :class:`~janus.qm_wrapper.QMWrapper` subclass
            Wrapper for performing the high-level computation. 
            Traditionally QM but user can define MM.
        ll_wrapper : :class:`~janus.mm_wrapper.MMWrapper` subclass
            Wrapper for performing the low-level computation
        sys_info : str 
            A string with the filename or a list with multiple filenames 
            that contain position and topology information. 
        sys_info_format : str 
            Describes what kind of input is contained in sys_info. Default is pdb.
        modified_variant : bool
            Whether to use the modified version mSAP, which disregards 
            gradient terms that come from the switching function.
            Default is False.
        qm_center: list 
            Atoms that define the qm center, default is [0].
            If more than one index is given, COM is used as qm_center
        partition_scheme: str 
            Scheme to use to define buffer groups,
            default is distance (only scheme available as of now)
        Rmin: float 
            Inner radius for distance partition in angstroms, default is 3.8
        Rmax: float 
            Outer radius for distance partition in angstroms, default is 4.5
        qmmm_param : dict
            A dictionary with any parameters to pass into the QMMM class.
            See QMMM class for specifics

    """

    def __init__(self, modified_variant=False, *args, **kwargs):

        self.modified_variant = modified_variant
        print(self.modified_variant)
        super().__init__('SAP', *args, **kwargs)

    def find_configurations(self): 
        """
        Finds the partitions as required by the SAP method 
        and saves each partition as a system object.
        Saves all systems in the dictionary self.systems

        Parameters
        ----------
        qm_center : list 
            atoms that define the qm center, default is None

        """

        qm = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm')
        qm.buffer_groups = self.buffer_groups

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            self.partitions = self.get_combos(list(self.buffer_groups))
            print('partitions', self.partitions)

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
            qm.aqmmm_energy = deepcopy(qm.qmmm_energy)
            qm.aqmmm_forces = deepcopy(qm.qmmm_forces)
            print('qm aqmmm energy', qm.aqmmm_energy)

            for i, buf in self.buffer_groups.items():
                qm.aqmmm_energy *= (1 - buf.phi_i)
                qm.aqmmm_forces.update((x, y*(1 - buf.phi_i)) for x,y in qm.aqmmm_forces.items())

            energy = deepcopy(qm.aqmmm_energy)
            qmmm_forces = deepcopy(qm.aqmmm_forces)

            # getting rest of the terms of sap energy and forces (w/o gradient of switching function)
            for i, part in enumerate(self.partitions):

                sys = self.systems[self.run_ID][i]
                sys.aqmmm_energy = deepcopy(sys.qmmm_energy)
                sys.aqmmm_forces = deepcopy(sys.qmmm_forces)

                for j, buf in self.buffer_groups.items():
                    if (j in part and buf.order == i):
                        sys.aqmmm_energy *= buf.phi_i
                        sys.aqmmm_forces.update((x, y*buf.phi_i) for x,y in sys.aqmmm_forces.items())
                    elif j not in part:
                        sys.aqmmm_energy *= (1 - buf.phi_i)
                        sys.aqmmm_forces.update((x, y*(1 - buf.phi_i)) for x,y in sys.aqmmm_forces.items())

                energy += sys.aqmmm_energy

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
                for j, force in forces.items():
                    if j in qmmm_forces:
                        qmmm_forces[j] += force
                    else:
                        qmmm_forces[j] = force

            self.systems[self.run_ID]['qmmm_energy'] = energy
            self.systems[self.run_ID]['qmmm_forces'] = qmmm_forces
            print('forces')
            print(qmmm_forces)
            
    def compute_sf_gradient(self):
        """
        Computes forces due to the gradient of the switching function
        
        Returns
        -------
        dict
            Forces due to gradient of switching function 

        """

        # computing forces due to gradient of switching function for SAP
        forces_sf = {}
        for idx in self.qm_center:
            forces_sf[idx] = np.zeros((3))

        for i, b_i in self.buffer_groups.items():

            # get scaler terms
            b_i.energy_scaler = -1 * self.systems[self.run_ID]['qm'].aqmmm_energy / (1 - b_i.s_i)

            for p, part in enumerate(self.partitions):
                aqmmm_energy = self.systems[self.run_ID][p].aqmmm_energy
                if (i in part and b_i.order == p):
                    b_i.energy_scaler += aqmmm_energy / b_i.s_i
                elif i not in part:
                    b_i.energy_scaler -= aqmmm_energy / (1 - b_i.s_i)
            
            # get d_phi
            for j, b_j in self.buffer_groups.items():

                force_j = b_i.energy_scaler * b_j.d_s_i * b_j.COM_coord * b_i.d_phi_i[j] * b_i.d_phi_i_scaler

                for idx, ratio in self.qm_center_weight_ratio.items():
                    forces_sf[idx] -= ratio * force_j

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
        items : list 
            indices to get combinations for
        buffer_distance: dict 
            indices from items and their distances from the qm center
    
        Returns
        -------
        list 
            combinations

        Examples    
        --------
        >>> combos = get_combos([1,2,3])

        If 1 was closer to the center than 2, and 2 closer to the center than 3,
        functions returns following:
        >>> [(1), (1,2), (1,2,3)]

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

                
