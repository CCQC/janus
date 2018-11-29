from janus.qmmm import AQMMM
from janus.system import System
from copy import deepcopy
import numpy as np
import itertools as it
from scipy.misc import logsumexp
from collections import Counter

class DAS(AQMMM):
    """
    Class for the Oniom-XS adaptive QM/MM method.
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

    def __init__(self, hl_wrapper, 
                       ll_wrapper, 
                       sys_info,
                       sys_info_format='pdb',
                       qm_center=[0],
                       partition_scheme='distance',
                       Rmin=3.8,
                       Rmax=4.5,
                       qmmm_param={},
                       **kwargs):

        self.qm_center = qm_center
        self.partition_scheme = partition_scheme
        self.Rmin = Rmin
        self.Rmax = Rmax

        super().__init__(hl_wrapper, ll_wrapper, sys_info, sys_info_format, qmmm_param, 'DAS')

    def partition(self, qm_center=None): 
        """
        Finds the partitions as required by the PAP method 
        and saves each partition as a system object.
        Saves all systems in the dictionary self.systems

        Parameters
        ----------
        qm_center : list 
            atoms that define the qm center, default is None

        """
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm')
        qm.buffer_groups = self.buffer_groups

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            self.partitions, sigmas = self.get_combos(list(self.buffer_groups))
            print('partitions', self.partitions)

            for i, part in enumerate(self.partitions):
                sys = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID=i)
                sys.sigma = sigmas[i]
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
        according to the PAP method
        """
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_forces

        else:

            # getting first term of ap energy and forces (w/o gradient of switching function)
            qm.aqmmm_energy = deepcopy(qm.qmmm_energy)
            qm.aqmmm_forces = deepcopy(qm.qmmm_forces)
            
            dis = sorted(self.buffer_distance, key=self.buffer_distance.get)
            sigma = self.buffer_groups[dis[0]].s_i
            qm.aqmmm_energy *= sigma
            qm.aqmmm_forces.update((x, y*sigma) for x,y in qm.aqmmm_forces.items())

            energy = deepcopy(qm.aqmmm_energy)
            qmmm_forces = deepcopy(qm.aqmmm_forces)

            # getting rest of the terms of ap energy and forces (w/o gradient of switching function)
            for i, part in enumerate(self.partitions):

                sys = self.systems[self.run_ID][i]
                sys.aqmmm_energy = sys.qmmm_energy * sys.sigma

                sys.aqmmm_forces = deepcopy(sys.qmmm_forces)
                sys.aqmmm_forces.update((x, y*sys.sigma) for x,y in sys.aqmmm_forces.items())

                energy += sys.aqmmm_energy

            # combining all forces
            for i, part in enumerate(self.partitions):
                forces = self.systems[self.run_ID][i].aqmmm_forces
                for j, force in forces.items():
                    if j in qmmm_forces:
                        qmmm_forces[j] += force
                    else:
                        qmmm_forces[j] = force


            # need to deal with bookkeeping term
            energy_bk = 0.0
            systems_f = []
            systems_match = []
            if self.run_ID != 0:
                for sys_f in self.systems[self.run_ID]:
                    if isinstance(sys_f, System()):
                        systems_f.append(Counter(sys_f.qm_residues))
                        for sys_i in self.systems[self.run_ID-1]:
                            if isinstance(sys_i, System()):
                                if Counter(sys_f.qm_residues) == Counter(sys_i.qm_residues):
                                    systems_match.append(Counter(sys_f.qm_residues))
                                    energy_bk += sys_f.aqmmm_energy * (sys_f.sigma-sys_i.sigma)

                for sys in systems_f:
                    if sys not in systems_match:
                        energy_bk += sys_f.aqmmm_energy * sys_f.sigma
                        
            energy -= energy_bk
            self.systems[self.run_ID]['qmmm_energy'] = energy
            self.systems[self.run_ID]['qmmm_forces'] = qmmm_forces
            

    def get_combos(self, items=None):
        """
        Gets all combinations of a given list of indices 
        according to the DAS formulation

        Parameters
        ----------
        items : list 
            indices to get combinations for
    
        Returns
        -------
        list 
            combinations

        """
        
        all_combo = []
        edited_combos = []
        sigmas = []

        for i in range(1, len(items) +1):
            all_combo += list(it.combinations(items, i))

        k = 500.0 
        thresh = 0.001
        
        n_qm_res = len(self.qm_residues)
        n_res = self.topology.n_residues

        for combo in all_combo:

            qm_lamda = np.zeros((n_qm_res))
            mm_lamda = np.zeros((n_res - n_qm_res - len(items)))
            for group in items:
                if group in combo:
                    qm_lamda = np.insert(qm_lamda, 0, self.buffer_groups[group].s_i)
                else:
                    mm_lamda = np.insert(mm_lamda, 0, 1-self.buffer_groups[group].s_i)
            
            # makes sure there are elements in mm_lamda
            if not mm_lamda:
                mm_lamda = np.insert(mm_lamda, 0, 0)

            min_mm_lamda = 1 - logsumexp(k*mm_lamda)/k
            max_qm_lamda = logsumexp(k*qm_lamda)/k
            sigma = logsumexp(k*np.array([min_mm_lamda - max_qm_lamda,0]))/k

            if sigma > thresh:
                edited_combos.append(combo)
                sigmas.append(sigma)

        return edited_combos, sigmas

    def compute_lamda_i(self, r_i):
        """
        Computes the switching function of the DAS method
        and overrides the compute_lamda_i function from the AQMMM class
        
        Parameters
        ----------
        r_i : float 
            the distance between the qm center and the COM in angstroms 

        Returns
        -------
        float
            lamda_i, unitless
        None
            placeholder for the derivative of r_i
        
        """

        if r_i <= self.Rmin:
            lamda_i = 0


        elif r_i > self.Rmax:
            lamda_i = 1

        else:
            lamda_i = (r_i - self.Rmin)**2 * (3*self.Rmax - self.Rmin - 2*r_i)
            lamda_i *= 1/(self.Rmax - self.Rmin)**3

        return lamda_i, None
