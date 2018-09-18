from .aqmmm import AQMMM
from .system import System
import numpy as np

class ONIOM_XS(AQMMM):

    def __init__(self, param, qm_wrapper, mm_wrapper):
        """
        Initializes the ONIOM_XS class object
    
        Parameters
        ----------
        See parameters for AQMMM class 

        Returns
        -------
        A ONIOM_XS class object

        Examples
        --------
        ox = ONIOM_XS(param, psi4_wrapper, openmm_wrapper)
        """
        
        super().__init__(param, qm_wrapper, mm_wrapper)

    def partition(self, qm_center=None, info=None): 
        """
        Finds the partitions as required by the ONIOM-XS method 
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

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm_bz')
            for i, buf in self.buffer_groups.items():
                for idx in buf.atoms:
                    qm_bz.qm_atoms.append(idx)

            # each partition has a copy of its buffer groups - 
            qm_bz.buffer_groups = self.buffer_groups

            self.systems[self.run_ID][qm_bz.partition_ID] = qm_bz

    def run_aqmmm(self):
        """
        Interpolates the energy and gradients from each partition
        according to the ONIOM-XS method
        """
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_forces

        else:
            qm_bz = self.systems[self.run_ID]['qm_bz']
            lamda = self.get_switching_function()
            

            self.systems[self.run_ID]['qmmm_energy'] = \
            (1- lamda)*qm.qmmm_energy + lamda*qm_bz.qmmm_energy

            # needs work!
            # computing gradients
            forces = {}
            for f, coord in qm_bz.qmmm_forces.items():
                if f in qm.qmmm_forces:
                    forces[f] = lamda*coord + (1-lamda)*qm.qmmm_forces[f] 
                else: 
                    forces[f] = lamda*coord

            # computing gradient of switching function
            scaler = (qm_bz.qmmm_energy - qm.qmmm_enery) / len(qm_bz.buffer_groups)

            for i, buf in qm_bz.buffer_groups.items():
                for idx, ratio in buf.weight_ratio.items():
                    forces[idx] += ratio * scaler * buf.d_s_i * buf.COM_coord
                forces[self.qm_center[0]] -= scaler * buf.d_s_i * buf.COM_coord 

            self.systems[self.run_ID]['qmmm_forces'] = forces


    def get_switching_function(self):
        """
        Averages the individual switching functions 
        of each buffer group

        Parameters
        ----------
        None

        Returns
        -------
        The average of the switching functions as a float

        Examples
        --------
        s = get_switching_function()
        """

        s = 0.0
                
        for i, buf in self.buffer_groups.items():
            s += buf.s_i
            
        s *= 1/len(self.buffer_groups)

        return s
            

        
