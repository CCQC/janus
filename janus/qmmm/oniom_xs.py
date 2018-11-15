from janus.qmmm import AQMMM
from janus.system import System
import numpy as np

class OniomXS(AQMMM):

    def __init__(self, param, hl_wrapper, ll_wrapper, md_simulation_program):
        """
        Initializes the OniomXS class object
    
        Parameters
        ----------
        See parameters for AQMMM class 

        """
        
        super().__init__(param, hl_wrapper, ll_wrapper,md_simulation_program,'ONIOM-XS')

    def partition(self, qm_center=None): 
        """
        Finds the partitions as required by the ONIOM-XS method 
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

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm_bz = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm_bz')
            for i, buf in self.buffer_groups.items():
                qm_bz.qm_residues.append(i)
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
            scaler = (qm_bz.qmmm_energy - qm.qmmm_energy) / len(qm_bz.buffer_groups)

            for i, buf in qm_bz.buffer_groups.items():
                for idx, ratio in buf.weight_ratio.items():
                    forces[idx] += ratio * scaler * buf.d_s_i * buf.COM_coord
                for idx, ratio in self.qm_center_weight_ratio.items():
                    forces[idx] -= ratio * scaler * buf.d_s_i * buf.COM_coord

            self.systems[self.run_ID]['qmmm_forces'] = forces


    def get_switching_function(self):
        """
        Averages the individual switching functions 
        of each buffer group

        Returns
        -------
        float
            The average of the switching functions

        """

        s = 0.0
                
        for i, buf in self.buffer_groups.items():
            s += buf.s_i
            
        s *= 1/len(self.buffer_groups)

        return s
            

        
