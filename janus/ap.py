from .aqmmm import AQMMM
from .system import System
import itertools as it
from copy import deepcopy
import numpy as np

class AP(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper)

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            self.partitions = self.get_combos(list(self.buffer_groups))

            for i, part in enumerate(self.partitions):
                sys = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID=i)
                for group in part:
                    for idx in self.buffer_groups[group]:
                        sys.qm_atoms.append(idx)
                
                # each partition has a copy of its buffer groups - 
                # don't know if this is actually needed
                sys.buffer_groups = {k: self.buffer_groups[k] for k in part}
                self.systems[self.run_ID][sys.partition_ID] = sys

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_forces

        else:

            if self.aqmmm_scheme == 'PAP': 
                switching_functions = self.buffer_switching_functions
            if self.aqmmm_scheme == 'SAP': 
                switching_functions = self.get_sap_switching_functions()

            # getting first term of ap energy
            energy = self.systems[self.run_ID]['qm'].qmmm_energy
            for buf, func in switching_functions.items():
                energy *= (1 - func[0])

            # getting rest of the terms of ap energy
            for i, part in enumerate(self.partitions):
                part_energy = self.systems[self.run_ID][i].qmmm_energy
                for buf, func in switching_functions.items():
                    if buf in part:
                        part_energy *= func[0]
                    else:
                        part_energy *= (1 - func[0])

                energy += part_energy

            # Need to do gradients! 
                        

    def get_combos(self, items=None, buffer_distance=None):

        if buffer_distance is None:
            buffer_distance = self.buffer_distance
        all_combo = []

        if self.aqmmm_scheme == 'PAP':
            for i in range(1, len(items) +1):
                all_combo += list(it.combinations(items, i))

        if self.aqmmm_scheme == 'SAP':
            groups = sorted(buffer_distance, key=buffer_distance.get)
            self.sap_order = groups
            combo = []
            for g in groups:
                combo.append(g)
                all_combo.append(deepcopy(combo))

        return all_combo


    def get_sap_switching_functions(self):

        switching_function = {}
        sf = self.buffer_switching_functions
    
        for i, b_i in enumerate(self.sap_order):
            chi = (1 - sf[b_i][0])/sf[b_i][0]
            for j, b_j in enumerate(self.sap_order):
                if j < i:
                    chi += (1 - sf[b_j][0])/(sf[b_j][0] - sf[b_i][0])
                elif j > i:
                    chi += ((1 - sf[b_i][0])/(sf[b_i][0] - sf[b_j][0])) * sf[b_j][0]

            switching_function[b_i] = [1/((1 + chi)**3)]
        
        return switching_function


            
