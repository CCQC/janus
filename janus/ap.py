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
                    for idx in self.buffer_group[group]:
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

            energy = self.systems[self.run_id]['qm'].qmmm_energy
            for buf, func in switching_functions:
                energy *= (1 - func[0])

            for i, part in enumerate(self.partitions):
                part_energy = self.systems[self.run_id][i].qmmm_energy

                for buf, func in switching_functions:
                    if buf in part:
                        part_energy *= func[0]
                    else:
                        part_energy *= (1 - func[0])

                energy += part_energy

            # Need to do gradients! 
                        

    def get_combos_pap(self, items=None):

        all_combo = []

        if self.aqmmm_scheme == 'PAP':
            for i in range(1, len(items) +1):
                all_combo += list(it.combinations(items, i))

        if self.aqmmm_scheme == 'SAP':
            groups = sorted(self.buffer_distance, key=self.buffer_distance.get)
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
            chi = (1 - sf[b_i])/sf[b_i]
            for j, b_j in enumerate(self.sap_order):
                if j < i:
                    chi += (1 - sf[b_j])/(sf[b_j] - sf[b_i])
                elif j > i:
                    chi += ((1 - sf[b_i])/(sf[b_i] - sf[b_j])) * sf[b_j]

            switching_function[b_i] = chi
        
        return switching_function


            
