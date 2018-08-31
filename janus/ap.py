from .aqmmm import AQMMM
from .system import System
import itertools as it
from copy import deepcopy
import numpy as np

class PAP(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper, 'PAP')

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            partitions = self.get_combos(list(self.buffer_groups))

            for i, part in enumerate(partitions):
                sys = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID=i)
                for group in part:
                    for idx in self.buffer_group[group]:
                        sys.qm_atoms.append(idx)
                
                # each partition has a copy of its buffer groups - 
                sys.buffer_groups = {k: self.buffer_groups[k] for k in part}
                self.systems[self.run_ID][sys.partition_ID] = sys

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_forces

        else:
            pass

    def get_combos(self, items):

        all_combo = []

        for i in range(1, len(items) +1):

            all_combo += list(it.combinations(items, i))

        return all_combo


class SAP(AQMMM):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        super().__init__(config, qm_wrapper, mm_wrapper, 'PAP')

    def partition(self, qm_center=None, info=None): 
    
        if qm_center is None:
            qm_center = self.qm_center

        self.define_buffer_zone(qm_center)

        qm = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID='qm')

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:

            partitions = self.get_all_combos(list(self.buffer_groups))

            for i, part in enumerate(partitions):
                sys = System(qm_indices=self.qm_atoms, run_ID=self.run_ID, partition_ID=i)
                for group in part:
                    for idx in self.buffer_group[group]:
                        sys.qm_atoms.append(idx)
                
                # each partition has a copy of its buffer groups - 
                sys.buffer_groups = {k: self.buffer_groups[k] for k in part}
                self.systems[self.run_ID][sys.partition_ID] = sys

    def run_aqmmm(self):
        
        qm = self.systems[self.run_ID]['qm']

        if not self.buffer_groups:
            self.systems[self.run_ID]['qmmm_forces'] = qm.qmmm_energy
            self.systems[self.run_ID]['qmmm_energy'] = qm.qmmm_forces

        else:
            pass


    def get_combos(self):

        all_combos = [] 
        groups = sorted(self.buffer_distance, key=self.buffer_distance.get)
        combo = []
        for g in groups:
            combo.append(g)
            all_combos.append(deepcopy(mid))

        return all_combos
        
