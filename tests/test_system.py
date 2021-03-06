import pytest
from janus import system
import mdtraj as md
import numpy as np
import os

sys = system.System(qm_indices=[0,1,2], qm_residues=[0], run_ID=0)
water = os.path.join(str('tests/files/test_openmm/water.pdb'))
traj = md.load(water)

def test_compute_scale_factor_g():
    
    g1 = system.System.compute_scale_factor_g('C', 'C', 'H')
    g2 = system.System.compute_scale_factor_g('C', 'N', 'H')
    
    assert g1 == 0.7133333333333334
    assert g2 == 0.7328767123287672



        
