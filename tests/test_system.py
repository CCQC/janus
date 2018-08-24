import pytest
from janus import system
import numpy as np

sys = system.System(qm_indices=[0,1,2], run_ID=0)

def test_compute_scale_factor_g():
    
    g1 = sys.compute_scale_factor_g('C', 'C', 'H')
    g2 = sys.compute_scale_factor_g('C', 'N', 'H')
    
    assert g1 == 0.7133333333333334
    assert g2 == 0.7328767123287672


def test_compute_COM():

    coord = [['O', [0.123, 3.593, 5.841]],
        ['H', [-0.022, 2.679, 5.599]],
        ['H', [0.059, 3.601, 6.796]]]

    xyz = sys.compute_COM(coord)
    com = np.array([0.11130575, 3.54230624, 5.88089475])
    
    assert np.allclose(xyz, com)

        
