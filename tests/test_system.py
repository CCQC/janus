import pytest
from janus import system
import numpy as np

sys1 = system.System()
sys2 = system.System(qm={'qm_basis_set': '6-31G',
                          'qm_scf_type' : 'pk',
                          'qm_e_convergence' : 1e-9})

def test_build_qm_param():
    
    assert sys1.qm_param['basis'] == 'STO-3G'
    assert sys1.qm_param['scf_type'] == 'df' 
    assert sys1.qm_param['guess'] == 'sad' 
    assert sys1.qm_param['reference'] == 'rhf'
    assert sys1.qm_param['e_convergence'] == 1e-8
    assert sys1.qm_param['d_convergence'] == 1e-8 

    assert sys2.qm_param['basis'] == '6-31G'
    assert sys2.qm_param['scf_type'] == 'pk' 
    assert sys2.qm_param['guess'] == 'sad' 
    assert sys2.qm_param['reference'] == 'rhf'
    assert sys2.qm_param['e_convergence'] == 1e-9
    assert sys2.qm_param['d_convergence'] == 1e-8 

def test_get_link_atom_positions():
    
    qm_position = np.array([0, 3, 4])
    mm_position = np.array([0, 2, -1])
    g = 0.7
    position = sys1.get_link_atom_position(qm_position, mm_position, g)

    assert np.allclose(position, np.array([0, 2.3, 0.5]))

def test_compute_scale_factor_g():
    
    g1 = sys1.compute_scale_factor_g('C', 'C', 'H')
    g2 = sys1.compute_scale_factor_g('C', 'N', 'H')
    
    assert g1 == 0.7133333333333334
    assert g2 == 0.7328767123287672


def test_compute_COM():

    part = system.Partition([0,1,2], 'qm')

    coord = [['O', [0.123, 3.593, 5.841]],
        ['H', [-0.022, 2.679, 5.599]],
        ['H', [0.059, 3.601, 6.796]]]

    xyz = part.compute_COM(coord)
    com = np.array([0.11130575, 3.54230624, 5.88089475])
    
    assert np.allclose(part.qm_atoms, np.array([0,1,2]))
    assert part.ID == 'qm'
    assert not part.switching_functions
    assert np.allclose(xyz, com)

        
