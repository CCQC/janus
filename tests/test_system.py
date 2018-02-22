import pytest
from janus import system

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

