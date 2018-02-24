import pytest
from janus import driver
from janus import system
import numpy as np
import os


json_file = os.path.join(str('tests/files/test_driver/input.json'))
water_pdb_file = os.path.join(str('tests/files/test_openmm/water.pdb'))

sys_mech = system.System(qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})


def test_load_system():
    sys = driver.load_system(json_file)

    assert sys.qm_basis_set == 'STO-3G'
    assert sys.qm_scf_type == 'pk' 
    assert sys.qm_guess == 'sad' 
    assert sys.qm_reference == 'rhf'
    assert sys.qm_e_convergence == 1e-9
    assert sys.qm_d_convergence == 1e-8 
    assert sys.qm_method == 'mp2'
    assert sys.qm_atoms == [0,1,2]
    assert sys.qm_residues == []
    assert sys.qm_charge_method == 'MULLIKEN_CHARGES'
    assert sys.mm_pdb_file == "water.pdb"
    assert sys.mm_ff == 'amber99sb.xml'
    assert sys.mm_ff_water == 'tip3p.xml'
    assert sys.mm_step_size == 0.005
    assert sys.embedding_method == 'Electrostatic'
    assert sys.qm_program == 'Psi4'

def test_additive():
    
    sys_mech.boundary['energy'] = -0.010631625405292172
    sys_mech.second_subsys['energy'] = 6.873303688617918e-05
    sys_mech.qm['energy'] = -74.9629737257 
    
    driver.additive(sys_mech)
    assert np.allclose(sys_mech.qmmm_energy, -74.9735366180801)

def test_subtractive():

    sys_mech.entire_sys['energy'] = -0.010562892368405992
    sys_mech.primary_subsys['energy'] = 0.0
    sys_mech.qm['energy'] = -74.9629737257 

    driver.subtractive(sys_mech)
    assert np.allclose(sys_mech.qmmm_energy, -74.9735366180801)
