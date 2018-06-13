"""
Testing for the openmm_wrapper module
"""
import pytest
from janus import openmm_wrapper 
from janus import system
import numpy as np
import os

ala_water_pdb_file = os.path.join(str('tests/files/test_openmm/ala_water.pdb'))
water_pdb_file = os.path.join(str('tests/files/test_openmm/water.pdb'))
ala_pdb_file = os.path.join(str('tests/files/test_openmm/ala_ala_ala.pdb'))

qmmm_elec = {"embedding_method" : "Electrostatic"}

sys_ala = system.System(qm={"qm_atoms": [0,1,2,3,4,5,6]}, mm={'mm_pdb_file' : ala_water_pdb_file})
sys_mech = system.System(qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_elec = system.System(qmmm=qmmm_elec,qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_ala_link = system.System(qm={"qm_atoms": [0, 1, 2, 3]}, mm={'mm_pdb_file' : ala_pdb_file})

openmm_ala = openmm_wrapper.OpenMM_wrapper(sys_ala)
openmm_mech = openmm_wrapper.OpenMM_wrapper(sys_mech)
openmm_elec = openmm_wrapper.OpenMM_wrapper(sys_elec)
openmm_ala_link = openmm_wrapper.OpenMM_wrapper(sys_ala_link)
openmm_ala_link_2 = openmm_wrapper.OpenMM_wrapper(sys_ala_link)


#@pytest.mark.datafiles('tests/examples/test_openmm/water.pdb')
#def test_get_state_info(datafiles):
#    """
#    Function to test that get_openmm_energy is getting energy correctly.
#    """
#    sys, pdb = create_system(datafiles, 'water.pdb')
#    sim = janus.openmm_wrapper.create_openmm_simulation(sys, 
#                                                        pdb.topology,
#                                                        pdb.positions)
#    state  = janus.openmm_wrapper.get_state_info(sim)
#    energy = state['potential'] + state['kinetic']
#    assert np.allclose(energy, -0.01056289236)
#
#
#
#def test_get_second_subsys():
#def test_get_entire_sys():
#def test_get_boundary():



def test_find_boundary_bonds():

    bonds1 = openmm_ala_link.find_boundary_bonds(qm_atoms=[])
    bonds2 = openmm_ala_link.find_boundary_bonds(qm_atoms=[x for x in range(12,20)])

    assert len(openmm_ala_link._boundary_bonds) == 1
    assert len(bonds1) == 0
    assert len(bonds2) == 2
    assert len(openmm_ala._boundary_bonds) == 4
    assert len(openmm_mech._boundary_bonds) == 0

def test_prepare_link_atom():

    assert openmm_ala_link.link_atoms[0]['qm_id'] == '1'
    assert openmm_ala_link.link_atoms[0]['mm_id'] == '5' 
    assert openmm_ala_link.link_atoms[0]['link_atom'] == 'H'
    assert openmm_ala_link.link_atoms[0]['g_factor'] == 0.7054794520547946 
    assert np.allclose(openmm_ala_link.link_atoms[0]['link_positions'], np.array([0.08868014, 0.02342192, 0.04189384]))
    

def test_get_entire_sys():
    openmm_mech._entire_sys['energy'] = -0.010562892368405992
    info = openmm_mech.get_entire_sys()
    assert np.allclose(info['energy'], -0.010562892368405992)

def test_get_second_subsys():
    info = openmm_mech.get_second_subsys()
    assert np.allclose(info['energy'], 6.873303688617918e-05)
    assert np.allclose(info['charges'], np.array([-0.834, 0.417, 0.417, -0.834, 0.417, 0.417]))

def test_get_primary_subsys():
    openmm_mech._primary_subsys['energy'] = 0.0
    info = openmm_mech.get_primary_subsys()

    assert np.allclose(info['energy'], 0.0)

def test_primary_subsys_info():

    # only testing link=False for now
    positions = np.array([[ 0.0024    , -0.0103    , -0.0101    ],
                          [ 0.0027    , -0.1132    , -0.0239    ],
                          [-0.0805    ,  0.0163    ,  0.0471    ],
                          [-0.0059    ,  0.0384    , -0.1019    ]])

    forces = np.array([[-1478.41748047, -367.96447754, -777.9967041 ],
                       [  528.53955078,  -31.95013428,  249.72685242],
                       [  374.38952637,  175.23849487,  461.11102295],
                       [  575.48840332,  224.6761322 ,   67.1587677 ]])

    openmm_ala_link_2.primary_subsys_info(link=False, coulomb=True)

    assert np.allclose(openmm_ala_link_2._primary_subsys['energy'], 0.008915438339083044)
    assert np.allclose(openmm_ala_link_2._primary_subsys['positions'], positions * 10)
    assert np.allclose(openmm_ala_link_2._primary_subsys['forces'], forces) 
    
def test_create_new_residue_template():
    mm = openmm_ala_link.create_modeller(keep_qm = False)
    openmm_ala_link.create_new_residue_template(mm.topology)
    
    assert  openmm_ala_link._ff._templates['Modified_ALA'].name == 'Modified_ALA'

def test_get_boundary():
    openmm_mech._boundary['primary_subsys']['energy'] = 0.0
    openmm_mech._boundary['second_subsys']['energy'] = 6.873303688617918e-05
    openmm_mech._boundary['entire_sys']['energy'] = -0.010562892368405992
    info = openmm_mech.get_boundary()
    assert np.allclose(info['energy'], -0.010631625405292172)

def test_get_boundary_elec():
    openmm_elec._boundary['primary_subsys']['energy'] = 0.0
    openmm_elec._boundary['second_subsys']['energy'] = -1.6972570321173706e-05 
    openmm_elec._boundary['entire_sys']['energy'] = -0.00015962087002150763 
    info = openmm_elec.get_boundary()
    assert np.allclose(info['energy'], -0.00014264829970033393)

def test_get_qm_positions():
    """
    Function to test get_qmmm_energy function
    of systems class given the mm energy and the mm energy
    of the qm region
    """
    qm_mol = """O     0.123   3.593   5.841 
 H    -0.022   2.679   5.599 
 H     0.059   3.601   6.796 
 """

    qm_link_mol = """N     0.024  -0.103  -0.101 
 H     0.027  -1.132  -0.239 
 H    -0.805   0.163   0.471 
 H    -0.059   0.384  -1.019 
 H     0.887   0.234   0.419 
 """
    
    pos = openmm_mech.get_qm_positions()
    pos_link = openmm_ala_link.get_qm_positions()

    assert pos == qm_mol
    assert pos_link == qm_link_mol

def test_keep_residues():
    """
    Function to test keep_residues function.
    """
    mod_str = openmm_ala.create_modeller()
    mod_int = openmm_ala.create_modeller()
    openmm_wrapper.OpenMM_wrapper.keep_residues(mod_str, ['ALA'])
    openmm_wrapper.OpenMM_wrapper.keep_residues(mod_int, [0,1,2])
    res_str = mod_str.topology.getNumResidues()
    res_int = mod_int.topology.getNumResidues()
    assert res_str == 3 
    assert res_int == 3 

def test_keep_atoms():
    """
    Function to test keep_atoms function.
    """
    qm_atm = openmm_ala._system.qm_atoms 
    mod_str = openmm_ala.create_modeller()
    mod_int = openmm_ala.create_modeller()
    openmm_wrapper.OpenMM_wrapper.keep_atoms(mod_int, qm_atm)
    openmm_wrapper.OpenMM_wrapper.keep_atoms(mod_str, ['N'])
    atom_int = mod_int.topology.getNumAtoms()
    atom_str = mod_str.topology.getNumAtoms()
    assert atom_int == len(qm_atm)
    assert atom_str == 3

def test_delete_residues():
    """
    Function to test delete_residues function.
    """
    mod_str = openmm_ala.create_modeller()
    mod_int = openmm_ala.create_modeller()
    openmm_wrapper.OpenMM_wrapper.delete_residues(mod_str, ['ALA'])
    openmm_wrapper.OpenMM_wrapper.delete_residues(mod_int, [0,1,2])
    res_str = mod_str.topology.getNumResidues()
    res_int = mod_int.topology.getNumResidues()
    assert res_str == 28 
    assert res_int == 28 


def test_delete_atoms():
    """
    Function to test delete_atoms function.
    """
    qm_atm = openmm_ala._system.qm_atoms 
    mod_str = openmm_ala.create_modeller()
    mod_int = openmm_ala.create_modeller()
    openmm_wrapper.OpenMM_wrapper.delete_atoms(mod_int, qm_atm)
    openmm_wrapper.OpenMM_wrapper.delete_atoms(mod_str, ['N'])
    atom_int = mod_int.topology.getNumAtoms()
    atom_str = mod_str.topology.getNumAtoms()
    assert atom_int == 117 - len(qm_atm)
    assert atom_str == 114


def test_delete_water():
    """
    Function to test delete_water function.
    """
    mod = openmm_ala.create_modeller()
    openmm_wrapper.OpenMM_wrapper.delete_water(mod)
    res = mod.topology.getNumResidues()
    atom = mod.topology.getNumAtoms()
    assert res == 3
    assert atom == 33 
