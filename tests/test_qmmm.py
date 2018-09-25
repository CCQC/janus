from janus import qmmm, psi4_wrapper, openmm_wrapper, system, initializer
import numpy as np
import pytest
import os

water = os.path.join(str('tests/files/test_openmm/water.pdb'))
ala = os.path.join(str('tests/files/test_openmm/ala_ala_ala.pdb'))

param_m = {"system" : {"mm_pdb_file" : water},
            "qmmm" : {"qm_atoms" : [0,1,2]}}
param_ala = {"system" : {"mm_pdb_file" : ala},
              "qmmm" : {"embedding_method" : "Electrostatic"},
              "qm" : {"charge" : +1}}


config_m = initializer.Initializer(param_m, as_file=False)
config_ala = initializer.Initializer(param_ala, as_file=False)

psi4 = psi4_wrapper.Psi4_wrapper(config_m.qm_param)
psi4_ala = psi4_wrapper.Psi4_wrapper(config_ala.qm_param)

om_m = openmm_wrapper.OpenMM_wrapper(config_m.mm_param)
om_m.initialize('Mechanical')
main_info_m = om_m.get_main_info()

om_ala = openmm_wrapper.OpenMM_wrapper(config_ala.mm_param)
om_ala.initialize('Electrostatic')
main_info_ala = om_ala.get_main_info()

mech = qmmm.QMMM(config_m.qmmm_param, psi4, om_m)
ala_link = qmmm.QMMM(config_ala.qmmm_param, psi4_ala, om_ala)
ala_RC = qmmm.QMMM(config_ala.qmmm_param, psi4_ala, om_ala)
ala_RCD = qmmm.QMMM(config_ala.qmmm_param, psi4_ala, om_ala)
ala_RC.boundary_treatment = 'RC'
ala_RCD.boundary_treatment = 'RCD'

sys_mech = system.System([0,1,2], [0],  0)
sys_ala_link = system.System([0,1,2,3],[0],  0)
sys_ala_RC = system.System([0,1,2,3,4,5], [0,1], 0)
sys_ala_RCD = system.System([0], [0],  0)

sys_mech.entire_sys = main_info_m
sys_ala_link.entire_sys = main_info_ala
sys_ala_RC.entire_sys = main_info_ala
sys_ala_RCD.entire_sys = main_info_ala

def test_convert_trajectory():

    top, pos = mech.convert_trajectory(mech.traj)
    assert type(top) is openmm_wrapper.OM_app.Topology
    assert type(pos) is openmm_wrapper.OM_unit.Quantity

def test_edit_qm_atoms():

    a = np.array([0,1,2])
    b = np.array([3,4,5])
    c = np.array([0,1,2,3,4,5])

    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[0,1,2])), a)
    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[2,3])), b)
    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[0,1,2,3])), c)
    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[0,1,2,3,4])), c)
    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[0,1,2,4,5])), a)
    assert np.allclose(np.array(mech.edit_qm_atoms(qm_atoms=[0,1,2,5])), a)

def test_find_boundary_bonds():

    mech.find_boundary_bonds(qm_atoms=[1,2,3,4])
    ala_link.find_boundary_bonds(qm_atoms=[0,1,2,3])
    ala_RC.find_boundary_bonds(qm_atoms=[0,1,2,3,4,5])
    ala_RCD.find_boundary_bonds(qm_atoms=[0])

    assert len(mech.qmmm_boundary_bonds) == 0
    assert len(ala_link.qmmm_boundary_bonds) == 1 
    assert len(ala_RC.qmmm_boundary_bonds) == 2
    assert len(ala_RCD.qmmm_boundary_bonds) == 4

def test_prepare_link_atom():

    ala_link.prepare_link_atom()
    ala_RC.prepare_link_atom()
    ala_RCD.prepare_link_atom()

    assert len(ala_link.link_atoms['all_outer_bonds']) == 0
    assert ala_link.link_atoms[0]['link_atom'] == 'H'
    assert np.allclose(np.array(ala_RC.link_atoms['all_outer_bonds'][1]), np.array([7, 8, 9]))
    assert np.allclose(np.array(ala_RCD.link_atoms['all_outer_bonds'][0]), np.array([10, 6, 5]))
    assert len(np.array(ala_RCD.link_atoms['all_outer_bonds'][1])) == 0

def test_get_redistributed_positions():

    positions = sys_ala_RC.entire_sys['positions']
    pos1 = ala_RC.get_redistributed_positions(positions, ala_RC.link_atoms['all_outer_bonds'][0], ala_RC.link_atoms['all_mm'][0])
    pos2 = ala_RC.get_redistributed_positions(positions, [], ala_RC.link_atoms['all_mm'][0])
    
    pos = np.array([[ 0.15875,  0.2052 , -0.00825],
                    [ 0.26225001, 0.1605, -0.0083]])

    assert len(pos2) == 0
    assert np.allclose(np.array(pos1), pos)

def test_get_external_charges():

    charges_mech = mech.get_external_charges(sys_mech)
    charges_ala_link = ala_link.get_external_charges(sys_ala_link)
    charges_ala_RC = ala_RC.get_external_charges(sys_ala_RC)
    charges_ala_RCD = ala_RCD.get_external_charges(sys_ala_RCD)
    
    assert charges_mech is None
    assert len(charges_ala_link) == 29
    assert len(charges_ala_RC) == 30
    assert len(charges_ala_RCD) == 31

def test_make_primary_subsys_trajectory():

    traj_mech, link_mech = mech.make_primary_subsys_trajectory()
    traj_ala, link_ala = ala_RC.make_primary_subsys_trajectory(qm_atoms=sys_ala_RC.qm_atoms)

    mech.traj_ps = traj_mech
    ala_RC.traj_ps = traj_ala

    assert len(link_mech) == 0
    assert len(link_ala) == 2
    assert len(traj_mech.xyz[0]) == 3
    assert len(traj_ala.xyz[0]) == 8
    
def test_get_qm_geometry():
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
 C     1.247   0.375   0.636 
 H     0.814   0.861   1.495 
 H     1.753   1.234   0.208 
 H     1.825  -0.443   1.102 
 """
    out, elec = mech.get_qm_geometry(mech.traj_ps)
    out2, elec2 = ala_RC.get_qm_geometry(ala_RC.traj_ps)

    assert out == qm_mol
    assert out2 == qm_link_mol
    assert elec == 10
    assert elec2 == 19

   
def test_make_second_subsys_trajectory():

    traj_mech = mech.make_second_subsys_trajectory()
    traj_ala = ala_RC.make_second_subsys_trajectory(qm_atoms=sys_ala_RC.qm_atoms)
    
    assert len(traj_mech.xyz[0]) ==6
    assert len(traj_ala.xyz[0]) == 27


def test_mechanical():
    mech.mechanical(sys_mech, main_info_m)
    ala_link.mechanical(sys_ala_link, main_info_ala)
    
    assert len(sys_mech.qmmm_forces) == 3
    assert len(sys_ala_link.qmmm_forces) == 5
    assert sys_mech.entire_sys['energy'] == -0.010562891563767518
    assert sys_mech.primary_subsys['mm']['energy'] == 0.0
    assert sys_mech.qm_info['energy'] == -74.96297372571573
    assert sys_ala_link.entire_sys['energy'] == 0.4511258275881568
    assert sys_ala_link.primary_subsys['mm']['energy'] == 0.02109665357308523
    assert sys_ala_link.qm_info['energy'] == -55.86812550986576

def test_electrostatic():
    ala_link.electrostatic(sys_ala_link, main_info_ala)
    ala_RCD.electrostatic(sys_ala_RCD, main_info_ala)

    assert len(sys_ala_link.qmmm_forces) == 33
    assert len(sys_ala_RCD.qmmm_forces) == 33
    assert sys_ala_link.entire_sys['energy'] == 0.4511258275881568
    assert sys_ala_link.primary_subsys['mm']['energy'] == 0.02109665357308523
    assert sys_ala_link.second_subsys['mm']['energy'] == -0.04867360210802038
    assert sys_ala_link.qm_info['energy'] == -55.81153697205511
    assert sys_ala_RCD.entire_sys['energy'] == 0.4511258275881568
    assert sys_ala_RCD.primary_subsys['mm']['energy'] == 0.02121477274357972
    assert sys_ala_RCD.second_subsys['mm']['energy'] == -0.04867360210802038
    assert sys_ala_RCD.qm_info['energy'] == -55.83972081949333

def test_update_traj():

    mech.traj.xyz[0] = np.zeros((9,3))
    pos1 = mech.traj.xyz[0]
    mech.update_traj(main_info_m['positions'], main_info_m['topology'])
    
    assert np.allclose(pos1, np.zeros((9,3)))
    assert np.allclose(mech.traj.xyz[0], main_info_m['positions'])
    
def test_run_qmmm():
    mech.qm_atoms = [0,1,2]
    ala_link.qm_atoms = [0,1,2,3]

    mech.run_qmmm(main_info_m)
    ala_link.run_qmmm(main_info_ala)

    assert 'qm' in mech.systems[0]
    assert 'qm' in ala_link.systems[0]
    assert len(mech.systems[0]['qmmm_forces']) == 3
    assert len(ala_link.systems[0]['qmmm_forces']) == 33
    assert mech.systems[0]['qmmm_energy'] == -74.97353661727949
    assert ala_link.systems[0]['qmmm_energy'] == -55.37058880492839
    assert mech.run_ID == 1
    assert ala_link.run_ID == 1
    

def test_get_forces():
    f_mech = mech.get_forces()
    f_link = ala_link.get_forces()

    assert len(f_mech) == 3
    assert len(f_link) == 33


#def test_compute_gradients():
# pass in a system with gradients 
    # tested in mechanical and electrostatic, but need for RCD testing
