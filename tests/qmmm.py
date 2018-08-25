
sys_ala = system.System(qm={"qm_atoms": [0,1,2,3,4,5,6]}, mm={'mm_pdb_file' : ala_water_pdb_file})
sys_mech = system.System(qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_elec = system.System(qmmm=qmmm_elec,qm={"qm_atoms": [0,1,2]}, mm={'mm_pdb_file' : water_pdb_file})
sys_ala_link = system.System(qm={"qm_atoms": [0, 1, 2, 3]}, mm={'mm_pdb_file' : ala_pdb_file})

def test_run_qmmm():
    pass

def test_update_traj():
    pass
def test_convert_trajectory():
    pass
def test_mechanical():
    pass
def test_compute_gradients():
    pass

def test_edit_qm_atoms():
    pass
def test_prepare_link_atom():
    pass
def test_make_primary_subsys_trajectory():
    pass
def test_make_second_subsys_trajectory():
    pass
def test_get_forces():
    pass

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
 H     0.887   0.234   0.419 
 """
    
    pos = openmm_mech.get_qm_positions()
    pos_link = openmm_ala_link.get_qm_positions()

    assert pos == qm_mol
    assert pos_link == qm_link_mol

def test_get_redistributed_positions():

    bonds = psi4_rc._system.boundary_info['bonds_to_mm']
    mm = psi4_rc._system.boundary_info['mm_index']
    pos = psi4_rc._system.entire_sys['positions']

    new_pos = np.array([[ 1.60150003,  0.97700004,  0.33600002],
                       [ 1.65199999, -0.19850001,  0.96250005],
                       [ 1.03050001,  0.61799999,  1.06550001]])

    position = psi4_elec.get_redistributed_positions(pos, bonds, mm)
    
    assert np.allclose(position, new_pos)

def test_get_external_charges():

    # test link atom
    charge = psi4_elec.get_external_charges(link=True)
    crg = np.array([[-0.834, 0.115 ,  0.31300001 , 6.14799976],
                        [0.417, 0.12  ,  0.241      , 5.19299984],
                        [0.417, 0.66  , -0.41499998 , 6.44599974]])

    # test RC 
    charge_rc = psi4_rc.get_external_charges(RC=True)

    crg_rc = np.array([[0.0823, 0.81399999558925629, 0.86099997162818909, 1.4949999749660492],
                        [-0.1825, 2.0569999516010284, -0.7720000296831131, 1.2890000641345978], 
                        [0.0603, 3.1360000371932983, -0.75199998915195465, 1.0320000350475311],
                        [0.0603, 1.9900000095367432, -0.64099997282028198, 2.3950000107288361],
                        [0.0603, 1.656000018119812, -1.7820000648498535, 1.062999963760376],
                        [0.5973, 1.956000030040741, 1.5790000557899475, 0.03599999938160181],
                        [-0.5679, 1.2189999967813492, 2.5249999761581421, -0.20099999383091927],
                        [0.011233333333333333, 1.60150003, 0.97700004, 0.33600002],
                        [0.011233333333333333, 1.65199999, -0.19850001, 0.96250005],
                        [0.011233333333333333, 1.03050001, 0.61799999, 1.06550001]])

    # test RC 
    charge_rcd = psi4_rc.get_external_charges(RCD=True)
    crg_rcd = np.array([[0.07106666666666667, 0.81399999999999995, 0.86099996999999995, 1.4949999700000001],
                        [-0.19373333333333334, 2.0569999499999998, -0.77200002999999995, 1.28900006],
                        [0.0603, 3.1360000399999999, -0.75199998999999995, 1.03200004],
                        [0.0603, 1.9900000099999999, -0.64099996999999997, 2.39500001],
                        [0.0603, 1.65600002, -1.7820000600000001, 1.06299996],
                        [0.5860666666666667, 1.95600003, 1.57900006, 0.035999999999999997],
                        [-0.5679, 1.2190000000000001, 2.52499998, -0.20099998999999999],
                        [0.022466666666666666, 1.601500025, 0.97700003499999999, 0.33600002000000001],
                        [0.022466666666666666, 1.6519999849999998, -0.19850000999999998, 0.96250005000000005],
                        [0.022466666666666666, 1.0305000099999999, 0.61799998999999994, 1.0655000050000001]])

    assert np.allclose(np.array(charge), crg)
    assert np.allclose(np.array(charge_rc), crg_rc)
    assert np.allclose(np.array(charge_rcd), crg_rcd)
