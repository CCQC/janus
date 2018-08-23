
def test_partition():

    oxs.partition([0])
    oxs_0.partition([0])
    oxs_1.partition([0])
    oxs_2.partition([0])

    assert np.allclose(oxs.partitions['qm'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7,8]))
    assert np.allclose(oxs_0.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_1.partitions['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5]))
    assert np.allclose(oxs_2.partitions['qm'].qm_atoms, np.array([0, 1, 2]))
    assert np.allclose(oxs_2.partitions['qm_bz'].qm_atoms, np.array([0, 1, 2, 3, 4, 5, 6, 7, 8]))


def test_compute_lamda_i():

    s, d = oxs_1.compute_lamda_i(0.333580476259)
    assert (np.allclose(s, 1.1588880014834282) and np.allclose(d, 2.3113804921728383))

def test_switching_function():

    oxs_2.save('qm', np.zeros((9,3)), 1.0)
    oxs_2.save('qm_bz', np.ones((9,3)), 1.2)

    s_1 = oxs_1.get_switching_function(oxs_1.partitions['qm_bz'])
    s_2 = oxs_2.get_switching_function(oxs_2.partitions['qm_bz'])

    assert np.allclose(s_1, 0.9670822167703623)
    assert np.allclose(s_2, 0.862915541444506)

def test_run_aqmmm():

    oxs.save('qm', np.zeros((9,3)), 1.0)
    forces = oxs.get_info()
    
    forces_1 = oxs_1.get_info()
    #forces_2 = oxs_2.get_info()

    assert (np.allclose(forces, np.zeros((9,3))) and oxs.energy == 1.0)
    assert (forces_1 is None and oxs_1.energy == 1.0065835566459276)
    #assert (forces_2 is None and oxs_2.energy == 1.0274168917110988)
