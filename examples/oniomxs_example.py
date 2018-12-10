from janus import qm_wrapper, mm_wrapper, qmmm

# instantiate a Psi4Wrapper object as the high level wrapper
hl_wrapper = qm_wrapper.Psi4Wrapper()

# instantiate an OpenMMWrapper object as the low level wrapper
ll_wrapper = mm_wrapper.OpenMMWrapper(sys_info='water.pdb')

# instantiate an OniomXS object, setting Rmin and Rmax to be different values
p1 = qmmm.PAP(hl_wrapper, ll_wrapper, sys_info='water.pdb', Rmin=3.5, Rmax=4.0)
p2 = qmmm.PAP(hl_wrapper, ll_wrapper, sys_info='water.pdb', Rmin=3.5, Rmax=4.5)
p3 = qmmm.PAP(hl_wrapper, ll_wrapper, sys_info='water.pdb', Rmin=4.0, Rmax=5.0)

# partition the QM and buffer zone atoms
p1.find_buffer_zone()
p2.find_buffer_zone()
p3.find_buffer_zone()

# find QM/MM configurations that arise from buffer zone partitioning
p1.find_configurations()
p2.find_configurations()
p3.find_configurations()

#print Rmin, Rmax, number of qm groups, buffer groups, and QM/MM configurations
print(p1.Rmin, p1.Rmax, p1.n_qm_groups, p1.n_buffer_groups, p1.n_configs)
print(p2.Rmin, p2.Rmax, p2.n_qm_groups, p2.n_buffer_groups, p2.n_configs)
print(p3.Rmin, p3.Rmax, p3.n_qm_groups, p3.n_buffer_groups, p3.n_configs)












