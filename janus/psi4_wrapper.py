import psi4

def get_psi4_energy(system):

    psi4.core.set_output_file('output.dat', False)
    mol = psi4.geometry(system.qm_molecule)
    psi4.set_options(system.qm_param)
    energy, wavefunction = psi4.energy(system.qm_method, return_wfn=True)
    system.qm_energy = energy
