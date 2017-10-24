from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""


def create_openmm_system(pdb_file, forcefield='amber99sb.xml', forcefield_water='tip3p.xml', nonbond=PME, nonbond_cut=1*nanometer, cnstrin=HBonds):
    pdb = PDBFile(pdb_file)
    ff = ForceField(forcefield, forcefield_water)
    system = ff.createSystem(pdb.topology, nonbondedMethod=nonbond, nonbondedCutoff=nonbond_cut, constraints=cnstrin)
    return system, pdb


def create_openmm_simulation(system, pdb, temp=300*kelvin):
    """
    Calls OpenMM with Amber input files
    prmptop and inpcrd to create a MM system
Parameters
    ----------

    Returns
    -------
    None

    Examples
    --------
    get_openmm_amber('input.prmtop', 'input.inpcrd', system)
    """

    #print(system.getNumParticles())
    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    return simulation


def get_openmm_energy(simulation):

    sim = simulation.context.getState(getEnergy=True)
    return sim.getPotentialEnergy() + sim.getKineticEnergy()


def create_openmm_modeller(pdb):
    return Modeller(pdb.topology, pdb.positions)


def keep_residue(model, residue_name):
    """
    --------
    keep_residue(mod, 'HOH')
    """
    qm_list = []
    for res in model.topology.residues():
        if res.name != residue_name:
            qm_list.append(res)
    model.delete(qm_list)


def delete_qm_residues(model, qm_residues):
    qm_list = []
    for res in model.topology.residues():
        if res.index in qm_residues:
            qm_list.append(res)
    model.delete(qm_list)


def delete_qm_atoms(model, qm_atoms):
    qm_list = []
    for atom in test.topology.atoms():
        if atom.index in qm_atoms:
            qm_list.append(atom)
    model.delete(qm_list)


def delete_water(model):
    model.deleteWater()

#sys, pdb = create_openmm_system('input.pdb')
#sim = create_openmm_simulation(sys, pdb)
#print(get_openmm_energy(sim))
