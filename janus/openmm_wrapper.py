"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


def create_openmm_system(system):
    """
    Calls OpenMM to create a OpenMM System object, saves
    OpenMM System object as system.mm_system and OpenMM PDB object as
    system.mm_pdb

    Parameters
    ----------
    system : a Janus system object

    Later:
    expand forcefield to take not openmm built in but customized as well

    Returns
    -------
    None

    Examples
    --------
    create_openmm_system(system)

    To get OpenMM system information, e.g., Number of particles:
        print(sys.getNumParticles())
    """

    pdb = PDBFile(system.mm_pdb_file)
    ff = ForceField(system.mm_forcefield, system.mm_forcefield_water)
    openmm_system = ff.createSystem(pdb.topology,
                                    nonbondedMethod=system.mm_nonbond_method,
                                    nonbondedCutoff=system.mm_nonbond_cutoff,
                                    constraints=system.mm_constraints)
    system.mm_system = openmm_system
    system.mm_pdb = pdb


def create_openmm_simulation(system):
    """
    Creates an OpenMM simulation object and saves
    it to the Janus system object as system.qm_sim

    Parameters
    ----------
    system : Janus system object

    Returns
    -------
    None

    Examples
    --------
    create_open_simulation(sys)
    """

    integrator = LangevinIntegrator(system.mm_temp,
                                    1/picosecond, 0.002*picoseconds)
    simulation = Simulation(system.mm_pdb.topology,
                            system.mm_system, integrator)
    simulation.context.setPositions(system.mm_pdb.positions)
    system.mm_sim = simulation


def get_openmm_energy(system):
    """
    Gets the potential and kinetic energy of a OpenMM state
    and save it as a OpenMM Quantity object in kcal/mol to
    the system as system.mm_tot_energy.

    Parameters
    ----------
    system : Janus system object

    Returns
    -------
    None

    Examples
    --------
    get_openmm_energy(sys)
    To get the value:
    energy._value
    TODO: put this in system class

    ***need way to specify the unit
    """

    state = system.mm_sim.context.getState(getEnergy=True)
    system.mm_potential_e = state.getPotentialEnergy()
    system.mm_kinetic_e = state.getKineticEnergy()


def create_openmm_modeller(system):
    """
    Creates an OpenMM Modeller object for changing the MM system

    Parameters
    ----------
    pdb: OpenMM PDBFile object

    Returns
    -------
    OpenMM Modeller object

    Examples
    --------
    model = create_openmm_modeller(pdb)
    """

    return Modeller(system.mm_pdb.topology, system.mm_pdb.positions)


def keep_residue(model, residue_name):
    """
    Acts on a OpenMM Modeller object to keep the specified
    residue in the MM system and deletes everything else

    Parameters
    ----------
    model : OpenMM Modeller object
    residue_name : str of residue to keep in an OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    keep_residue(mod, 'HOH')
    """
    qm_list = []
    for res in model.topology.residues():
        if res.name != residue_name:
            qm_list.append(res)
    model.delete(qm_list)


def delete_qm_residues(model, qm_residues):
    """
    Delete specified residues from a OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object
    qm_residues : list of qm_residue IDs (int) to delete from
                  the OpenMM Modeller object

    In future: expand to take residue names, other forms of id

    Returns
    -------
    None

    Examples
    --------
    delete_qm_residues(model, [0, 3, 5])
    """
    qm_list = []
    for res in model.topology.residues():
        if res.index in qm_residues:
            qm_list.append(res)
    model.delete(qm_list)


def delete_qm_atoms(model, qm_atoms):
    """
    Delete specified atoms from a OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object
    qm_residues : list of qm_residue IDs (int) to delete
                  an OpenMM Modeller object

    In future: expand to take atom names, other forms of id

    Returns
    -------
    None

    Examples
    --------
    delete_qm_atoms(model, [0, 3, 5])
    """
    qm_list = []
    for atom in model.topology.atoms():
        if atom.index in qm_atoms:
            qm_list.append(atom)
    model.delete(qm_list)


def delete_water(model):
    """
    Delete all waters from an OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    delete_water(model)
    """
    model.deleteWater()


def openmm_units():
    pass
# PDBFile.writeFile(mod.topology, mod.positions, open('output.dat', 'w'))
