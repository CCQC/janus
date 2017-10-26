"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


def create_openmm_system(pdb_file, forcefield='amber99sb.xml', forcefield_water='tip3p.xml', nonbond=PME, nonbond_cut=1*nanometer, cnstrin=HBonds):
    """
    Creates an OpeMM system. More information about the input parameters
    can be found in the OpenMM documentation for systems

    Parameters
    ----------
    pdb_file : an input .pdb file
    forcefield : forcefield information in .xml format. Default is 'amber99sb.xml'.
    forcefield_water : forcefield used for water molecules in .xml format. Default is 'tip3p.xml'
    nonbond : NonbondedMethod to compute cutoffs for intermolecular interactions. Default is PME.
    nonbond_cut : Nonbonded interaction cutoff. Default is 1*nanometer.
    cnstrin : OpenMM system constraints. Default is HBonds.

    Later:
    expand forcefield to take not openmm built in but customized as well

    Returns
    -------
    OpenMM System object, pdb object created by PDBFile

    Examples
    --------
    sys, pdb = create_openmm_system('input.pdb')
    sys, pdb = create_openmm_system('input.pdb', nonbond=NoCutoff)

    To get OpenMM system information, e.g., Number of particles:
        print(sys.getNumParticles())
    """

    pdb = PDBFile(pdb_file)
    ff = ForceField(forcefield, forcefield_water)
    system = ff.createSystem(pdb.topology, nonbondedMethod=nonbond, nonbondedCutoff=nonbond_cut, constraints=cnstrin)
    return system, pdb


def create_openmm_simulation(mm_system, pdb, temp=300*kelvin):
    """
    Creates a OpenMM simulation object

    Parameters
    ----------
    system : OpenMM System object
    pdb : OpenMM PDBFile object
    temp : Simulation temperature. Default is 300*kelvin

    Returns
    -------
    OpenMM Simulation object

    Examples
    --------
    sim = create_open_simulation(sys, pdb)
    sim = create_open_simulation(sys, pdb, temp=0*kelvin)
    """

    integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    simulation = Simulation(pdb.topology, mm_system, integrator)
    simulation.context.setPositions(pdb.positions)
    return simulation


def get_openmm_energy(simulation):
    """
    Gets the total energy of a OpenMM state

    Parameters
    ----------
    simulation : OpenMM Simulation object

    Returns
    -------
    Total energy = Potential + Kinetic as a OpenMM Quantity object in kcal/mol

    Examples
    --------
    energy = get_openmm_energy(sim)
    To get the value:
    energy._value

    ***need way to specify the unit
    """

    sim = simulation.context.getState(getEnergy=True)
    return sim.getPotentialEnergy() + sim.getKineticEnergy()


def create_openmm_modeller(pdb):
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

    return Modeller(pdb.topology, pdb.positions)


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
    qm_residues : list of qm_residue IDs (int) to delete from the OpenMM Modeller object

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
    qm_residues : list of qm_residue IDs (int) to delete an OpenMM Modeller object 

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
