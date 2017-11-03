"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout

def create_openmm_pdb(mm_pdb_file):
    """
    Creates an OpenMM PDB object 

    Parameters
    ----------
    mm_pdb_file: string of pdb file name 

    Returns
    -------
    OpenMM PDB object

    Examples
    --------
    model = create_openmm_pdb('input.pdb')
    """
    
    pdb = PDBFile(mm_pdb_file)
    return pdb

def write_pdb(mod, filename):
    """ 
    Write a pdb file from a OpenMM modeller

    Parameters
    ----------
    mod : OpenMM modeller object
    filename : string of file to write to

    Returns
    -------
    None

    Examples
    --------
    write_pdb(mod, 'input.pdb')
    """
    PDBFile.writeFile(mod.topology, mod.positions, open(filename, 'w'))

def create_openmm_system(topology, forcefield='amber99sb.xml',
                         forcefield_water='tip3p.xml',
                         #nonbond=NoCutoff, nonbond_cutoff=1*nanometer,
                         cnstrnts= HBonds):
    """
    Calls OpenMM to create a OpenMM System object give a topology,
    forcefield, and other paramters

    Parameters
    ----------
    topology : a OpenMM topology 
    forcefield : string of forcefield name to use. Default is amber99sb.xml
    forcefield_water : string of forcefield name to use for water
                       application for amber forcefields that do no
                       define water. Default is tip3p.xml
    cnstrnts : contraints on the system. Default is HBonds

    TODO: need to put nonbond and nonbond_cutoff back but not doing for
          now because need non-periodic system. Other parameters are also needed

          also, expand forcefield to take not openmm built in but customized as well

    Returns
    -------
    An OpenMM system object

    Examples
    --------
    openmm_sys = create_openmm_system(topology)
    openmm_sys = create_openmm_system(pdb.topology, constrnts=None)

    To get OpenMM system information, e.g., Number of particles:
        print(sys.getNumParticles())
    Question - for things like this - do I need a wrapper? since I am technically still
               using openmm -> instead of saving an "OpenMM" object - should I define my
               own objects
    """

    ff = ForceField(forcefield, forcefield_water)

    openmm_system = ff.createSystem(topology,
#                                    nonbondedMethod=nonbond,
#                                    nonbondedCutoff=nonbond_cutoff,
                                    constraints=cnstrnts)
    return openmm_system



def create_openmm_simulation(openmm_system, topology, positions):
    """
    Creates an OpenMM simulation object given
    an OpenMM system, topology, and positions 

    Parameters
    ----------
    openmm_system : OpenMM system object
    topology : an OpenMM topology     
    positions : OpenMM positions

    Returns
    -------
    an OpenMM simulation object
    
    Examples
    --------
    create_open_simulation(openmm_sys, pdb.topology, pdb.positions)
    """
    integrator = LangevinIntegrator(300*kelvin,
                                    1/picosecond, 0.002*picoseconds)

    simulation = Simulation(topology, openmm_system, integrator)
    simulation.context.setPositions(positions)
    return simulation


def get_state_info(simulation, energy=True):
    """
    Gets information like the kinetic 
    and potential energy from an OpenMM state

    Parameters
    ----------
    simulation : an OpenMM simulation object
    energy : a bool for specifying whether to get the energy
    
    TODO: add other state information besides energy 

    Returns
    -------
    OpenMM Quantity objects in kcal/mol 

    Examples
    --------
    get_openmm_energy(sim)
    """
    state = simulation.context.getState(getEnergy=energy)
    potential = state.getPotentialEnergy()
    kinetic = state.getKineticEnergy()
    return potential, kinetic


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


def keep_residues(model, residues):
    """
    Acts on a OpenMM Modeller object to keep the specified
    residue in the MM system and deletes everything else

    Parameters
    ----------
    model : OpenMM Modeller object
    residues : list of residues names to keep in an OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    keep_residue(mod, ['HOH'])
    """
    lis = []
    for res in model.topology.residues():
        if res.name not in residues:
            lis.append(res)
    model.delete(lis)


def keep_atoms(model, atoms):
    """
    Acts on a OpenMM Modeller object to keep the specified
    atoms in the MM system and deletes everything else

    Parameters
    ----------
    model : OpenMM Modeller object
    atoms : list of atoms to keep in an OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    keep_atom(mod, [0,1])
    """
    lis = []
    for atom in model.topology.atoms():
        if atom.index not in atoms:
            lis.append(atom)
    model.delete(lis)


def delete_residues(model, residues):
    """
    Delete specified residues from a OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object
    residues : list of residue IDs (int) to delete from
                  the OpenMM Modeller object

    In future: expand to take residue names, other forms of id

    Returns
    -------
    None

    Examples
    --------
    delete_residues(model, [0, 3, 5])
    """
    lis = []
    for res in model.topology.residues():
        if res.index in residues:
            lis.append(res)
    model.delete(lis)


def delete_atoms(model, atoms):
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
    delete_atoms(model, [0, 3, 5])
    """
    lis = []
    for atom in model.topology.atoms():
        if atom.index in atoms:
            lis.append(atom)
    model.delete(lis)


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
