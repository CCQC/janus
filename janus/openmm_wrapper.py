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
    Write a pdb file from an OpenMM modeller

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
                         # nonbond=NoCutoff, nonbond_cutoff=1*nanometer,
                         cnstrnts=HBonds):
    """
    Calls OpenMM to create an OpenMM System object give a topology,
    forcefield, and other paramters

    Parameters
    ----------
    topology : an OpenMM topology
    forcefield : string of forcefield name to use. Default is amber99sb.xml
    forcefield_water : string of forcefield name to use for water
                       application for amber forcefields that do no
                       define water. Default is tip3p.xml
    cnstrnts : contraints on the system. Default is HBonds

    TODO: need to put nonbond and nonbond_cutoff back but not doing for now
          because need non-periodic system. Other parameters are also needed

          also, expand forcefield to take not openmm built in
            but customized as well

    Returns
    -------
    An OpenMM system object

    Examples
    --------
    openmm_sys = create_openmm_system(topology)
    openmm_sys = create_openmm_system(pdb.topology, constrnts=None)

    To get OpenMM system information, e.g., Number of particles:
        print(sys.getNumParticles())
    Question - for things like this - do I need a wrapper?
                since I am technically still
               using openmm -> instead of saving an "OpenMM" object -
                should I define my own objects
    """

    ff = ForceField(forcefield, forcefield_water)

    openmm_system = ff.createSystem(topology,
                                    # nonbondedMethod=nonbond,
                                    # nonbondedCutoff=nonbond_cutoff,
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


def get_state_info(simulation,
                   energy=True,
                   positions=False,
                   velocity=False,
                   forces=False,
                   parameters=False,
                   param_deriv=False,
                   periodic_box=False,
                   groups_included=-1):
    """
    Gets information like the kinetic
    and potential energy from an OpenMM state

    Parameters
    ----------
    simulation : an OpenMM simulation object
    energy : a bool for specifying whether to get the energy
    positions : a bool for specifying whether to get the positions
    velocity : a bool for specifying whether to get the velocities
    forces : a bool for specifying whether to get the forces acting
             on the system
    parameters : a bool for specifying whether to get the parameters
                 of the state
    param_deriv : a bool for specifying whether to get the parameter
                  derivatives of the state
    periodic_box : a bool for whether to translate the positions so the
                   center of every molecule lies in the same periodic box
    groups : a set of indices for which force groups to include when computing
             forces and energies. Default is all groups

    TODO: add how to get other state information besides energy

    Returns
    -------
    OpenMM Quantity objects in kcal/mol

    Examples
    --------
    get_openmm_energy(sim)
    get_openmm_energy(sim, groups=set{0,1,2})
    """
    state = simulation.context.getState(getEnergy=energy,
                                        getPositions=positions,
                                        getVelocities=velocity,
                                        getForces=forces,
                                        getParameters=parameters,
                                        getParameterDerivatives=param_deriv,
                                        enforcePeriodicBox=periodic_box,
                                        groups=groups_included)

    if energy is True:
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
    Acts on an OpenMM Modeller object to keep the specified
    residue in the MM system and deletes everything else

    Parameters
    ----------
    model : OpenMM Modeller object
    residues : list of residues names (str) or IDs (int)
               to keep in an OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    keep_residue(mod, ['HOH'])
    keep_residue(mod, [0, 1])
    """
    lis = []
    for res in model.topology.residues():
        if type(residues[0]) is int:
            if res.index not in residues:
                lis.append(res)
        elif type(residues[0]) is str:
            if res.name not in residues:
                lis.append(res)
    model.delete(lis)


def keep_atoms(model, atoms):
    """
    Acts on an OpenMM Modeller object to keep the specified
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
    keep_atom(mod, ['O', 'H'])
    """
    lis = []

    for atom in model.topology.atoms():
        if type(atoms[0]) is int:
            if atom.index not in atoms:
                lis.append(atom)
        elif type(atoms[0]) is str:
            if atom.name not in atoms:
                lis.append(atom)
    model.delete(lis)


def delete_residues(model, residues):
    """
    Delete specified residues from an OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object
    residues : list of residue IDs (int) or residue names (str)
               to delete from the OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    delete_residues(model, [0, 3, 5])
    delete_residues(model, ['HOH'])
    """
    lis = []
    for res in model.topology.residues():
        if type(residues[0]) is int:
            if res.index in residues:
                lis.append(res)
        elif type(residues[0]) is str:
            if res.name in residues:
                lis.append(res)
    model.delete(lis)


def delete_atoms(model, atoms):
    """
    Delete specified atoms from an OpenMM Modeller object

    Parameters
    ----------
    model : OpenMM Modeller object
    atoms : list of atom IDs (int) or atom names (str) to delete
                  an OpenMM Modeller object

    Returns
    -------
    None

    Examples
    --------
    delete_atoms(model, [0, 3, 5])
    delete_atoms(model, ['Cl'])
    """
    lis = []
    for atom in model.topology.atoms():
        if type(atoms[0]) is int:
            if atom.index in atoms:
                lis.append(atom)
        elif type(atoms[0]) is str:
            if atom.name in atoms:
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


