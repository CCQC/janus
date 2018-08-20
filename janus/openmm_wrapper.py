import simtk.openmm.app as OM_app
import simtk.openmm as OM
import simtk.unit as OM_unit
from .mm_wrapper import MM_wrapper
import numpy as np
from copy import deepcopy

"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""

class OpenMM_wrapper(MM_wrapper):

    def __init__(self, config):

        # OpenMM_wrapper class inherits from MM_wrapper super class
        super().__init__(config, "OpenMM")

        if 'mm_forcefield' in config:
            self.ff = config['mm_forcefield']
        else:
            self.ff = 'amber99sb.xml' 

        if 'mm_forcefield_water' in config:
            self.ff_water = config['mm_forcefield_water']
        else:
            self.ff_water = 'tip3p.xml'

        #need to tinker with these and figure out if specific to openmm
        if 'mm_nonbond_method' in config:
            self.nonbond_method = config['mm_nonbond_method']
        else:
            self.nonbond_method=OM_app.NoCutoff
            
        if 'mm_nonbond_cutoff' in config:
            self.nonbond_cutoff = config['mm_nonbond_cutoff']
        else:
            self.nonbond_cutoff = 1*OM_unit.nanometer

        if 'mm_constraints' in config:
            self.constraints = config['mm_constraints']
        else:
            self.constraints = OM_app.HBonds

        if 'is_periodic' in config:
            self.is_periodic = config['is_periodic']
        else:
            self.is_periodic is False

        if 'mm_fric_coeff' in config:
            self.fric_coeff = config['mm_fric_coeff']
        else:
            self.fric_coeff = 1/OM_unit.picosecond

        if 'embedding_method' in config:
            self.embedding_method = config['embedding_method']
        else:
            self.embedding_method = 'Mechanical'

        self.temp*OM_unit.kelvin
        self.step_size*OM_unit.picoseconds


        self.pdb = OpenMM_wrapper.create_pdb(self.pdb_file)
        self.positions = self.pdb.getPositions(asNumpy=True)/OM_unit.nanometer
        # self._positions *= MM_wrapper.nm_to_angstrom

        self.primary_subsys_modeller = None
        self.second_subsys_modeller = None

        self.boundary['entire_sys'] = {}
        self.boundary['second_subsys'] = {}
        self.boundary['primary_subsys'] = {}

        # save forcefield object
        self.forcefield = OM_app.ForceField(self.ff, self.ff_water)

    def initialize_system(self):

        self.main_simulation, self.main_info =\
        self.compute_mm(self.pdb, initialize=True, return_simulation=True, charges=True, get_coulomb=True)


    def take_step(self, force):

        for f, coord in force.items():
            # need to figure out if the first 2 parameters always the same or not
            # convert this back to openmm units
            coord *= MM_wrapper.au_bohr_to_kjmol_nm
            self.qmmm_force.setParticleParameters(f, f, coord)

        self.qmmm_force.updateParametersInContext(self.main_simulation.context)
        self.main_simulation.step(1)
        self.main_info = self.get_main_info()
        self.positions = self.main_info['positions']
    

    def get_main_info(self):
        
        return OpenMM_wrapper.get_state_info(self.main_simulation, main_info=True)

    def compute_mm(self, pdb, initialize=False, return_system=False, return_simulation=False, set_link_charge=False):
        """
        Gets information about a set of molecules as defined in the pdb, including energy, positions, forces

        Parameters
        ----------
        pbd : a OpenMM pdb object or OpenMM modeller object that contains the relevant system
        charges : a bool to specify whether to get the charges on the MM molecules. 
                  Default is False.
        return_system: a bool to specify whether to return the OpenMM system object. 
                       Default is True. 
        return_simulation: a bool to specify whether to return the OpenMM simulation object.
                       Default is True. 
        get_coulomb: Whether to include coulomic interactions in the MM computation, default is true
        set_link_charge: a bool to specify whether to set the charge of the link atom to zero.
                         Default is False.

        Returns
        -------
        system, simulation, state
        state: A dictionary with state information
        system: OpenMM system object returned unless return_system=False
        simulation: OpenMM simulation object returned unless return_simulation=False


        Examples
        --------
        system, simulation, state = System.get_info(mm_pdb)
        state = System.get_info(mm_pdb, charges=True, return_simulation=False, return_system=False)
        """

        # Create an OpenMM system from an object's topology

        if initialize is True:
            OM_system = self.create_openmm_system(pdb, initialize=True)
            self.main_charges = [OM_system.getForce(3).getParticleParameters(i)[0]/OM_unit.elementary_charge for i in range(OM_system.getNumParticles())]
        else:
            OM_system = self.create_openmm_system(pdb)


        '''
        NEED OUTER LAYER FOR THIS
        '''

        # Create an OpenMM simulation from the openmm system, topology, and positions.
        simulation = self.create_openmm_simulation(OM_system,pdb)

        # Calls openmm wrapper to get information specified
        state = OpenMM_wrapper.get_state_info(simulation,
                                      energy=True,
                                      positions=True,
                                      forces=True)

        if return_system is True and return_simulation is True:
            return OM_system, simulation, state
        elif return_system is True and return_simulation is False:
            return OM_system, state
        elif return_system is False and return_simulation is True:
            return simulation, state
        else:
            return state


    def create_openmm_system(self, pdb, initialize=False,
                             residue={}):
        """
        Calls OpenMM to create an OpenMM System object give a topology,
        forcefield, and other paramters

        Parameters
        ----------
        pdb: a pdb object 
        nonbond: What kind of long range interactions to include. 
        nonbond_cutoff: The cutoff for including long range interactions. Default is 1 nm.
        periodic: A bool to specify whether the system is periodic. Default is False.
        cnstrnts : contraints on the system. Default is HBonds
        residue: a dictionary with any custom residue templates. Default is empty.

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
        """

        # check to see if there are unmatched residues in pdb, create residue templates if there are
        unmatched = self.forcefield.getUnmatchedResidues(pdb.topology)
        if unmatched:
            self.create_new_residue_template(pdb.topology)

        if self.is_periodic is True:
            openmm_system = self.forcefieldcreateSystem(pdb.topology,
                                            constraints=self.constraints)
        else:
            openmm_system = self.forcefield.createSystem(pdb.topology,
                                            nonbondedMethod=self.nonbond_method,
                                            nonbondedCutoff=self.nonbond_cutoff,
                                            constraints=self.contraints,
                                            residueTemplates=residue,
                                            ignoreExternalBonds=False)



        if initialize is True:
            # this is for the initialization of the entire system
            # define a custom force for adding qmmm gradients
            self.qmmm_force = OM.CustomExternalForce("-x*fx-y*fy-z*fz")
            self.qmmm_force.addPerParticleParameter('fx')
            self.qmmm_force.addPerParticleParameter('fy')
            self.qmmm_force.addPerParticleParameter('fz')
            
            for i in range(openmm_system.getNumParticles()):
                self.qmmm_force.addParticle(i, np.array([0.0, 0.0, 0.0]))
            
            openmm_system.addForce(self.qmmm_force)

        # If in electrostatic embedding scheme need to get a system without coulombic interactions
        if self.embedding_method=='Electrostatic' and include_coulomb is False:
            # get the nonbonded force
            force = OM_system.getForce(3)
            for i in range(force.getNumParticles()):
                a = force.getParticleParameters(i)
                Sig, Eps = a[1]/OM_unit.nanometer, a[2]/OM_unit.kilojoule_per_mole
                # set the charge to 0 so the coulomb energy is zero
                force.setParticleParameters(i, charge=0, sigma=Sig, epsilon = Eps)

        # Why is this only Mechanical - need to check!
        if self.embedding_method=='Mechanical' and set_link_charge is True:
            # set charge of link atom to be zero, assumes link atom is last
            force = OM_system.getForce(3)
            idx = OM_system.getNumParticles() - 1
            a = force.getParticleParameters(idx)
            Sig, Eps = a[1]/OM_unit.nanometer, a[2]/OM_unit.kilojoule_per_mole
            # set the charge to 0 so the coulomb energy is zero
            force.setParticleParameters(idx, charge=0, sigma=Sig, epsilon = Eps)


        return openmm_system

    def create_modeller(self, keep_qm=None):
        """
        Makes a OpenMM modeller object based on given geometry

        Parameters
        ----------
        keep_qm : a bool of whether to keep the qm atoms in the
                  modeller or delete them.
                  The default is to make a modeller without the qm atoms

        Returns
        -------
        A OpenMM modeller object

        Examples
        --------
        modeller = self.make_modeller()
        modeller = self.make_modeller(keep_qm=True)
        """

        modeller = OM_app.Modeller(self._pdb.topology, self._positions)
        if keep_qm is False:
            OpenMM_wrapper.delete_atoms(modeller, self._system.qm_atoms)
        elif keep_qm is True:
            OpenMM_wrapper.keep_atoms(modeller, self._system.qm_atoms)
        return modeller

    def create_new_residue_template(self, topology):

        """
        Create a new OpeMM residue template when there is no matching residue and registers it into self._ff
        forcefield object.
        Note: currently, if there is unmatched name, currently only checks original unmodified
              residue, N-terminus form, and C-terminus form. This may not be robust.

        Parameters
        ----------
        topology: an OpenMM topology object

        Returns
        -------
        None

        Examples
        --------
        create_new_residue_template(topology)
        """
        template, unmatched_res = self.ff.generateTemplatesForUnmatchedResidues(topology)

        # Loop through list of unmatched residues
        print('Loop through list of unmatched residues')
        for i, res in enumerate(unmatched_res):
            res_name = res.name                             # get the name of the original unmodifed residue
            n_res_name = 'N' + res.name                     # get the name of the N-terminus form of original residue
            c_res_name = 'C' + res.name                     # get the name of the C-terminus form of original residue
            name = 'Modified_' + res_name                   # assign new name
            template[i].name = name

        # loop through all atoms in modified template and all atoms in orignal template to assign atom type
        print('loop through all atoms in modified template and all atoms in orignal template to assign atom type')
        for atom in template[i].atoms:
            for atom2 in self._ff._templates[res_name].atoms:
                if atom.name == atom2.name:
                    atom.type = atom2.type
            # the following is for when there is a unmatched name, check the N and C terminus residues
            if atom.type == None:
                print('check n')
                for atom3 in self._ff._templates[n_res_name].atoms:
                    if atom.name == atom3.name:
                        atom.type = atom3.type
            if atom.type == None:
                print('check c')
                for atom4 in self._ff._templates[c_res_name].atoms:
                    if atom.name == atom4.name:
                        atom.type = atom4.type

        # override existing modified residues with same name
        print(name)
        if name in self.ff._templates:
            print('override existing modified residues with same name')
            template[i].overrideLevel = self.ff._templates[name].overrideLevel + 1

        # register the new template to the forcefield object
        print('register the new template to the forcefield object')
        self.ff.registerResidueTemplate(template[i])


    def create_openmm_simulation(self, openmm_system, pdb):
        """
        Creates an OpenMM simulation object given
        an OpenMM system and pdb
        Note: currently options are all default, need way to specify 

        Parameters
        ----------
        openmm_system : OpenMM system object
        pdb: An OpenMM pdb object containing topology and positions

        Returns
        -------
        an OpenMM simulation object

        Examples
        --------
        create_open_simulation(openmm_sys, pdb)
        """

        # The following need to be set as writable options 

        integrator = OM.LangevinIntegrator(self.temp, self.fric_coefficient, self.step_size)
        simulation = OM_app.Simulation(pdb.topology, openmm_system, integrator)
        simulation.context.setPositions(pdb.positions)

        return simulation

    def get_state_info(simulation,
                       main_info=False,
                       energy=True,
                       positions=True,
                       velocity=False,
                       forces=True,
                       parameters=False,
                       param_deriv=False,
                       periodic_box=False,
                       groups_included=-1):
        """
        Gets information like the kinetic and potential energy from an OpenMM state

        Parameters
        ----------
        simulation : an OpenMM simulation object
        energy : a bool for specifying whether to get the energy,
                 returns in hartrees, default is true.
        positions : a bool for specifying whether to get the positions,
                    returns in angstroms, default is true
        velocity : a bool for specifying whether to get the velocities, default is false
        forces : a bool for specifying whether to get the forces acting
                on the system, returns as numpy array in jk/mol/nm, default is true
        parameters : a bool for specifying whether to get the parameters
                    of the state, default is false.
        param_deriv : a bool for specifying whether to get the parameter
                    derivatives of the state, default is false
        periodic_box : a bool for whether to translate the positions so the
                    center of every molecule lies in the same periodic box, default is false
        groups : a set of indices for which force groups to include when computing
                forces and energies. Default is all groups

        Returns
        -------
        A dictionary with information specified by parameters

        Examples
        --------
        get_state_info(sim)
        get_state_info(sim, groups_included=set{0,1,2})
        get_state_info(sim, positions=True, forces=True)
        """

        state = simulation.context.getState(getEnergy=energy,
                                            getPositions=positions,
                                            getVelocities=velocity,
                                            getForces=forces,
                                            getParameters=parameters,
                                            getParameterDerivatives=param_deriv,
                                            enforcePeriodicBox=periodic_box,
                                            groups=groups_included)

        values = {}
        # divide by unit to give value without units, then convert value to atomic units
        if energy is True:
            values['potential'] = state.getPotentialEnergy()/OM_unit.kilojoule_per_mole
            values['potential'] *= MM_wrapper.kjmol_to_au
            values['kinetic'] = state.getKineticEnergy()/OM_unit.kilojoule_per_mole
            values['kinetic'] *= MM_wrapper.kjmol_to_au

        if positions is True:
            values['positions'] = state.getPositions(asNumpy=True)/OM_unit.nanometer
            #values['positions'] *= MM_wrapper.nm_to_angstrom

        if forces is True:
            values['forces'] = state.getForces(asNumpy=True)/(OM_unit.kilojoule_per_mole/OM_unit.nanometer)
            values['gradients'] = (-1) * values['forces'] * MM_wrapper.kj_mol_nm_to_au_bohr   

        if main_info is True:
            # need to check if the topology actually updates 
            values['topology'] = simulation.topology


        values['energy'] = values['potential'] + values['kinetic']

        return values

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

    def create_pdb(mm_pdb_file):
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

        pdb = OM_app.PDBFile(mm_pdb_file)
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
        OM.PDBFile.writeFile(mod.topology, mod.positions, open(filename, 'w'))



