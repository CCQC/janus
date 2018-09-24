import simtk.openmm.app as OM_app
import simtk.openmm as OM
import simtk.unit as OM_unit
from .mm_wrapper import MM_wrapper
import numpy as np
from copy import deepcopy


class OpenMM_wrapper(MM_wrapper):
    """
    A wrapper class that calls OpenMM
    to obtain molecular mechanics information and 
    take steps in a molecular dynamics simulation. 
    Class inherits from MM_wrapper.
    """

    def __init__(self, param):
        """
        Initializes an OpenMM wrapper class with a set of parameters for 
        running OpenMM. Creates an OpenMM pdb object using the given pdb_file
        and forcefield object using the given forcefield

        Parameters
        ----------
        param : dict containing parameters for molecular mechanics computations.
                Individual parameters include:
                - mm_pdb_file: a pdb file that contains the system of interest
                - mm_forcefield: the name of the forcefield to use, default is amber99sb.xml
                - mm_water_forcefield: the name of the forcefield to use for water, default is tip3p.xml
                - is_periodic: whether to treat the system periodically, default is False
                - step_size: step size to integrate system in picoseconds, 
                             default is 0.002*OM_unit.picoseconds
                - integrator" : which integrator to use for simulation, default is Langevin
                - fric_coeff" : friction coefficient to couple the system to heat bath in inverse
                                picoseconds, default is 1/OM_unit.picosecond
                - temp: the temperature at which the simulation runs in kelvin, default is 300*OM_unit.kelvin
                - nonbondedMethod: method for nonbonded interactions, default is OM_app.NoCutoff
                - nonbondedCutoff: cutoff distance for nonbonded interactions in nanometers,
                                   default is "1*OM_unit.nanometer",
                - constraints: which bonds and angles implemented with constraints,
                               default is OM_app.HBonds
                - rigid_water: whether water is treated as rigid, default is True
                - removeCMMotion: whether to include a CMMotionRemover, default is True
                - ignoreExternalBonds: whether to ignore external bonds when matching residues to templates,
                                       default is True
                - flexibleConstraints: whether to add parameters for constrained parameters,
                                       default is False
                - hydrogenMass: the mass to use for hydrogen atoms bonded to heavy atoms,
                                default is False
                - residueTemplates: allows user to specify a template for a residue,
                                    default is empty dict {}
                - switchDistance: the distance to turn on potential energy switching function for 
                                  Lennard-Jones interactions. Default is None

                For more information about these pararmeters and 
                other possible parameter values consult docs.openmm.org
                

        Returns
        -------
        An OpenMM_wrapper object

        Examples
        --------
        mm_wrapper = OpenMM_wrapper(param)
        """

        super().__init__(param, "OpenMM")

        self.ff = param['mm_forcefield']
        self.ff_water = param['mm_water_forcefield']
        self.is_periodic = param['is_periodic']

        self.temp = eval(param['temp'])
        self.step_size = eval(param['step_size'])
        self.fric_coeff = eval(param['fric_coeff'])

        # instantiate OpenMM pdb object
        self.pdb = OpenMM_wrapper.create_pdb(self.pdb_file)
        self.positions = None

        # instantiate OpenMM forcefield object
        self.forcefield = OM_app.ForceField(self.ff, self.ff_water)

    def initialize(self, embedding_method):
        """
        Calls compute_mm to get information for the system
        of interest in its initial state and saves the simulation 
        object and information dictionary returned by compute_mm
    
        Parameters
        ----------
        embedding_method: what embedding method to use for initialization.
                          If 'Mechanical', all forces are included
                          If 'Electrostatic', all coulomic forces are excluded

        Returns
        -------
        None

        Examples
        --------
        initialize('Mechanical')
        initialize('Electrostatic')
        """

        # should I minimize energy here? If so, need to return new positions
        if embedding_method == 'Mechanical':
            self.main_simulation, self.main_info =\
            self.compute_mm(self.pdb.topology, self.pdb.positions, initialize=True, return_simulation=True)

        elif embedding_method == 'Electrostatic':
            self.main_simulation, self.main_info =\
            self.compute_mm(self.pdb.topology, self.pdb.positions, include_coulomb=None, initialize=True, return_simulation=True)
        else:
            print('only mechanical and electrostatic embedding schemes implemented at this time')

    def take_step(self, force):
        """
        Updates the system with forces from qmmm 
        and takes a simulation step

        Parameters
        ----------
        forces: dict of forces(particle index: forces) in au/bohr to 
                be updated in custom qmmm force and fed into simulation

        Returns
        -------
        None

        Examples
        --------
        take_step(forces)
        take_step({0: [0.0, 0.0, 0.0]})
        """

        for f, coord in force.items():
            coord *= MM_wrapper.au_bohr_to_kjmol_nm             # convert this back to openmm units
            self.qmmm_force.setParticleParameters(f, f, coord)  # need to figure out if the first 2 parameters always the same or not

        self.qmmm_force.updateParametersInContext(self.main_simulation.context)  # update forces with qmmm force
        self.main_simulation.step(1)                                             # take a step
        self.main_info = self.get_main_info()                                    # get the energy and gradients after step
        self.positions = self.main_info['positions']                             # get positions after step
    

    def get_main_info(self):
        """
        Gets the information for the system of interest
        
        Parameters
        ----------
        None

        Returns
        -------
        dictionary with information like energy and gradient 
        for the system of interest
    
        Examples
        --------
        info = get_main_info()
        """
        
        return OpenMM_wrapper.get_state_info(self.main_simulation, main_info=True)

    def compute_mm(self, topology, positions, include_coulomb='all', initialize=False, return_system=False, return_simulation=False, link_atoms=None, minimize=False):
        """
        Gets information about a set of molecules as defined in the pdb, including energy, positions, forces

        Parameters
        ----------
        topology: an OpenMM topology object
        positions: an OpenMM Vec3 vector containing 
                   the positions of the system in nm
        include_coulomb: whether to include coulombic interactions. 
                         'all' (default) includes coulombic forces for all particles,
                         'no_link' excludes coulombic forces for link atoms,
                         'only' excludes all other forces for all atoms,
                         None excludes coulombic forces for all particles.
        initialize: bool to specifying whether the main system is being initialized.
        return_system: a bool to specify whether to return the OpenMM system object. 
                       Default is True. 
        return_simulation: a bool to specify whether to return the OpenMM simulation object.
                       Default is True. 
        link_atoms: if included as a list with include_coulomb='no_link', specifies which 
                    atoms to remove coulombic forces from. Default is None.

        Returns
        -------
        system, simulation, state 
        state: A dictionary with state information
        system: OpenMM system object returned unless return_system=False
        simulation: OpenMM simulation object returned unless return_simulation=False


        Examples
        --------
        system, simulation, state = compute_mm(top, pos)
        state = compute_mm(top, pos, return_simulation=False, return_system=False)
        """

        # Create an OpenMM system from an object's topology
        if initialize is True:
            OM_system = self.create_openmm_system(topology, include_coulomb, link_atoms,initialize=True)
        else:
            OM_system = self.create_openmm_system(topology, include_coulomb, link_atoms)

        # Create an OpenMM simulation from the openmm system, topology, and positions.
        simulation = self.create_openmm_simulation(OM_system, topology, positions)

        if minimize is True:
            simulation.minimizeEnergy()

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


    def create_openmm_system(self, topology, include_coulomb='all', link_atoms=None, initialize=False):
        """
        Calls OpenMM to create an OpenMM System object give a topology,
        forcefield, and other paramters as given in input

        Parameters
        ----------
        topology: an OpenMM topology object
        include_coulomb: whether to include coulombic interactions. 
                         'all' (default) includes coulombic forces for all particles,
                         'no_link' excludes coulombic forces for link atoms,
                         'only' excludes all other forces for all atoms,
                         None excludes coulombic forces for all particles.
        link_atoms: if included as a list with include_coulomb='no_link', specifies which 
                    atoms to remove coulombic forces from. Default is None.
        initialize: bool to specifying whether the main system is being initialized.
                    If True, this will create a custom force to store qmmm forces. 
                    default is false

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
        openmm_sys = create_openmm_system(top, include_coulomb='no_link', link_atoms=[0,1,2])
        openmm_sys = create_openmm_system(top, include_coulomb='only')
        openmm_sys = create_openmm_system(top, initialize=True)
        """

        # check to see if there are unmatched residues in pdb, create residue templates if there are
        unmatched = self.forcefield.getUnmatchedResidues(topology)
        if unmatched:
            self.create_new_residue_template(topology)

        if self.is_periodic is True:
            print("periodic")
            openmm_system = self.forcefield.createSystem(topology,
                                            constraints=self.constraints)
        else:
            openmm_system = self.forcefield.createSystem(topology,
                                            nonbondedMethod=eval(self.param['nonbondedMethod']),
                                            nonbondedCutoff=eval(self.param['nonbondedCutoff']),
                                            constraints=eval(self.param['constraints']),
                                            residueTemplates=self.param['residueTemplates'],
                                            hydrogenMass=eval(self.param['hydrogenMass']),
                                            switchDistance=eval(self.param['switchDistance']),
                                            rigid_water=self.param['rigid_water'],
                                            removeCMMotion=self.param['removeCMMotion'],
                                            flexibleConstraints=self.param['flexibleConstraints'],
                                            ignoreExternalBonds=self.param['ignoreExternalBonds'])


        if initialize is True:                                             # this is for the initialization of the entire system
            self.qmmm_force = OM.CustomExternalForce("-x*fx-y*fy-z*fz")    # define a custom force for adding qmmm gradients
            self.qmmm_force.addPerParticleParameter('fx')
            self.qmmm_force.addPerParticleParameter('fy')
            self.qmmm_force.addPerParticleParameter('fz')
            
            for i in range(openmm_system.getNumParticles()):
                self.qmmm_force.addParticle(i, np.array([0.0, 0.0, 0.0]))
            
            openmm_system.addForce(self.qmmm_force)

            self.main_charges = [openmm_system.getForce(3).getParticleParameters(i)[0]/OM_unit.elementary_charge for i in range(openmm_system.getNumParticles())]

        # If in electrostatic embedding scheme need to get a system without coulombic interactions
        if include_coulomb is None:
            # get the nonbonded force
            self.set_charge_zero(openmm_system)

        if (include_coulomb == 'no_link' and link_atoms is not None):
            self.set_charge_zero(openmm_system, link_atoms)

        if include_coulomb == 'only':
        # Remove Bond, Angle, and Torsion forces to leave only nonbonded forces
            for i in range(openmm_system.getNumForces()):             
                if type(openmm_system.getForce(0)) is not OM.NonbondedForce:     
                    openmm_system.removeForce(0)                              
            self.set_LJ_zero(openmm_system)


        return openmm_system

    def set_charge_zero(self, OM_system, link_atoms=None):
        """
        Removes the coulombic forces by setting charges of 
        specified atoms to zero

        Parameters
        ----------
        OM_system : an OpenMM system object
        link_atoms : list of link_atoms to set the charge to zero,
                     if link_atoms is None (default), the charge of 
                     all particles in the system will be set to zero
    
        Returns
        -------
        None

        Examples
        --------
        set_charge_zero(system)
        set_charge_zero(system, link_atoms=[0,1,2])
        """

        for force in OM_system.getForces():
            if type(force) is OM.NonbondedForce:
                if link_atoms:
                    # set the charge of link atoms to 0 so the coulomb energy is zero
                    for i in link_atoms:
                        a = force.getParticleParameters(i)
                        force.setParticleParameters(i, charge=0.0, sigma=a[1], epsilon=a[2])
                else:
                    # set the charge of all particles to 0 so the coulomb energy is zero
                    for i in range(force.getNumParticles()):
                        a = force.getParticleParameters(i)
                        force.setParticleParameters(i, charge=0.0, sigma=a[1], epsilon=a[2])

    def set_LJ_zero(self, OM_system):
        """
        Removes the Lennard-Jones (van der Waals) force from the system

        Parameters
        ----------
        OM_system : an OpenMM system object
    
        Returns
        -------
        None

        Examples
        --------
        set_LJ_zero(system)
        """

        for force in OM_system.getForces():
            if type(force) is OM.NonbondedForce:
                for i in range(force.getNumParticles()):
                    a = force.getParticleParameters(i)
                    force.setParticleParameters(i, charge=a[0], sigma=0.0, epsilon=0.0)
        

    def create_new_residue_template(self, topology):
        """
        Create a new OpeMM residue template when there is no matching residue 
        and registers it into self.forcefield forcefield object.
        Note: currently, if there is unmatched name, currently only checks original 
              unmodified residue, N-terminus form, and C-terminus form. 
              This may not be robust.

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
        template, unmatched_res = self.forcefield.generateTemplatesForUnmatchedResidues(topology)

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
                for atom2 in self.forcefield._templates[res_name].atoms:
                    if atom.name == atom2.name:
                        atom.type = atom2.type
                # the following is for when there is a unmatched name, check the N and C terminus residues
                if atom.type == None:
                    print('check n')
                    for atom3 in self.forcefield._templates[n_res_name].atoms:
                        if atom.name == atom3.name:
                            atom.type = atom3.type
                if atom.type == None:
                    print('check c')
                    for atom4 in self.forcefield._templates[c_res_name].atoms:
                        if atom.name == atom4.name:
                            atom.type = atom4.type

            # override existing modified residues with same name
            print(name)
            if name in self.forcefield._templates:
                print('override existing modified residues with same name')
                template[i].overrideLevel = self.forcefield._templates[name].overrideLevel + 1

            # register the new template to the forcefield object
            print('register the new template to the forcefield object')
            self.forcefield.registerResidueTemplate(template[i])


    def create_openmm_simulation(self, openmm_system, topology, positions):
        """
        Creates an OpenMM simulation object given
        an OpenMM system, topology, and positions

        Parameters
        ----------
        openmm_system: OpenMM system object
        topology: an OpenMM topology object
        positions: an OpenMM Vec3 vector containing 
                   the positions of the system in nm

        Returns
        -------
        an OpenMM simulation object

        Examples
        --------
        create_open_simulation(openmm_sys, top, pos) 
        create_open_simulation(openmm_sys, pdb.topology. pdb.positions)
        """

        if self.param['integrator'] == 'Langevin':
            integrator = OM.LangevinIntegrator(self.temp, self.fric_coeff, self.step_size)
        else:
            print('only Langevin integrator supported currently')

        integrator.setRandomNumberSeed(1)

        simulation = OM_app.Simulation(topology, openmm_system, integrator)
        simulation.context.setPositions(positions)

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
        Gets information like the kinetic and potential energy,
        positions, forces, and topology from an OpenMM state.
        Some of these may need to be made accessible to user.

        Parameters
        ----------
        simulation : an OpenMM simulation object
        main_info : a bool for specifying whether to return the topology of the system
        energy : a bool for specifying whether to get the energy,
                 returns in hartrees(a.u.), default is true.
        positions : a bool for specifying whether to get the positions,
                    returns in nanometers, default is true
        velocity : a bool for specifying whether to get the velocities, default is false
        forces : a bool for specifying whether to get the forces acting
                on the system, returns as numpy array in jk/mol/nm, as well as the gradients,
                in au/bohr, default is true
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
        A dictionary with information specified by parameters.
        Keys include 'energy', 'potential', 'kinetic', 'forces',
        'gradients', 'topology'

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
            values['energy'] = values['potential'] + values['kinetic']

        if positions is True:
            values['positions'] = state.getPositions(asNumpy=True)/OM_unit.nanometer

        if forces is True:
            values['forces'] = state.getForces(asNumpy=True)/(OM_unit.kilojoule_per_mole/OM_unit.nanometer)
            values['gradients'] = (-1) * values['forces'] * MM_wrapper.kjmol_nm_to_au_bohr   

        if main_info is True:
            # need to check if the topology actually updates 
            values['topology'] = simulation.topology

        return values

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


    def create_modeller(self, qm_atoms, keep_qm=None):
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

        modeller = OM_app.Modeller(self.pdb.topology, self.pdb.getPositions())
        if keep_qm is False:
            OpenMM_wrapper.delete_atoms(modeller, qm_atoms)
        elif keep_qm is True:
            OpenMM_wrapper.keep_atoms(modeller, qm_atoms)
        return modeller

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

    def get_main_charges(self):
        """
        Gets the MM point charges for the system of interest
        
        Parameters
        ----------
        None

        Returns
        -------
        list of charges
    
        Examples
        --------
        charges = get_main_charges()
        """

        return self.main_charges

