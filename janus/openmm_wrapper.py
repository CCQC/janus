import simtk.openmm.app as OM_app
import simtk.openmm as OM
import simtk.unit as OM_unit
from .mm_wrapper import MM_wrapper
import mdtraj as md
import numpy as np
from copy import deepcopy

"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""

class OpenMM_wrapper(MM_wrapper):

    def __init__(self, system):

        # OpenMM_wrapper class inherits from MM_wrapper super class
        super().__init__(system, "OpenMM")

        self._pdb = OpenMM_wrapper.create_pdb(self._system.mm_pdb_file)
        self._positions = self._pdb.getPositions(asNumpy=True)/OM_unit.nanometer
        # self._positions *= MM_wrapper.nm_to_angstrom

        self._primary_subsys_modeller = None
        self._second_subsys_modeller = None

        self._boundary['entire_sys'] = {}
        self._boundary['second_subsys'] = {}
        self._boundary['primary_subsys'] = {}

        # save forcefield object
        self._ff = OM_app.ForceField(self._system.mm_ff, self._system.mm_ff_water)

        # find if there are any bonds that are cut across the QM/MM boundary
        self._boundary_bonds = self.find_boundary_bonds()

        if self._boundary_bonds:
            # Get info for adding link atom to primary subsystem
            if self._system.boundary_treatment == 'link_atom':
                self.prepare_link_atom()
            # Get info for adding link atom according to RC or RCD scheme
            if self._system.boundary_treatment == 'RC' or self._system.boundary_treatment == 'RCD':
                self.prepare_link_atom(RC=True)


    def initialize_system(self):

        self.main_OM_system, self.main_simulation, self.main_info =\
        self.get_info(self._pdb, initialize=True, charges=True, get_coulomb=True)

       # # return a MDtraj trajectory object
       # # I don't know if putting these positions and topology is okay or wait until minimize system- 
       # # should I always minimize or always not minimize?
       # return md.Trajectory(self._positions, self._pdb.topology)


    def take_step(self, force):

        for f, coord in force.items():
            # need to figure out if the first 2 parameters always the same or not
            self.qmmm_force.setParticleParameters(f, f, coord)

        self.qmmm_force.updateParametersInContext(self.main_simulation.context)
        self.main_simulation.step(1)
        self.main_info = self.get_main_info()
        self._positions = self.main_info['positions']
    

    def get_main_info(self):
        
        return OpenMM_wrapper.get_state_info(self.main_simulation)
        

    def find_boundary_bonds(self, qm_atoms=None):
        """
        Identified any covalent bonds that the QM/MM boundary cuts across

        Parameters
        ----------
        qm_atoms: A list of atom indicies corresponding to the atoms in
                  the primary subsystem. Default list is taken from qm_atoms
                  stored in the System object

        Returns
        -------
        A list of tuples corresponding to the qm and mm atoms (as OpenMM atom object) involved in every bond
        that need to be cut.

        Examples
        --------
        bonds = find_boundary_bonds()
        bonds = find_boundary_bonds(qm_atoms=[0,1,2,3])
        """

        if qm_atoms is None:
            qm_atoms = self._system.qm_atoms
        bonds = []
        # determining if there are bonds that need to be cut
        for bond in self._pdb.topology.bonds():
            # find any bonds that involve the qm atoms
            if bond.atom1.index in qm_atoms or bond.atom2.index in qm_atoms:
                # isolate bonds that involve one in the qm atoms and one outside
                if bond.atom1.index not in qm_atoms or bond.atom2.index not in qm_atoms:
                    qm = {}
                    mm = {}
                    if bond.atom1.index in qm_atoms:
                        qm = bond.atom1
                        mm = bond.atom2
                    else:
                        qm = bond.atom2
                        mm = bond.atom1

                    bonds.append((qm, mm))
        return bonds


    def entire_sys_info(self, coulomb=True):
        """
        Gets the information for an entire system and saves the OpenMM system object,
        OpenMM simulation object, and a dictionary of relevant information to self

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        None

        Examples
        --------
        entire_sys_info()
        entire_sys_info(coulomb=False)
        """
        self._entire_sys_system, self._entire_sys_simulation, self._entire_sys =\
        self.get_info(self._pdb, charges=True, get_coulomb=coulomb)

    def second_subsys_info(self, coulomb=True):
        """
        Gets the information for the secondary subsystem (anything to be treated with MM)
        and saves the OpenMM system object, OpenMM simulation object, and a dictionary of
        relevant information to self

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        None

        Examples
        --------
        second_subsys_info()
        second_subsys_info(coulomb=False)
        """
        self._second_subsys_modeller = self.create_modeller(keep_qm=False)
        self._second_subsys_system, self._second_subsys_simulation, self._second_subsys =\
        self.get_info(self._second_subsys_modeller, charges=True, get_coulomb=coulomb)

    def primary_subsys_info(self, link=False, coulomb=True):
        """
        Gets the information for the primary subsystem (anything to be treated with QM)
        and saves the OpenMM system object, OpenMM simulation object, and a dictionary of
        relevant information to self
        Note: This function currently only works systems with one link atom!

        Parameters
        ----------
        link: Whether to include a link atom in the MM computation. 
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        None

        Examples
        --------
        primary_subsys_info(link=True)
        primary_subsys_info(coulomb=False)
        """

        self._primary_subsys_modeller = self.create_modeller(keep_qm=True)

        if link is True and self._boundary_bonds:
            #print('link is true')
            if self._system.boundary_treatment == 'link_atom':
                # this structure only working for adding one link atom for now. 
                # NOT FUNCTIONAL FOR more than 1 link atom!!!!!!!!
                for atom in self.link_atoms:
                    #print('getting link modeller')
                    self._primary_subsys_modeller_link = self.create_link_atom_modeller(self._primary_subsys_modeller, self.link_atoms[atom])
                    #print('getting  modeller info')
                    self._primary_subsys_system, self._primary_subsys_simulation, self._primary_subsys =\
                    self.get_info(self._primary_subsys_modeller_link, get_coulomb=coulomb, set_link_charge=True)
            
        else:
            #print('not doing link')
            self._primary_subsys_system, self._primary_subsys_simulation, self._primary_subsys =\
            self.get_info(self._primary_subsys_modeller, get_coulomb=coulomb)


    def boundary_info(self, coulomb=True):
        """
        Gets the information for the interaction energy between the primary and secondary subsystem (QM-MM)
        and saves the OpenMM system object, OpenMM simulation object, and a dictionary of
        relevant information to self. Current implemented by subtracting the energy of the second subsystem
        and primary subsystem from the energy of the entire system.

        Parameters
        ----------
        coulomb: Whether to include coulombic interactions in MM computation.
                 Default is True.

        Returns
        -------
        None

        Examples
        --------
        boundary_info()
        boundary_info(coulomb=False)
        """

        if not self._boundary['entire_sys']:
            self._boundary['entire_sys'] = self.get_info(self._pdb, return_system=False, return_simulation=False,get_coulomb=coulomb)
        if self._second_subsys_modeller is None:
            self._second_subsys_modeller = self.create_modeller(keep_qm=False)
        if not self._boundary['second_subsys']:
            self._boundary['second_subsys'] = self.get_info(self._second_subsys_modeller, return_system=False, return_simulation=False,get_coulomb=coulomb)
        if self._primary_subsys_modeller is None:
            self._primary_subsys_modeller = self.create_modeller(keep_qm=True)
        if not self._boundary['primary_subsys']:
            self._boundary['primary_subsys'] = self.get_info(self._primary_subsys_modeller, return_system=False, return_simulation=False, get_coulomb=coulomb)

        self._boundary['energy'] = self._boundary['entire_sys']['energy'] \
                                - self._boundary['second_subsys']['energy'] \
                                - self._boundary['primary_subsys']['energy']


    def qm_positions(self):
        """
        Grabs the positions of the atoms in the primary subsystem from self._positions
        and makes a string with the element and xyz geometry coordinates. Adds the link atom positions
        when link atoms are needed.
        Note: In a MD time step, does the position update or not? Need to make sure this updates

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        qm_positions()
        """

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        qm_atoms = self._system.qm_atoms

        for idx in qm_atoms:
            # when I phase function out won't need the for and if statements following this
            # see implementation in AQMMM
            for atom in self._pdb.topology.atoms():
                if atom.index == idx:
                    x, y, z =   self._positions[idx][0]*MM_wrapper.nm_to_angstrom,\
                                self._positions[idx][1]*MM_wrapper.nm_to_angstrom,\
                                self._positions[idx][2]*MM_wrapper.nm_to_angstrom
                    out += line.format(atom.element.symbol, x, y, z)

        # if there are bonds that need to be cut
        if self._boundary_bonds:
            # Need to add if statement for any treatments that don't need link atoms
            #if self._system.boundary_treatment !=
            for atom in self.link_atoms:
                pos = self.link_atoms[atom]['link_positions']*MM_wrapper.nm_to_angstrom
                x, y, z = pos[0], pos[1], pos[2]
                out += line.format(self.link_atoms[atom]['link_atom'], x, y, z)

        self._qm_positions = out


    def get_info(self, pdb, initialize=False, charges=False, return_system=True, return_simulation=True, get_coulomb=True, set_link_charge=False):
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
        else:
            OM_system = self.create_openmm_system(pdb)


        # If in electrostatic embedding scheme need to get a system without coulombic interactions
        if self._system.embedding_method=='Electrostatic' and get_coulomb is False:
            # get the nonbonded force
            force = OM_system.getForce(3)
            for i in range(force.getNumParticles()):
                a = force.getParticleParameters(i)
                Sig, Eps = a[1]/OM_unit.nanometer, a[2]/OM_unit.kilojoule_per_mole
                # set the charge to 0 so the coulomb energy is zero
                force.setParticleParameters(i, charge=0, sigma=Sig, epsilon = Eps)

        # Why is this only Mechanical - need to check!
        if self._system.embedding_method=='Mechanical' and set_link_charge is True:
            # set charge of link atom to be zero, assumes link atom is last
            force = OM_system.getForce(3)
            idx = OM_system.getNumParticles() - 1
            a = force.getParticleParameters(idx)
            Sig, Eps = a[1]/OM_unit.nanometer, a[2]/OM_unit.kilojoule_per_mole
            # set the charge to 0 so the coulomb energy is zero
            force.setParticleParameters(idx, charge=0, sigma=Sig, epsilon = Eps)



        # Create an OpenMM simulation from the openmm system, topology, and positions.
        simulation = self.create_openmm_simulation(OM_system,pdb)

        # Calls openmm wrapper to get information specified
        state = OpenMM_wrapper.get_state_info(simulation,
                                      energy=True,
                                      positions=True,
                                      forces=True)

        state['energy'] = state['potential'] + state['kinetic']

        if charges is True:
            state['charges'] = [OM_system.getForce(3).getParticleParameters(i)[0]/OM_unit.elementary_charge for i in range(OM_system.getNumParticles())]

        if return_system is True and return_simulation is True:
            return OM_system, simulation, state
        elif return_system is True and return_simulation is False:
            return OM_system, state
        elif return_system is False and return_simulation is True:
            return simulation, state
        else:
            return state


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


    def create_openmm_system(self, pdb, initialize=False,
                             nonbond=OM_app.NoCutoff, 
                             nonbond_cutoff=1*OM_unit.nanometer,
                             periodic=False,
                             cnstrnts=OM_app.HBonds,
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
        unmatched = self._ff.getUnmatchedResidues(pdb.topology)
        if unmatched:
            self.create_new_residue_template(pdb.topology)

        if periodic is True:
            openmm_system = self._ff.createSystem(pdb.topology,
                                            constraints=cnstrnts)
        else:
            openmm_system = self._ff.createSystem(pdb.topology,
                                            nonbondedMethod=nonbond,
                                            nonbondedCutoff=nonbond_cutoff,
                                            constraints=cnstrnts,
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


        return openmm_system

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
        template, unmatched_res = self._ff.generateTemplatesForUnmatchedResidues(topology)

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
        if name in self._ff._templates:
            print('override existing modified residues with same name')
            template[i].overrideLevel = self._ff._templates[name].overrideLevel + 1

        # register the new template to the forcefield object
        print('register the new template to the forcefield object')
        self._ff.registerResidueTemplate(template[i])


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
        temp = 300*OM_unit.kelvin
        fric = 1/OM_unit.picosecond
        size = 0.002*OM_unit.picoseconds

        integrator = OM.LangevinIntegrator(temp, fric, size)
        simulation = OM_app.Simulation(pdb.topology, openmm_system, integrator)
        simulation.context.setPositions(pdb.positions)

        return simulation

    def get_state_info(simulation,
                       main = False,
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



    def prepare_link_atom(self, RC=False):
        """
        Saves the qm and mm atom associated with the bond being cut across the QM/MM boundary
        and computes the g scaling factor for the link atom.

        Parameters
        ----------
        RC: a bool specifying whether to find the indices of the atoms bonded to mm atom of 
            bond being cut. Default is false

        Returns
        -------
        None

        Examples
        --------
        prepare_link_atom(RC=True)
        """

        for i, bond in enumerate(self._boundary_bonds):

            qm = bond[0]
            mm = bond[1]

            self.link_atoms[i] = {}

            # saving id because id does not change between sys and modeller
            self.link_atoms[i]['qm_id'] = qm.id

            self.link_atoms[i]['mm_id'] = mm.id
            self.link_atoms[i]['mm_index'] = mm.index

            self.link_atoms[i]['link_atom'] = self._system.link_atom
            g = self._system.compute_scale_factor_g(qm.element.symbol, mm.element.symbol, self._system.link_atom)
            self.link_atoms[i]['g_factor'] = g
            # this is in nm
            self.link_atoms[i]['link_positions'] = self._system.get_link_atom_position(self._positions[qm.index], self._positions[mm.index], g)

            # MAYBE CAN USE MDTRAJ FIND NEIGHBORS FOR THIS!!!!!!!!
            if RC is True:
                bonds = []
                # find index of atoms bonded to mm atom
                for bond in self._pdb.topology.bonds():
                    if bond.atom1.id == mm.id or bond.atom2.id == mm.id:
                        if bond.atom1.id != qm.id and bond.atom2.id != qm.id:
                            if bond.atom1.id != mm.id:
                                bonds.append(bond.atom1.index)
                            elif bond.atom2.id != mm.id:
                                bonds.append(bond.atom2.index)
                self.link_atoms[i]['bonds_to_mm'] = bonds


    def create_link_atom_modeller(self, mod, atom):
        '''
        Creates an OpenMM modeller object that includes any link atoms.
        Note: Currently adds a very specfic link atom, need to expand 

        Parameters
        ----------
        mod: modeller object of primary system without link atom added
        atom: dictionary containing information for the link atom 

        Returns
        -------
        A OpenMM modeller object

        Examples
        --------
        create_link_atom_modeller(mod=modeller, atom=link)
        '''
        # get element object
        link = OM_app.element.Element.getBySymbol(atom['link_atom'])

        # get residue where qm atom is
        for atm in mod.topology.atoms():
            if atm.id == atom['qm_id']:
                qm_res = atm.residue

       # add link atom
       # this is a VERY specific case with H1 - need to determine what type of H in the future
        mod.topology.addAtom(name='H1', element=link, residue=qm_res, id='link')

        # add bond between link atom and qm atom
        for atom1 in mod.topology.atoms():
            for atom2 in mod.topology.atoms():
                if atom1.id == atom['qm_id'] and atom2.id == 'link':
                    mod.topology.addBond(atom2, atom1)

        # add link atom position
        positions = mod.getPositions()/OM_unit.nanometer
        pos = OM.vec3.Vec3(atom['link_positions'][0], atom['link_positions'][1], atom['link_positions'][2])
        positions.append(pos)
        mod.positions = positions*OM_unit.nanometer

        return mod
