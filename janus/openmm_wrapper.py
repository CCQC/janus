"""
This module is a wrapper that calls OpenMM
to obtain MM information
"""
import simtk.openmm.app as OM_app
import simtk.openmm as OM 
import simtk.unit as OM_unit
from .mm_wrapper import MM_wrapper

class OpenMM_wrapper(MM_wrapper):

    def __init__(self, system):
        super().__init__(system, "OpenMM")

        self._pdb = OpenMM_wrapper.create_pdb(self._system.mm_pdb_file)
        self._positions = self._pdb.getPositions(asNumpy=True)/OM_unit.nanometer
#        self._positions *= MM_wrapper.nm_to_angstrom

        self._primary_subsys_modeller = None
        self._second_subsys_modeller = None

        self._boundary['entire_sys'] = {}
        self._boundary['second_subsys'] = {}
        self._boundary['primary_subsys'] = {}

        self._ff = OM_app.ForceField(self._system.mm_ff, self._system.mm_ff_water)

        self._boundary_bonds = self.find_boundary_bonds()

        if self._boundary_bonds:
            if self._system.boundary_treatment == 'link_atom':
                self.link_atoms = {}
                self.prepare_link_atom() 
    

    def find_boundary_bonds(self, qm_atoms=None):

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


    def entire_sys_info(self):
        self._entire_sys_system, self._entire_sys_simulation, self._entire_sys = self.get_info(self._pdb)
        
    def second_subsys_info(self):
        self._second_subsys_modeller = self.create_modeller(keep_qm=False)
        self._second_subsys_system, self._second_subsys_simulation, self._second_subsys = self.get_info(self._second_subsys_modeller, charges=True)

    def primary_subsys_info(self, link=False):

        self._primary_subsys_modeller = self.create_modeller(keep_qm=True)

        if link is True: 
            if self._boundary_bonds and self._system.boundary_treatment == 'link_atom':  
                # this structure only working for adding one link atom for now
                for atom in self.link_atoms:
                    self._primary_subsys_modeller_link = self.create_link_atom_modeller(self._primary_subsys_modeller, self.link_atoms[atom])
                    self._primary_subsys_system, self._primary_subsys_simulation, self._primary_subsys = self.get_info(self._primary_subsys_modeller_link)
        else:
            self._primary_subsys_system, self._primary_subsys_simulation, self._primary_subsys = self.get_info(self._primary_subsys_modeller)



    def boundary_info(self):

        if not self._boundary['entire_sys']:
            self._boundary['entire_sys'] = self.get_info(self._pdb, return_system=False, return_simulation=False)
        if self._second_subsys_modeller is None:
            self._second_subsys_modeller = self.create_modeller(keep_qm=False)
        if not self._boundary['second_subsys']: 
            self._boundary['second_subsys'] = self.get_info(self._second_subsys_modeller, return_system=False, return_simulation=False)
        if self._primary_subsys_modeller is None:
            self._primary_subsys_modeller = self.create_modeller(keep_qm=True)
        if not self._boundary['primary_subsys']:
            self._boundary['primary_subsys'] = self.get_info(self._primary_subsys_modeller, return_system=False, return_simulation=False)

        self._boundary['energy'] = self._boundary['entire_sys']['energy'] \
                                - self._boundary['second_subsys']['energy'] \
                                - self._boundary['primary_subsys']['energy']


    def qm_positions(self):

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        qm_atoms = self._system.qm_atoms
    
        for idx in qm_atoms:
            for atom in self._pdb.topology.atoms():
                if atom.index == idx:
                    x, y, z =   self._positions[idx][0]*MM_wrapper.nm_to_angstrom,\
                                self._positions[idx][1]*MM_wrapper.nm_to_angstrom,\
                                self._positions[idx][2]*MM_wrapper.nm_to_angstrom
                    out += line.format(atom.element.symbol, x, y, z)

        # if there are bonds that need to be cut
        if self._boundary_bonds:
            if self._system.boundary_treatment == 'link_atom':
                for atom in self.link_atoms:
                    pos = self.link_atoms[atom]['link_positions']*MM_wrapper.nm_to_angstrom
                    x, y, z = pos[0], pos[1], pos[2]
                    out += line.format(self.link_atoms[atom]['link_atom'], x, y, z)

        self._qm_positions = out

        
    def get_info(self, pdb, charges=False, return_system=True, return_simulation=True):
        """
        Gets information about a system, e.g., energy, positions, forces

        Parameters
        ----------
        pbd : a OpenMM pdb object or OpenMM modeller object
        charges : a bool to specify whether to get the charges on the
                  MM molecules. Default is false

        Returns
        -------
        A OpenMM system object, a dictionary containing information
        for the system


        Examples
        --------
        system, state = System.get_info(mm_pdb)
        system, state = System.get_info(mm_pdb, force='nonbonded', charge=True)
        """

        # Create an OpenMM system from an object's topology
        OM_system = self.create_openmm_system(pdb)

        # Remove Bond, Angle, and Torsion forces to leave only nonbonded forces

        # this part no longer needed 
       # if forces == 'nonbonded':
       #     for i in range(3):
       #         OM_system.removeForce(0)

        if self._system.embedding_method=='Electrostatic':
            # get the nonbonded force
            force = OM_system.getForce(4)
            for i in range(force.getNumParticles()):
                a = force.getParticleParameters(i)
                Sig, Eps = a[1]/OM_unit.nanometer, a[2]/OM_unit.kilojoule_per_mole
                # set the charge to 0 so the coulomb energy is zero
                force.setParticleParameters(i, charge=0, sigma=Sig, epsilon = Eps)

        # Create an OpenMM simulation from the openmm system,
        # topology and positions.
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

        modeller = OM_app.Modeller(self._pdb.topology, self._pdb.positions)
        if keep_qm is False:
            OpenMM_wrapper.delete_atoms(modeller, self._system.qm_atoms)
        elif keep_qm is True:
            OpenMM_wrapper.keep_atoms(modeller, self._system.qm_atoms)
        return modeller




    def create_openmm_system(self, pdb, 
                             nb_forces_only=False,
                             nonbond=OM_app.NoCutoff, nonbond_cutoff=1*OM_unit.nanometer,
                             periodic=False,
                             cnstrnts=OM_app.HBonds,
                             residue={}):
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

        # check to see if there are unmatched residues, create residue templates if there are
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


        return openmm_system

    def create_new_residue_template(self, topology):
        
        template, unmatched_res = self._ff.generateTemplatesForUnmatchedResidues(topology)
 
        # Loop through list of unmatched residues  
        for i, res in enumerate(unmatched_res):
            res_name = res.name                             # get the name of the original unmodifed residue
            n_res_name = 'N' + res.name                     # get the name of the N-terminus form of original residue
            c_res_name = 'C' + res.name                     # get the name of the C-terminus form of original residue
            template[i].name = 'Modified_' + res_name       # assign new name
 
        # loop through all atoms in modified template and all atoms in orignal template to assign atom type
        for atom in template[i].atoms:
            for atom2 in self._ff._templates[res_name].atoms:
                if atom.name == atom2.name:
                    atom.type = atom2.type
            # the following is for when there is a unmatched name, chech the N and C terminus residues
            if atom.type == None:
                for atom3 in self._ff._templates[n_res_name].atoms:
                    if atom.name == atom3.name:
                        atom.type = atom3.type
            if atom.type == None:
                for atom4 in self._ff._templates[c_res_name].atoms:
                    if atom.name == atom4.name:
                        atom.type = atom4.type

        # register the new template to the forcefield object
        self._ff.registerResidueTemplate(template[i])


    def create_openmm_simulation(self, openmm_system, pdb):
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
        sys = self._system
        temp = 300*OM_unit.kelvin
        fric = 1/OM_unit.picosecond
        size = 0.002*OM_unit.picoseconds
        integrator = OM.LangevinIntegrator(temp, fric, size)

        simulation = OM_app.Simulation(pdb.topology, openmm_system, integrator)
        simulation.context.setPositions(pdb.positions)
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
                returns in kj/mol
        positions : a bool for specifying whether to get the positions
                    returns in nm
        velocity : a bool for specifying whether to get the velocities
        forces : a bool for specifying whether to get the forces acting
                on the system
                    returns in jk/mol/nm
        parameters : a bool for specifying whether to get the parameters
                    of the state
        param_deriv : a bool for specifying whether to get the parameter
                    derivatives of the state
        periodic_box : a bool for whether to translate the positions so the
                    center of every molecule lies in the same periodic box
        grouprimary_subsys : a set of indices for which force grouprimary_subsys to include when computing
                forces and energies. Default is all grouprimary_subsys

        TODO: add how to get other state information besides energy

        Returns
        -------
        OpenMM Quantity objects in kcal/mol

        Examples
        --------
        get_state_info(sim)
        get_state_info(sim, grouprimary_subsys_included=set{0,1,2})
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
        # divide by unit to give value without units
        # then convert value to atomic units
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

 

    def prepare_link_atom(self):

        for i, bond in enumerate(self._boundary_bonds):
            qm = bond[0]
            mm = bond[1]

            self.link_atoms[i] = {} 
            #self.link_atoms[i]['qm_atom'] = qm.element.symbol
            #self.link_atoms[i]['qm_index'] = qm.index 
            ## this is in nm
            #self.link_atoms[i]['qm_positions'] = self._positions[qm.index] 

            # saving id because id does not change between sys and modeller
            self.link_atoms[i]['qm_id'] = qm.id

            #self.link_atoms[i]['mm_atom'] = mm.element.symbol
            #self.link_atoms[i]['mm_index'] = mm.index 
            ## this is in nm
            #self.link_atoms[i]['mm_positions'] = self._positions[mm.index] 
            self.link_atoms[i]['mm_id'] = mm.id
            
            self.link_atoms[i]['link_atom'] = self._system.link_atom
            g = self._system.compute_scale_factor_g(qm.element.symbol, mm.element.symbol, self._system.link_atom)
            self.link_atoms[i]['g_factor'] = g 
            # this is in nm
            self.link_atoms[i]['link_positions'] = self._system.get_link_atom_position(self._positions[qm.index], self._positions[mm.index], g) 


    def create_link_atom_modeller(self, mod, atom):
        '''
        mod is modeller of primary system without link atom added
        atom is link_atom dictionary 
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
