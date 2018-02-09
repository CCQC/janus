from . import psi4_wrapper as pw
from . import openmm_wrapper as ow


class System:
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    """
    kjmol_to_au = 1/2625.5
    nm_to_angstrom = 10.0

    def __init__(self, qm_param=None, qm_method='scf', qm_molecule=None,
                 qm_charge_method='MULLIKEN_CHARGES',
                 qm_atoms=None, qm_residues=None, mm_pdb_file=None,
                 mm_forcefield='amber99sb.xml',
                 mm_forcefield_water='tip3p.xml',
                 mm_nonbond_method=None, mm_nonbond_cutoff=None,
                 mm_constraints=None, mm_temp=None,
                 embedding_method='Mechanical'):
        """
        Initializes system with None values for all parameters

        Parameters
        ----------
        qm_param: dict of qm parameters
        qm_method: str of desired qm method
        qm_molecule: str of geometry

        I think I should just leave the parameters blank and
        define everything as None..

        Returns
        -------
        A built system

        Examples
        --------
        qm_parameters = System.qm_param()
        method = System.qm_method()
        """

        self.qm_param = qm_param
        self.qm_method = qm_method
        self.qm_molecule = qm_molecule
        self.qm_atoms = qm_atoms
        self.qm_residues = qm_residues
        self.qm_charge_method = qm_charge_method

        self.mm_pdb_file = mm_pdb_file
        self.mm_forcefield = mm_forcefield
        self.mm_forcefield_water = mm_forcefield_water
        self.mm_nonbond_method = mm_nonbond_method
        self.mm_nonbond_cutoff = mm_nonbond_cutoff
        self.mm_constraints = mm_constraints
        self.mm_temp = mm_temp

        self.mm_pdb = None
        self.qm_energy = None
        self.mm_Te = None
        self.mm_Ke = None
        self.mm_tot_energy = None
        self.mod_Te = None
        self.mod_Ke = None
        self.mm_qm_energy = None
        self.mm_positions = None

        self.embedding_method = embedding_method

        if self.mm_pdb_file is not None:
            self.mm_pdb = ow.create_openmm_pdb(self.mm_pdb_file)

    def make_modeller(self, keep_qm=False):
        if self.mm_pdb is not None:
            modeller =  ow.create_openmm_modeller(self.mm_pdb)
        if keep_qm is False: 
            ow.delete_atoms(modeller, self.qm_atoms)
        else: 
            ow.keep_atoms(modeller, self.qm_atoms)

        return modeller

    def make_zero_energy():
        pass

#    def get_modeller_state_info(self):
#        """
#        Get the MM energy of a modeller object created with OpenMM
#        *currently only gives energy 
#        *but can now get positions and forces if needed
#        """
#        # Create an OpenMM system from object's modeller topology
#        self.mod_openmm_sys = ow.create_openmm_system(self.modeller.topology)
#
#        # Create an OpenMM simulation from the openmm system,
#        # modeller topology and positions.
#        self.mod_openmm_sim = ow.create_openmm_simulation(self.mod_openmm_sys,
#                                                          self.modeller.topology,
#                                                          self.modeller.positions)
#
#        # Calls openmm wrapper to get the kinetic and potential
#        # energy of the state
#        state = ow.get_state_info(self.mod_openmm_sim, energy=True, positions=True, forces=True)
#
## Converts the energy values from kj mol^-1 to au and stores in self
#self.mod_Te = state['potential'] 
#self.mod_Ke = state['kinetic'] 
#
#self.mod_positions = state['positions'] 
#
## forces currently in kJ/(mol nm) might need to convert to another unit later
#self.mod_forces = state['forces'] 


#    def get_openmm_state_info(self):
#        """
#        Get state information of a MM system described in a given pdb
#        *currently only gives energy and positions
#        """
#
#        # Create an OpenMM system from object's pdb topology
#        self.openmm_sys = ow.create_openmm_system(self.mm_pdb.topology)
#        
#        # Create an OpenMM simulation from the openmm system,
#        # pdb topology and positions
#        self.openmm_sim = ow.create_openmm_simulation(self.openmm_sys,
#                                                      self.mm_pdb.topology,
#                                                      self.mm_pdb.positions)
#
#        # Calls openmm wrapper to get specified info 
#        state = ow.get_state_info(self.openmm_sim, energy=True,positions=True)
#
#        # Converts the energy values from kj mol^-1 to au and stores in self
#        self.mm_Te = state['potential'] 
#        self.mm_Ke = state['kinetic'] 
#
#        self.mm_tot_energy = self.mm_Te + self.mm_Ke
#
#        # Converts the energy values from nm to angstroms and stores in self
#        self.mm_positions = state['positions'] 



#    def get_mm_qm_energy(self):
#        """
#        Get the MM energy of a user-defined QM system
#        """
#
#        # Create a modeller object of only the qm atoms
#        ow.keep_atoms(self.modeller, self.qm_atoms)
#
#        # Get the energy of the modeller system
#        self.get_modeller_state_info()
#
#        # Save the energy of the modeller object
#        self.mm_qm_energy = self.mod_Te + self.mod_Ke

    def additive(self):
        # TODO: need to work on this...
        # I don't know if having the if state about energy is None is smart...have to rethink
        """
        Gets energies of needed components and computes
        a qm/mm energy with a specified embedding scheme
        """

        # Get MM energy on MM region 
        # Create a modeller object of only the mm atoms
        mm_mm = self.make_modeller()
        # Get the energy of the modeller system
        mm_mm_sys, mm_mm_state = System.get_info(mm_mm, charges=True)

        # Get QM energy
        if self.embedding_method == 'Mechanical':
            self.qm_energy = pw.get_psi4_energy(self.qm_molecule,
                                                self.qm_param,
                                                self.qm_method)
            #get only nonbonded energy of whole system
            all_nb_sys, all_nb_state = System.get_info(self.mm_pdb, forces='nonbonded') 
        
            #get only nonbonded energy of qm region 
            mm_qm = self.make_modeller(keep_qm=True) 
            qm_nb_sys, qm_nb_state = System.get_info(mm_qm, forces='nonbonded')

            #get only nonbonded energy of mm region 
            mm_nb_sys, mm_nb_state = System.get_info(mm_mm, forces='nonbonded') 
    
            self.nb_energy = all_nb_state['energy'] - qm_nb_state['energy'] - mm_nb_state['energy']

        # at somepoint put if qm_energy is None qualifier? 
        if self.embedding_method == 'Electrostatic':
            self.qm_energy = pw.get_psi4_energy(self.qm_molecule,
                                                self.qm_param,
                                                self.qm_method,
                                                'Electrostatic',
                                                mm_mm_state['charge'],
                                                mm_mm_state['positions'])

            # need to implement no LJ stuff
            self.nb_energy = 0 


        self.qmmm_energy = mm_mm_state['energy'] + self.qm_energy + self.nb_energy


    def subtractive(self):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a specified embedding scheme
        """
        # Get MM energy on whole system
        all_sys, all_state = System.get_info(self.mm_pdb)
        print(all_state['energy'])

        # Get MM energy on QM region
        mm_qm = self.make_modeller(keep_qm=True) 
        mm_qm_sys, mm_qm_state = System.get_info(mm_qm)
        print(mm_qm_state['energy'])

        # Get QM energy
        if self.qm_energy is None and self.embedding_method == 'Mechanical':
            self.qm_energy = pw.get_psi4_energy(self.qm_molecule,
                                                self.qm_param,
                                                self.qm_method)
        # Compute the total QM/MM energy based on
        # subtractive Mechanical embedding
        if self.embedding_method == 'Mechanical':
            self.qmmm_energy = self.qm_energy \
                            + all_state['energy'] \
                            - mm_qm_state['energy']

    def make_qm_molecule(self):

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        if self.mm_positions is None:
            sys, state = System.get_info(self.mm_pdb)
            self.mm_positions = state['positions']

        for idx in self.qm_atoms:
            for atom in self.mm_pdb.topology.atoms():
                if atom.index == idx:
                    x,y,z = self.mm_positions[idx][0], self.mm_positions[idx][1], self.mm_positions[idx][2]
                    out += line.format(atom.element.symbol, x, y, z)
        out += 'no_reorient \n '
        out += 'no_com \n '
        return out

    def get_info(pdb, forces=None, charges=False):

        # Create an OpenMM system from an object's topology
        system = ow.create_openmm_system(pdb.topology)

        if forces=='nonbonded':
            for i in range(3):
                system.removeForce(0)


        # Create an OpenMM simulation from the openmm system,
        # topology and positions.
        simulation = ow.create_openmm_simulation(system,
                                                pdb.topology, pdb.positions)

        # Calls openmm wrapper to get information specified
        state = ow.get_state_info(simulation,
                                energy=True, 
                                positions=True, 
                                forces=True)
        state['energy'] = state['potential'] + state['kinetic']
        if charges is True:
            state['charge'] = ow.get_sys_info(system)
            
        return system, state
