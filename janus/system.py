from . import psi4_wrapper as pw
from . import openmm_wrapper as ow


class System:
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    """

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

        if self.mm_pdb is not None:
            modeller = ow.create_openmm_modeller(self.mm_pdb)
        if keep_qm is False:
            ow.delete_atoms(modeller, self.qm_atoms)
        else:
            ow.keep_atoms(modeller, self.qm_atoms)

        return modeller

    def make_zero_energy():
        pass

    def additive(self):
        # TODO: need to work on this...
        # I don't know if having the if state about energy is
        # None is smart...have to rethink
        """
        Gets energies of needed components and computes
        a qm/mm energy with a specified embedding method using
        an additive scheme
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
            # get only nonbonded energy of whole system
            all_nb_sys, all_nb_state = System.get_info(self.mm_pdb,
                                                       forces='nonbonded')

            # get only nonbonded energy of qm region
            mm_qm = self.make_modeller(keep_qm=True)
            qm_nb_sys, qm_nb_state = System.get_info(mm_qm, forces='nonbonded')

            # get only nonbonded energy of mm region
            mm_nb_sys, mm_nb_state = System.get_info(mm_mm, forces='nonbonded')

            self.nb_energy = all_nb_state['energy'] \
                             - qm_nb_state['energy'] \
                             - mm_nb_state['energy']

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

        self.qmmm_energy = mm_mm_state['energy']\
                           + self.qm_energy\
                           + self.nb_energy

    def subtractive(self):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a subtractive mechanical embedding scheme
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

        """
        Extracts positions of qm atoms from a pdb file to make a string of
        qm atom geometry

        Parameters
        ----------
        None

        Returns
        -------
        A string that can be feed into Psi4 to make a molecule

        Examples
        --------
        qm = self.make_qm_molecule()
        """

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        if self.mm_positions is None:
            sys, state = System.get_info(self.mm_pdb)
            self.mm_positions = state['positions']

        for idx in self.qm_atoms:
            for atom in self.mm_pdb.topology.atoms():
                if atom.index == idx:
                    x, y, z = self.mm_positions[idx][0],\
                              self.mm_positions[idx][1],\
                              self.mm_positions[idx][2]
                    out += line.format(atom.element.symbol, x, y, z)
        out += 'no_reorient \n '
        out += 'no_com \n '
        return out

    def get_info(pdb, forces=None, charges=False):
        """
        Gets information about a system, e.g., energy, positions, forces

        Parameters
        ----------
        pbd : a OpenMM pdb object or OpenMM modeller object
        forces : if forces=='nonbonded', any bonded forces are not considered
                 Default is None
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
        system = ow.create_openmm_system(pdb.topology)

        # Remove Bond, Angle, and Torsion forces to leave only nonbonded forces
        if forces == 'nonbonded':
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
