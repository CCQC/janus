from . import psi4_wrapper as pw
from . import openmm_wrapper as ow

class System:
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    """
    kjmol_to_au = 1/2625.5

    def __init__(self, qm_param=None, qm_method='scf', qm_molecule=None,
                 qm_charge_method='MULLIKEN_CHARGES',
                 qm_atoms=None, qm_residues=None, mm_pdb_file=None,
                 mm_forcefield = 'amber99sb.xml', mm_forcefield_water='tip3p.xml',
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
        self.mm_potential_e = None
        self.mm_kinetic_e = None
        self.mm_system = None

        self.embedding_method = embedding_method

        if self.mm_pdb_file is not None:
            self.mm_pdb = ow.create_openmm_pdb(self.mm_pdb_file)
        if self.mm_pdb is not None:
            self.modeller = ow.create_openmm_modeller(self.mm_pdb)

    def make_zero_energy():
        pass

    
    def get_openmm_energy_from_modeller(self):
    """
    Get the MM energy of a modeller object created with OpenMM
    """
        # Create an OpenMM system from object's modeller topology 
        self.mod_openmm_sys = ow.create_openmm_system(self.modeller.topology)

        # Create an OpenMM simulation from the openmm system, modeller topology and positions.
        self.mod_openmm_sim = ow.create_openmm_simulation(self.mod_openmm_sys,
                                                          self.modeller.topology,
                                                          self.modeller.positions)

        # Calls openmm wrapper to get the kinetic and potential energy of the state
        T, K = ow.get_state_info(self.mod_openmm_sim)

        # Converts the energy values from kj mol^-1 to au and stores in self
        self.mod_Te = T._value * System.kjmol_to_au
        self.mod_Ke = K._value * System.kjmol_to_au 

    def get_openmm_energy(self):
        """
        Get MM energy of a MM system described in a given pdb 
        """

        # Create an OpenMM system from object's pdb topology 
        self.openmm_sys = ow.create_openmm_system(self.mm_pdb.topology)

        # Create an OpenMM simulation from the openmm system, pdb topology and positions.
        self.openmm_sim = ow.create_openmm_simulation(self.openmm_sys,
                                                          self.mm_pdb.topology,
                                                          self.mm_pdb.positions)

        # Calls openmm wrapper to get the kinetic and potential energy of the state
        T, K = ow.get_state_info(self.openmm_sim)

        # Converts the energy values from kj mol^-1 to au and stores in self
        self.mm_Te = T._value * System.kjmol_to_au
        self.mm_Ke = K._value * System.kjmol_to_au 
        
        self.mm_tot_energy = self.mm_Te + self.mm_Ke

    def get_mm_qm_energy(self): 
        """
        Get the MM energy of a user-defined QM system
        """

        # Create a modeller object of only the qm atoms
        ow.keep_atoms(self.modeller, self.qm_atoms)
        
        # Get the energy of the modeller system
        self.get_openmm_energy_from_modeller()

        # Save the energy of the modeller object
        self.mm_qm_energy = self.mod_Te + self.mod_Ke 

    def get_qmmm_energy(self):
        """
        Gets energies of needed components and computes 
        a qm/mm energy with a specified embedding scheme
        """

        # Get MM energy on whole system
        self.get_openmm_energy()

        # Get MM energy on QM region
        self.get_mm_qm_energy()

        # Get QM energy 
        self.qm_energy = pw.get_psi4_energy(self.qm_molecule, self.qm_param, self.qm_method)

        # Compute the total QM/MM energy based on subtractive Mechanical embedding
        if self.embedding_method == 'Mechanical':
            self.qmmm_energy = self.qm_energy + self.mm_tot_energy - self.mm_qm_energy

    def make_qm_molecule():
        pass
