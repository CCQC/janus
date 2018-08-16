"""
QM wrapper super class
"""
from abc import ABC, abstractmethod
class QM_wrapper(ABC):

    def __init__(self, config, program):

        self._program = program
        self.config = config

        if 'qm_basis_set' in config:
            self.basis_set = config['qm_basis_set']
        else: 
            self.basis_set = 'STO-3G'
        if 'qm_scf_type' in config:
            self.scf_type = config['qm_scf_type']
        else:
            self.scf_type = 'df'
        if 'qm_guess' in config:
            self.guess_orbitals = config['qm_guess']
        else:
            self.guess_orbitals = 'sad'
        if 'qm_reference' in config:
            self.reference = config['qm_reference']
        else: 
            self.reference = 'rhf'
        if 'e_convergence' in config: 
            self.e_convergence = config['qm_e_convergence']
        else:
            self.e_convergence = 1e-8
        if 'd_convergence' in config: 
            self.d_convergence = config['qm_d_convergence']
        else:
            self.d_convergence = 1e-8
        if 'qm_method' in config:
            self.method = config['qm_method']
        else:
            self.method = 'scf'
        if 'qm_atoms' in config:
            self.qm_atoms = config['qm_atoms']
        else:
            self.qm_atoms = []
        if 'qm_residues' in config:
            self.qm_residues = config['qm_residues']
        else:
            self.qm_residues = []
        if 'charge_method' in config:
            self.charge_method = config['qm_charge_method']
        else: 
            self.charge_method = 'MULLIKEN_CHARGES'
    
    @abstractmethod
    def compute_qm(self):
        pass

    @abtractmethod
    def build_qm_param(self):
        pass

    def get_qm(self):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 
        NEED TO RENAME

        Parameters
        ----------
        None

        Returns
        -------
        A dictionary with energy and gradient information
        """
        self.build_qm_param()
        self.qm_info()
        return self.energy, self.gradient

    def set_qm_geometry(self, qm_geometry):

        self.qm_geometry = qm_geometry

    def set_external_charges(self, charges):
        
        self.charges = charges

