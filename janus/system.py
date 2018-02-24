class System(object):
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    psi4_mp2._system.
    """

    def __init__(self, qmmm={}, qm={}, mm={}):

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

        if 'qm_basis_set' in qm:
            self.qm_basis_set = qm['qm_basis_set']
        else: 
            self.qm_basis_set = 'STO-3G'
        if 'qm_scf_type' in qm:
            self.qm_scf_type = qm['qm_scf_type']
        else:
            self.qm_scf_type = 'df'
        if 'qm_guess' in qm:
            self.qm_guess = qm['qm_guess']
        else:
            self.qm_guess = 'sad'
        if 'qm_reference' in qm:
            self.qm_reference = qm['qm_reference']
        else: 
            self.qm_reference = 'rhf'
        if 'qm_e_convergence' in qm: 
            self.qm_e_convergence = qm['qm_e_convergence']
        else:
            self.qm_e_convergence = 1e-8
        if 'qm_d_convergence' in qm: 
            self.qm_d_convergence = qm['qm_d_convergence']
        else:
            self.qm_d_convergence = 1e-8
        if 'qm_method' in qm:
            self.qm_method = qm['qm_method']
        else:
            self.qm_method = 'scf'
        if 'qm_atoms' in qm:
            self.qm_atoms = qm['qm_atoms']
        else:
            self.qm_atoms = []
        if 'qm_residues' in qm:
            self.qm_residues = qm['qm_residues']
        else:
            self.qm_residues = []
        if 'qm_charge_method' in qm:
            self.qm_charge_method = qm['qm_charge_method']
        else: 
            self.qm_charge_method = 'MULLIKEN_CHARGES'

        if 'mm_pdb_file' in mm:
            self.mm_pdb_file = mm['mm_pdb_file']
        if 'mm_forcefield' in mm:
            self.mm_ff = mm['mm_forcefield']
        else:
            self.mm_ff = 'amber99sb.xml' 
        if 'mm_forcefield_water' in mm:
            self.mm_ff_water = mm['mm_forcefield_water']
        else:
            self.mm_ff_water = 'tip3p.xml'

        #need to tinker with these and figure out if specific to openmm
        if 'mm_nonbond_method' in mm:
            self.mm_nonbond_method = mm['mm_nonbond_method']
        if 'mm_nonbond_cutoff' in mm:
            self.mm_nonbond_cutoff = mm['mm_nonbond_cutoff']
        if 'mm_constraints' in mm:
            self.mm_constraints = mm['mm_constraints']

        if 'mm_temp' in mm:
            self.mm_temp = mm['mm_temp']
        else:
            self.mm_temp = float(300)
        if 'mm_fric_coeff' in mm:
            self.mm_fric_coeff = mm['mm_fric_coeff']
        else:
            self.mm_fric_coeff = float(1)
        if 'mm_step_size' in mm:
            self.mm_step_size = mm['mm_step_size']
        else: 
            self.mm_step_size = float(0.002)

        if 'embedding_method' in qmmm:
            self.embedding_method = qmmm['embedding_method']
        else:
            self.embedding_method = 'Mechanical'
        if 'qm_program' in qmmm:
            self.qm_program = qmmm['qm_program']
        else:
            self.qm_program = "Psi4"
        if 'mm_program' in qmmm:
            self.mm_program = qmmm['mm_program']
        else:
            self.mm_program = "OpenMM"
        # default integrator in openmm is the langevin integrator
    

        self.build_qm_param()

        self.entire_sys = {}
        self.primary_subsys = {}
        self.second_subsys = {}
        self.boundary = {}
        self.qm = {}
        self.qm_positions = None

    def build_qm_param(self):
        qm_param = {}
        qm_param['basis'] = self.qm_basis_set
        qm_param['scf_type'] = self.qm_scf_type
        qm_param['guess'] = self.qm_guess
        qm_param['reference'] = self.qm_reference
        qm_param['e_convergence'] = self.qm_e_convergence
        qm_param['d_convergence'] = self.qm_d_convergence
        
        self.qm_param = qm_param


        
    

    # need to add all the other parameters written into system e.g. energy

