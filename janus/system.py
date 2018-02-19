class System(object):
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    """

    def __init__(self, qmmm, qm, mm):

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

        self.qm_param = qm['qm_param']
        self.qm_method = qm['qm_method']
        self.qm_atoms = qm['qm_atoms']
        self.qm_residues = qm['qm_residues']
        self.qm_charge_method = qm['qm_charge_method']

    # need to add other mm parameters 
        self.mm_pdb_file = mm['mm_pdb_file']
        self.mm_forcefield = mm['mm_forcefield']
        self.mm_forcefield_water = mm['mm_forcefield_water']
        self.mm_nonbond_method = mm['mm_nonbond_method']
        self.mm_nonbond_cutoff = mm['mm_nonbond_cutoff']
        self.mm_constraints = mm['mm_constraints']
        self.mm_temp = mm['mm_temp']

        self.embedding_method = qmmm['embedding_method']

    # need to add all the other parameters written into system e.g. energy

