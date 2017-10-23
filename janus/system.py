class System:
    """
    A system class that stores QM, MM, and QM/MM 
    input parameters, as well as system information 
    such as geometry, and energy
    """ 
    
    def __init__(self, qm_param=None,qm_method=None,qm_molecule=None):
        """
        Initializes system with None values for all parameters 

        Parameters
        ----------
        qm_param: dict of qm parameters
        qm_method: str of desired qm method 
        qm_molecule: str of geometry

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
        self.qm_energy = None
    
