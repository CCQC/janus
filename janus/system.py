class System:
    """
    A system class that stores QM, MM, and QM/MM 
    input parameters, as well as system information 
    such as geometry, and energy
    """ 
    
    def __init__(self):
        """
        Initializes system with None values for all parameters 

        Parameters
        ----------
        None

        Returns
        -------
        A built system
        
        Examples
        --------
        qm_parameters = System.qm_param()
        method = System.qm_method()
        """
        
        self.qm_param = None
        self.qm_method = None
        self.qm_molecule = None
        self.qm_energy = None
    
