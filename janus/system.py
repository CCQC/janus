import numpy as np
from copy import deepcopy
from mendeleev import element

class System(object):
    """
    A system class that stores QM, MM, and QM/MM
    input parameters, as well as system information
    such as geometry, and energy
    psi4_mp2._system.
    """

    def __init__(self, qm_indices, run_ID, partition_ID='qm'):

        """
        Initializes system with None values for all parameters

        Parameters
        ----------

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
        self.qm_atoms = deepcopy(qm_indices)
        self.run_ID = run_ID
        self.partition_ID = partition_ID
        self.qm_positions = None
        self.buffer_groups = None
        self.switching_functions = []
        self.forces = None
        self.energy = None

        self.entire_sys = {}
        self.primary_subsys = {}
        self.second_subsys = {}
        self.boundary = {}


    def compute_scale_factor_g(self, qm, mm, link):
        '''
        Computes scale factor g for link atom, RC, and RCD schemes. 
        Note: r given in pmm but don't need to convert because it is a ratio
              need to get other ways to compute g factor
              need functionality for different link atoms

        Parameters
        ----------
        qm: string of element symbol of the QM atom involved in broken bond 
        mm: string of element symbol of the MM atom involved in broken bond 
        link: string of element symbol for link atom

        Returns
        -------
        g, the scaling factor

        Examples
        --------
        compute_scale_factor(qm='C', mm='C', link='H')
        """
        '''
        
        r_qm = element(qm).covalent_radius_pyykko
        r_mm = element(mm).covalent_radius_pyykko
        r_link = element(link).covalent_radius_pyykko

        g = (r_qm + r_link)/(r_qm + r_mm)
        
        return g
        
    # need to add all the other parameters written into system e.g. energy


    def compute_COM(self, positions):
    
        xyz = np.zeros(3)
        M = 0

        for pos in positions:
            m = element(pos[0]).atomic_weight
            M += m
            xyz += m * np.array(pos[1])

        xyz *= 1/M
        
        return xyz



