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

    def __init__(self, aqmmm = {}, qmmm={}, qm={}, mm={}, simulation = {}):

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

        if 'scheme' in qmmm:
            self.qmmm_scheme = qmmm['scheme']
        else: 
            self.qmmm_scheme = 'subtractive'

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
    
        if 'boundary_treatment' in qmmm:
            self.boundary_treatment = qmmm['boundary_treatment']
        else:
            self.boundary_treatment = 'link_atom'

        if 'link_atom' in qmmm:
            self.link_atom = qmmm['link_atom']
        elif self.boundary_treatment == 'link_atom':
            self.link_atom = 'H'

        if 'aqmmm_scheme' is aqmmm:
            self.aqmmm_scheme = aqmmm['aqmmm_scheme']
        else:
            self.aqmmm_scheme = 'ONIOM-XS'

        if 'steps' in simulation:
            self.steps = simulation['steps']
        else:
            self.steps = 1
        
        self.aqmmm = aqmmm

        self.build_qm_param()

        self.entire_sys = {}
        self.primary_subsys = {}
        self.second_subsys = {}
        self.boundary = {}
        self.qm = {}
        self.qm_positions = None
        self.boundary_info = {}


        
    def get_link_atom_position(self, qm_position, mm_position, g):

        return qm_position + g*(mm_position - qm_position)


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


class Partition(object):

    def __init__(self, indices, ID):

        self.qm_atoms = deepcopy(indices)
        self.ID = ID
        self.qm_positions = None
        self.forces = None
        self.energy = None
        self.buffer_groups = None
        self.switching_functions = []


    def compute_COM(self, positions):
    
        xyz = np.zeros(3)
        M = 0

        for pos in positions:
            m = element(pos[0]).atomic_weight
            M += m
            xyz += m * np.array(pos[1])

        xyz *= 1/M
        
        return xyz



            
        





        
