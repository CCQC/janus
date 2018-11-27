from abc import ABC, abstractmethod

class MMWrapper(ABC):
    """
    A super class for all molecular mechanics wrappers

    Note
    ----
    Since MMWrapper is a super class and has abstract methods
    the user cannot actually instantiate a MMWrapper object, but only its child objects
    """

    kjmol_to_au = 1/2625.5002 
    nm_to_angstrom = 10.0000000
    nm_to_bohr = 0.052917999999 
    kjmol_nm_to_au_bohr = kjmol_to_au*nm_to_bohr 
    au_bohr_to_kjmol_nm = 1/kjmol_nm_to_au_bohr

    def __init__(self, class_type,
                       sys_info=None, 
                       sys_info_format=None, 
                       **kwargs):
        self.system_info = sys_info
        self.system_info_format = sys_info_format
        self.class_type = class_type
        self.main_info = None
        self.main_charges = None

        self.start_qmmm =  0
        self.end_qmmm   =  0
        self.md_steps = 0
        self.md_ensemble = 'NVE'

        self.return_trajectory_interval = 0
        self.return_trajectory_filename = 'output.nc'
        self.trajectory_format = 'NetCDF'                                                              
        self.return_checkpoint_interval = 0                                                                     
        self.return_checkpoint_filename = 'checkpoint.chk'                                                      
        self.return_system = True,                                                                               
        self.return_system_filename = 'final.pdb'                                                               
        self.return_info = ["potentialEnergy", "kineticEnergy", "totalEnergy", "temperature"]
        self.return_info_interval = 0                                                                         
        self.return_info_filename = 'info.dat'                                                                         
        self.return_forces_filename = 'forces.pkl'                                                              
        self.return_forces_interval = 0                                                                         

        for k, v in kwargs.items():
            setattr(self, k, v)

        self.qmmm_steps = self.end_qmmm - self.start_qmmm

        if (type(self.md_steps) is list and type(self.md_ensemble) is list):
            self.md_ensemble = self.md_ensemble[-1]
            self.other_md_ensembles = self.md_ensemble[0:-1]
            self.other_ensemble_steps = self.md_steps[0:-1]
        elif (type(self.md_steps) is int and type(self.md_ensemble) is str):
            self.other_md_ensembles = None
            self.other_ensemble_steps = None

        if type(self.md_steps) is int:
            self.end_steps = self.md_steps - self.end_qmmm
        elif type(self.md_steps) is list:
            self.end_steps = self.md_steps[-1] - self.end_qmmm

        super().__init__()

    def get_energy_and_gradient(self, traj, geometry=None, include_coulomb='all', link_atoms=None, minimize=False, charges=None):
        """
        Gets the energy and gradient from a MM computation

        Parameters
        ----------
        traj : MDtraj trajectory object
        geometry : str
            A string containing geometry information as XYZ coordinates. Not applicable for MM programs
        include_coulomb : str

            whether to include coulombic interactions. 
            'all' (default) includes coulombic forces for all particles,
            'no_link' excludes coulombic forces for link atoms,
            'only' excludes all other forces for all atoms,
            'none' excludes coulombic forces for all particles.

        link_atoms : list
            indices of link_atoms
        minimize : bool
            whether to return the geometry optimized energy 
        charges : list
            charges and corresponding positions in angstroms as xyz coordinates

        Returns
        -------
        dict
            A dictionary with energy('energy') and gradient('gradients') information
             
        """

        topology, positions = self.convert_trajectory(traj)

        if charges is not None:
            self.set_external_charges(charges)

        info = self.compute_info(topology, positions, include_coulomb=include_coulomb, link_atoms=link_atoms, minimize=minimize)

        return info


    @abstractmethod
    def compute_info(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def set_external_charges(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def initialize(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def take_step(self, force):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def get_main_info(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def get_main_charges(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def convert_trajectory(self):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def convert_input(self):
        """
        Function implemented in individual child classes
        """
        pass
    
    @abstractmethod
    def set_up_reporters(self):
        """
        Function implemented in individual child classes
        """
        pass

    def get_geom_from_trajectory(self):
        """
        Function not implemented for MM wrappers
        """
        raise Exception('method not implemented for class')

    def set_qm_geometry(self):
        """
        Function not implemented for MM wrappers
        """
        raise Exception('method not implemented for class')

    def build_qm_param(self):
        """
        Function not implemented for MM wrappers
        """
        raise Exception('method not implemented for class')

    def optimize_geometry(self):
        """
        Function not implemented for MM wrappers
        """
        raise Exception('method not implemented for class')
