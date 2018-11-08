from abc import ABC, abstractmethod
import mendeleev as mdlv

class QMWrapper(ABC):

    def __init__(self, param, program):
        """
        QM wrapper super class

        Parameters
        ----------
        param : dict 
            parameters for QM computations

        program : str
            what program to use for QM computations

        """
        self.program = program
        self.param = param

        self.qm_param = None
        self.external_charges = None
        self.charges = None
        self.is_open_shelled = False


    def get_energy_and_gradient(self, traj, include_coulomb='all', link_atoms=None, minimize=False, charges=None):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 

        Parameters
        ----------
        traj : MDtraj trajectory object
        include_coulomb : str
            whether to include coulombic interactions. Not applicable for QM programs

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

        Examples
        --------
        >>> run_qm(geom, 10)
        """
        self.get_qm_geometry(traj)

        if charges is not None:
            self.external_charges = charges

        if not self.qm_param:
            self.build_qm_param()

        if minimize is True:
            geom = self.optimize_geometry()
        else:
            self.compute_info()

        self.info = {}
        self.info['energy'] = self.energy
        self.info['gradients'] = self.gradient
        
        return self.info

            
    def get_qm_geometry(self, qm_traj=None):
        """
        Uses the atoms and positions from a MDtraj trajectory object
        with just the qm region to obtain the geometry information

        Parameters
        ----------
        qm_traj : MDtraj object
             describes just the primary subsystem, default is None

        Returns
        -------
        str
        geometry information in angstroms
        int
        total number of electrons in the primary subsystem

        """

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        self.total_elec = 0.0

        for i in range(qm_traj.n_atoms):
            x, y, z =   qm_traj.xyz[0][i][0],\
                        qm_traj.xyz[0][i][1],\
                        qm_traj.xyz[0][i][2]

            symbol = qm_traj.topology.atom(i).element.symbol
            n = mdlv.element(symbol).atomic_number
            self.total_elec += n
            
            out += line.format(symbol, x*10, y*10, z*10)

        self.qm_geometry = out

        if self.total_elec % 2 != 0:
            self.total_elec += self.charge   # takes charge into account
            if self.total_elec % 2 != 0:
                self.is_open_shelled = True

    @abstractmethod
    def compute_info(self):
        pass

    @abstractmethod
    def build_qm_param(self):
        pass
    @abstractmethod
    def optimize_geometry(self):
        pass

    def get_main_info(self):
        raise Exception('method not implemented for class')

    def set_external_charges(self):
        raise Exception('method not implemented for class')

    def initialize(self):
        raise Exception('method not implemented for class')

    def take_step(self, force):
        raise Exception('method not implemented for class')

    def get_main_charges(self):
        raise Exception('method not implemented for class')

    def convert_trajectory(self):
        raise Exception('method not implemented for class')

    def convert_input(self):
        raise Exception('method not implemented for class')

    def set_up_reporters(self):
        raise Exception('method not implemented for class')
    
