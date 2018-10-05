from abc import ABC, abstractmethod
import mendeleev as mdlv
"""
QM wrapper super class
"""

class QM_wrapper(ABC):

    def __init__(self, param, program):
        self.program = program
        self.param = param

        self.qm_param = None
        self.external_charges = None
        self.charges = None
        self.is_open_shelled = False

    @abstractmethod
    def build_qm_param(self):
        pass
    @abstractmethod
    def optimize_geometry(self):
        pass
    @abstractmethod
    def compute_energy_and_gradient(self):
        pass

    @abstractmethod
    def equilibrate(self):
        raise Exception('method not implemented for class')

    @abstractmethod
    def get_main_info(self):
        raise Exception('method not implemented for class')

    def get_energy_and_gradient(self, traj, include_coulomb='all', link_atoms=None, minimize=False, charges=None):
        """
        Gets the energy and gradient from a QM computation of the primary subsystem 

        Parameters
        ----------
        geometry : str
            Geometry information for the QM region
        total_elec : int 
            The total number of electrons present in the QM region

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
            self.set_external_charges(charges)

        if not self.qm_param:
            self.build_qm_param()

        if minimize is True:
            geom = self.optimize_geometry()
        else:
            self.compute_energy_and_gradient()

        self.info = {}
        self.info['energy'] = self.energy
        self.info['gradients'] = self.gradient
        
        return self.info

            
    def set_external_charges(self, charges):
        """
        Sets the charges due to the MM point charges
        as self.external_charges

        Parameters
        ----------
        charges: a list with the charge and position(in angstroms) for all 
                 particles outside the QM region 

        Returns
        -------
        None

        Examples
        --------
        set_external_charges([[-.45, 0.0, 0.0, 0.0]])
        set_external_charges(charges)
        """
        
        self.external_charges = charges

    def get_qm_geometry(self, qm_traj=None):
        """
        Uses the atoms and positions from a MDtraj trajectory object
        with just the qm region to obtain the geometry information

        Parameters
        ----------
        qm_traj: a MDtraj object describing just the primary subsystem,
                 default is None

        Returns
        -------
        out, total
        out: the str with geometry information in angstroms
        total: total number of electrons in the primary subsystem

        Examples
        --------
        geom, total_elec = qm_positions()
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

