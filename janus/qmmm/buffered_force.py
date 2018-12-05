from janus.qmmm import AQMMM
from janus.system import System
import numpy as np
from copy import deepcopy

class BufferedForce(AQMMM):
    """
    Class for the Buffered-Force adaptive QM/MM method.
    Inherits from AQMMM class

    Parameters
    ----------
        hl_wrapper : :class:`~janus.mm_wrapper.MMWrapper` subclass or :class:`~janus.qm_wrapper.QMWrapper` subclass
            Wrapper for performing the high-level computation. 
            Traditionally QM but user can define MM.
        ll_wrapper : :class:`~janus.mm_wrapper.MMWrapper` subclass
            Wrapper for performing the low-level computation
        sys_info : str 
            A string with the filename or a list with multiple filenames 
            that contain position and topology information. 
        sys_info_format : str 
            Describes what kind of input is contained in sys_info. Default is pdb.
        qm_center: list 
            Atoms that define the qm center, default is [0].
            If more than one index is given, COM is used as qm_center
        partition_scheme: str 
            Scheme to use to define buffer groups,
            default is distance (only scheme available as of now)
        Rmin: float 
            Inner radius for distance partition in angstroms, default is 3.8
        Rmax: float 
            Outer radius for distance partition in angstroms, default is 4.5
        qmmm_param : dict
            A dictionary with any parameters to pass into the QMMM class.
            See QMMM class for specifics

    """

    def __init__(self, *args, **kwargs):

        super().__init__('Buffered-Force', *args, **kwargs)

    def find_configurations(self): 
        """
        Finds the partitions as required by the buffered-force method 
        and saves each partition as a system object.
        Saves all systems in the dictionary self.systems

        Parameters
        ----------
        qm_center : list 
            atoms that define the qm center, default is None

        """
    
        qm = System(qm_indices=self.qm_atoms, qm_residues=self.qm_residues, run_ID=self.run_ID, partition_ID='qm')
        qm.buffer_groups = self.buffer_groups

        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][qm.partition_ID] = qm

        # the following only runs if there are groups in the buffer zone
        if self.buffer_groups:
            qm.original_qm_residues = deepcopy(qm.qm_residues)
            for i, buf in self.buffer_groups.items():
                qm.qm_residues.append(i)
                for idx in buf.atoms:
                    qm.qm_atoms.append(idx)

            # qm has a copy of its buffer groups - 
            qm.buffer_groups = self.buffer_groups

    def run_qmmm(self, main_info, wrapper_type):
        """
        """

        self.update_traj(main_info['positions'], main_info['topology'], wrapper_type)
        self.find_buffer_zone()
        self.find_configurations()
        
        qm = self.systems[self.run_ID]['qm']

        self.qm_atoms = deepcopy(qm.qm_atoms)

        self.electrostatic(qm, main_info)

        print('Interpolating QM/MM partitions')
        self.run_aqmmm(qm)
        self.systems[self.run_ID]['kinetic_energy'] = main_info['kinetic']

        # updates current step count
        self.run_ID += 1

        # delete the information of 2 runs before, only save current run and previous run information at a time
        if self.run_ID > 1:
            del self.systems[self.run_ID - 2]

    def electrostatic(system, main_info):

        if self.qmmm_scheme == 'subtractive':
            # Get MM energy on whole system
            system.entire_sys = deepcopy(main_info)
            
            # Get QM energy and gradient
            charges = self.get_external_charges(system)
            system.primary_subsys['hl'] = self.hl_wrapper.get_energy_and_gradient(traj_ps, charges=charges)

    def run_aqmmm(system)

        qm_grad = system.primary_subsys['hl']['gradients']
        qmmm_force = {}
        # print('sys qm_atom', system.qm_atoms)
            
        # iterate over list of qm atoms
        for i, atom in enumerate(self.qm_atoms):

            # compute the qmmm gradient for the qm atoms: 
            # mm_entire - mm_entire + qm
            qmmm_force[atom] = np.zeros(3)
            # these are in units of au_bohr, convert to openmm units in openmm wrapper
            qmmm_force[atom] += -1 * (-system.entire_sys['gradients'][i] + qm_grad[i])
            
            # treating gradients for link atoms
            if self.qmmm_boundary_bonds:
                raise Exception('Buffered-Force method currently cannot treat link atoms')
                                
        system.qmmm_forces = qmmm_force

            

        
