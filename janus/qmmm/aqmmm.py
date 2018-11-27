from abc import ABC, abstractmethod
from copy import deepcopy
import mdtraj as md
import numpy as np
import mendeleev as mdlv
from janus.qmmm import QMMM
from janus.system import Buffer

class AQMMM(ABC, QMMM):
    """
    AQMMM super class for adaptive QMMM computations.
    Inherits from QMMM class.

    Note
    ----
    Since AQMMM is a super class and has abstract methods
    cannot actually instantiate AQMMM object, but only its child objects
    """

    nm_to_angstrom = 10.0000000

    def __init__(self, hl_wrapper, 
                       ll_wrapper, 
                       sys_info,
                       sys_info_format,
                       qmmm_param,
                       class_type):

        super().__init__(hl_wrapper, ll_wrapper, sys_info, sys_info_format=sys_info_format, **qmmm_param)

        self.class_type = class_type
        self.buffer_groups = {}
        self.compute_zero_energy()

    def run_qmmm(self, main_info, wrapper_type):
        """
        Drives QM/MM computation.
        Updates the positions and topology given in main_info,
        determines the partitions to be computed, for each partition,
        determines the QM/MM energy and gradients, then interpolates all
        partitions using specified adaptive QM/MM scheme

        Parameters
        ----------
        main_info : dict 
            Contains the energy, forces, topology, and position information 
            for the whole system
        wrapper_type : str
            Defines the program used to obtain main_info
        """

        self.update_traj(main_info['positions'], main_info['topology'], wrapper_type)
        self.partition()
            
        counter = 0
        for i, system in self.systems[self.run_ID].items():
            print('Running QM/MM partition {}'.format(counter))
            print('Number of QM atoms for partition {} is {}'.format(counter,len(system.qm_atoms)))

            self.qm_atoms = deepcopy(system.qm_atoms)

            print(self.embedding_method)
            if self.embedding_method =='Mechanical':
                self.mechanical(system, main_info)
            elif self.embedding_method =='Electrostatic':
                self.electrostatic(system, main_info)
            else:
                print('only mechanical and electrostatic embedding schemes implemented at this time')
            counter += 1

        print('QM/MM partitions done. Getting zero energies')
        self.get_zero_energy()
        print('Interpolating QM/MM partitions')
        self.run_aqmmm()
        self.systems[self.run_ID]['kinetic_energy'] = main_info['kinetic']
        #print('!qmmm_energy', self.systems[self.run_ID]['qmmm_energy'])
        #if self.run_ID % 10 == 0:
        print('!', self.run_ID, self.systems[self.run_ID]['qmmm_energy'] + self.systems[self.run_ID]['kinetic_energy'])

        # updates current step count
        self.run_ID += 1

        # delete the information of 2 runs before, only save current run and previous run information at a time
        if self.run_ID > 1:
            del self.systems[self.run_ID - 2]

    def define_buffer_zone(self, qm_center):
        """
        Determines buffer group atoms.
        Gets the buffer groups in the buffer zone based on a distance 
        partitioning scheme, saves each buffer group as a 
        :class:`~janus.system.Buffer` object,
        and saves all buffer groups in the dictionary self.buffer_groups.
        For water as a solvent, considers the whole water molecule as a buffer group.

        Note
        ----
        Not generalized for all solvents

        Parameters
        ----------
        qm_center : list 
            the indicies that make up the qm center

        """

        self.find_buffer_atoms(qm_center)
        
        self.edit_qm_atoms()
        print('edited qm atoms:')
        print(self.qm_atoms)

        # for adding identifying water buffer groups
        top = self.topology

        self.buffer_groups = {}
        print('buffer atoms',self.buffer_atoms)
        for i in self.buffer_atoms:
            print('processing buffer atom {}'.format(i))
            # since if a hydrogen is in buffer zone with link atoms the qm would be the same as qm_bz, and the center 
             # of mass would not be in the buffer zone
            # only if oxygen in buffer zone
            if (top.atom(i).residue.is_water and top.atom(i).element.symbol == 'O'):
                print('{} is O, adding to buffer groups'.format(i))
                idx = top.atom(i).residue.index
                if idx not in self.buffer_groups.keys():
                    buf = Buffer(ID=idx)
                    for a in top.residue(idx).atoms:
                        if a.index not in buf.atoms:
                            buf.atoms.append(a.index)
                    self.buffer_groups[idx] = buf

        if self.buffer_groups:
            self.get_buffer_info()

    def find_buffer_atoms(self, qm_center):
        """
        Find the buffer groups whose COM falls in between Rmin and Rmax

        Parameters
        ----------
        qm_center : list 
            the indicies that make up the qm center

        """

        if len(qm_center) == 1:
            self.COM_as_qm_center = False
            self.qm_center_xyz = self.traj.xyz[0][qm_center[0]]
            qm_center_idx = qm_center
            temp_traj = self.traj
            self.qm_center_weight_ratio = {qm_center[0] : 1}
        else:
            self.COM_as_qm_center = True
            self.qm_center_xyz, self.qm_center_atom_weights, self.qm_center_weight_ratio = self.compute_COM(qm_center)

            t = deepcopy(self.traj)
            t.topology.add_atom('DUM', md.element.Element.getBySymbol('H'), t.topology.atom(0).residue, 1)
            for atom in t.topology.atoms:
                if atom.name == 'DUM':
                    qm_center_idx = [atom.index]
            t.xyz = np.append(t.xyz[0], [self.qm_center_xyz], axis=0)
            temp_traj = t

        if self.partition_scheme == 'distance': 

            rmin_atoms = md.compute_neighbors(temp_traj, self.Rmin/10, qm_center_idx)
            rmax_atoms = md.compute_neighbors(temp_traj, self.Rmax/10, qm_center_idx)
            self.buffer_atoms = np.setdiff1d(rmax_atoms, rmin_atoms)
            print('buffer atoms identified by find_buffer_atom function:')
            print(self.buffer_atoms)
            self.qm_atoms = rmin_atoms[0].tolist()

            if self.COM_as_qm_center is False:
                self.qm_atoms.append(qm_center[0])

            print('qm_atoms identified by the find_buffer_atom function: ' )
            print(self.qm_atoms)
        
    def get_buffer_info(self):
        """
        Computes switching function parameters for a particular buffer group.
        Computes the distance r_i for each buffer group i
        between the qm_center and the COM of buffer group i
        and then computes lamda_i and d_lamda_i for that buffer group

        """

        self.buffer_distance = {}

        for i, buf in self.buffer_groups.items():

            xyz, buf.atom_weights, buf.weight_ratio = self.compute_COM(buf.atoms)
            buf.COM_coord = xyz
            buf.r_i = np.linalg.norm(buf.COM_coord - np.array(self.qm_center_xyz))*AQMMM.nm_to_angstrom
            self.buffer_distance[i] = buf.r_i
            buf.s_i, buf.d_s_i = self.compute_lamda_i(buf.r_i)

    def compute_COM(self, atoms):
        """
        Computes the center of mass of a specified group

        Parameters
        ----------
        atoms : list 
            indices defining the group to compute the COM for

        Returns
        -------
        numpy array
            COM xyz coordinates 
        dict
            the indices of each atom in the group and the atomic weight for each atom
        dict
            the indices of each atom in the group and the weight ratio of each atom to
            the total weight of the group (sum of all atomic weights)

        """
        
        xyz = np.zeros(3)
        M = 0

        atom_weight = {}
        weight_ratio = {}
        for i in atoms:

            symbol = self.traj.topology.atom(i).element.symbol
            m = mdlv.element(symbol).atomic_weight
            # this gives positions in nm
            position =  np.array(self.traj.xyz[0][i])

            M += m
            xyz += m * position
            atom_weight[i] = m

        for i in atoms:
            weight_ratio[i] = atom_weight[i]/M
            
        xyz *= 1/M
        
        return xyz, atom_weight, weight_ratio

    def compute_lamda_i(self, r_i):
        """
        Computes the switching function and the derivative 
        of the switching function defined as a 5th order spline:
        
        .. math::
            \lambda_i = -6x^5 + 15x^4 - 10x^3 + 1 
            d_{\lambda_i} = -30x^4 + 60x^3 - 30x^2

        where x is the reduced distance (r_i - rmin)/(rmax - rmin)
        of buffer group i

        Parameters
        ----------
        r_i : float 
            the distance between the qm center and the COM (in angstroms)

        Returns
        -------
        float
            lamda_i, unitless
        float
            derivative of lamda_i, unitless

        """

        print(r_i)
        print(self.Rmin)
        print(self.Rmax)
        x_i = float((r_i - self.Rmin) / (self.Rmax - self.Rmin))

        if (x_i < 0 or x_i > 1):
            raise ValueError("reduced distance x_i has to be between 0 and 1")

        lamda_i = -6*((x_i)**5) + 15*((x_i)**4) - 10*((x_i)**3) + 1

        d_lamda_i = (-30*(x_i)**4  + 60*(x_i)**3 - 30*(x_i)**2)
        d_lamda_i *= 1/(r_i * (self.Rmax - self.Rmin))

        return lamda_i, d_lamda_i

    def get_Rmin(self):
        """
        Function to return self.Rmin

        Returns
        -------
        float
            the distance from qm center to inner limit of buffer zone in angstroms

        """
        return self.Rmin

    def get_Rmax(self):
        """
        Function to return self.Rmin

        Returns
        -------
        float
            the distance from qm center to outer limit of buffer zone in angstroms

        """

        return self.Rmax

    def set_Rmin(self, Rmin):
        """
        Function to set self.Rmin
        
        Parameters
        ----------
        Rmin : float 
            the distance from qm center to inner limit of buffer zone in angstroms
        
        """

        self.Rmin = Rmin

    def set_Rmax(self, Rmax):
        """
        Function to set self.Rmax
        
        Parameters
        ----------
        Rmax : float 
            the distance from qm center to outer limit of buffer zone in angstroms
        
        """
        self.Rmax = Rmax

    def compute_zero_energy(self):
        """
        Compute the energy of the isolated groups at their minimum geometry
        """
 
        # for explicit solvent systems can just do once, but for bond forming/breaking processes
        # need to update??
        # this is only functional for explicitly solvated systems
        
        # get all the unique groups 
        residues = {}
        for res in self.topology.residues:
            if res.name not in residues:
                residues[res.name] = []
                for atom in res.atoms:
                    residues[res.name].append(atom.index)

        self.qm_zero_energies = {}
        self.mm_zero_energies = {}
        for res in residues:

            traj = self.traj.atom_slice((residues[res]))

            mm = self.ll_wrapper.get_energy_and_gradient(traj, minimize=True)
            self.mm_zero_energies[res] = mm['energy']

            qm = self.hl_wrapper.get_energy_and_gradient(traj, minimize=True)
            self.qm_zero_energies[res] = qm['energy']
    
    def get_zero_energy(self):
        """
        Incorporates the zero energy of groups to the total qmmm energy
        """

        for i, sys in self.systems[self.run_ID].items():
            for res in self.topology.residues:
                if res.index in sys.qm_residues:
                    sys.zero_energy += self.qm_zero_energies[res.name]
                else:
                    sys.zero_energy += self.mm_zero_energies[res.name]

            # maybe I should save a separate copy of qmmm energy somewhere
            sys.qmmm_energy -= sys.zero_energy
                    
    @abstractmethod
    def partition(self, info):
        """
        Function implemented in individual child classes
        """
        pass

    @abstractmethod
    def run_aqmmm(self):
        """
        Function implemented in individual child classes
        """
        pass


