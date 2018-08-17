from copy import deepcopy
import numpy as np
import mdtraj as md
from .system import System
"""
QMMM class for QMMM computations
"""
class QMMM(object):

    def __init__(self, config, qm_wrapper):
        
        self.system = System()
        self.qm_wrapper = qm_wrapper
        self.qm_geometry = None

        self.traj = md.load(config['mm_pdb_file'])
        self.topology = self.traj.topology
        self.positions = self.traj.xyz

        if 'qm_atoms' in config:
            self.qm_atoms = config['qm_atoms']
        else:
            self.qm_atoms = []

        if 'scheme' in config:
            self.qmmm_scheme = config['scheme']
        else: 
            self.qmmm_scheme = 'subtractive'

        if 'embedding_method' in config:
            self.embedding_method = config['embedding_method']
        else:
            self.embedding_method = 'Mechanical'

        if 'boundary_treatment' in config:
            self.boundary_treatment = config['boundary_treatment']
        else:
            self.boundary_treatment = 'link_atom'

        if 'link_atom' in config:
            self.link_atom_element = config['link_atom']
        else:
            self.link_atom_element = 'H'

        
    def additive(self, mm_wrapper):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a specified embedding method using
        an additive scheme
        """

        #need to add if these things are none then do the following?
        # maybe not because already checks in mm_wrapper functions

        # Get MM energy on MM region
        self.second_subsys = mm_wrapper.get_second_subsys()

        # Get non coulomb MM energy on PS-SS interaction
        self.boundary = mm_wrapper.get_boundary(coulomb=False)

        # Get any link atom information
        self.boundary_info = mm_wrapper.get_boundary_info()

        # Get QM energy
        # get QM positions from pdb
        if self.qm_positions is None:
            self.qm_positions = mm_wrapper.get_qm_positions() 
        self.qm = self.qm_wrapper.get_qm(self.qm_positions)

        # Compute total QM/MM energy based on additive scheme
        self.qmmm_energy = self.second_subsys['energy']\
                      + self.boundary['energy']\
                      + self.qm['energy']

        # Compute QM/MM gradients 
        qmmm_gradients = self.compute_gradients(scheme='additive')

    def subtractive(self, mm_wrapper):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a subtractive mechanical embedding scheme
        """

        # Get MM energy on whole system
        # should be updating as we go
        self.entire_sys = mm_wrapper.main_info

        # Get MM energy on QM region
        self.primary_subsys = mm_wrapper.get_primary_subsys(link=True)

        # Get position and identity of link atom for QM computation if relevant
        self.boundary_info = mm_wrapper.get_boundary_info()

        # Get QM energy
        self.qm_geometry = self.get_qm_positions() 
        self.qm = self.qm_wrapper.get_qm(self.qm_positions)

        # Compute the total QM/MM energy based on
        # subtractive Mechanical embedding
        self.qmmm_energy = self.entire_sys['energy']\
                      - self.primary_subsys['energy']\
                      + self.qm['energy']


    def compute_gradients(self, scheme='subtractive'):
        # NEED TO MAKE SURE: am I working with GRADIENTS or FORCES? NEED TO MAKE SURE CONSISTENT!
        # NEED TO MAKE SURE UNITS CONSISTENT

        if scheme == 'subtractive':

            ps_mm_grad, qm_grad = self.primary_subsys['gradients'], self.qm['gradients']
            print('mm gradients', ps_mm_grad)
            print('qm gradients', qm_grad)
            #qmmm_grad = np.zeros((len(all_mm_grad),3))
            qmmm_force = {}
                
            # iterate over list of qm atoms
            for i, atom in enumerate(self.qm_atoms):

                # compute the qmmm gradient for the qm atoms: 
                # mm_entire - mm_primary - qm
                qmmm_force[atom] = np.zeros(3)
                # these are in units of au_bohr, convert to openmm units in openmm wrapper
                qmmm_force[atom] += -1 * (- ps_mm_grad[i] + qm_grad[i])
                
                # treating gradients for link atoms
                if self.boundary_info:
                    q1 = int(self.boundary_info[0]['qm_id']) - 1
                    m1 = int(self.boundary_info[0]['mm_id']) - 1
                    g = self.boundary_info[0]['g_factor']
                    if atom == q1:
                        if self.boundary_treatment == 'link_atom':
                            # Project forces of link atoms onto the mm and qm atoms of the link atom bond
                            qmmm_force[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                            qmmm_force[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]
                            
                        # Forces on M2 requires forces on point charges which I'm not sure about so need to double check
                        if self.boundary_treatment == 'RC' or self.boundary_treatment == 'RCD':
                            qmmm_force[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                            qmmm_force[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]

            self.qmmm_forces = qmmm_force
        
    def get_info(self, scheme, mm_wrapper, partition=None):

        if not partition:
            if scheme =='subtractive':

                self.subtractive(mm_wrapper)
                self.compute_gradients(scheme='subtractive')
                
            if scheme == 'additive':
                print("Additive scheme needs some work and is not available yet") 

        if partition:
            # update relevant info for each partition
            self.qm_positions = partition.qm_positions
            mm_wrapper._system.qm_atoms = partition.qm_atoms
            self.qm_atoms = partition.qm_atoms

            if scheme =='subtractive':

                self.subtractive(mm_wrapper)
                self.compute_gradients(scheme='subtractive')
                
            if scheme == 'additive':
                print("Additive scheme needs some work and is not available yet") 
            

    def get_qm_geometry(self, qm_atoms=None, as_string=True):
        """
        TODO:
        1. need to phase out getting qm_positions in the openmm wrapper
        2. need to phase out getting link atom stuff through the openmm wrapper
        - MDtraj can do ALL - just need to convert to openmm trajectory
        - this way would be more general and robust - the only thing is to make sure the qmmm 
         only things still work - not just with aqmmm
        """
        """
        Grabs the positions of the atoms in the primary subsystem from self._positions
        and makes a string with the element and xyz geometry coordinates. Adds the link atom positions
        when link atoms are needed.
        Note: In a MD time step, does the position update or not? Need to make sure this updates

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        qm_positions()
        """
        if qm_atoms is None:
            qm_atoms = self.qm_atoms
            

        if as_string is True:
            out = ""
            line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '
        if as_string is False:
            out = []

        for idx in qm_atoms:
            x, y, z =   self.traj.xyz[0][idx][0],\
                        self.traj.xyz[0][idx][1],\
                        self.traj.xyz[0][idx][2]

            symbol = self.traj.topology.atom(idx).element.symbol
            
            if as_string is True:
                # convert to angstroms
                out += line.format(symbol, x*10, y*10, z*10)

            else:
                out.append([symbol, [x, y, z]])
        ## if there are bonds that need to be cut
        #if self._boundary_bonds:
        #    # Need to add if statement for any treatments that don't need link atoms
        #    #if self._system.boundary_treatment !=
        #    for atom in self.link_atoms:
        #        pos = self.link_atoms[atom]['link_positions']*MM_wrapper.nm_to_angstrom
        #        x, y, z = pos[0], pos[1], pos[2]
        #        out += line.format(self.link_atoms[atom]['link_atom'], x, y, z)
        return out


    def find_boundary_bonds(self, qm_atoms=None):
        """
        Identified any covalent bonds that the QM/MM boundary cuts across

        Parameters
        ----------
        qm_atoms: A list of atom indicies corresponding to the atoms in
                  the primary subsystem. Default list is taken from qm_atoms
                  stored in the System object

        Returns
        -------
        A list of tuples corresponding to the qm and mm atoms (as OpenMM atom object) involved in every bond
        that need to be cut.

        Examples
        --------
        bonds = find_boundary_bonds()
        bonds = find_boundary_bonds(qm_atoms=[0,1,2,3])
        """

        if qm_atoms is None:
            qm_atoms = self.qm_atoms

        qm_atoms = self.edit_qm_atoms(qm_atoms)

        self.qmmm_boundary_bonds = []
        # determining if there are bonds that need to be cut
        for bond in self.topology.bonds:
            # find any bonds that involve the qm atoms
            # iterates over tuple of mdtraj atom objects
            if bond[0].index in qm_atoms or bond[1].index in qm_atoms:
                # isolate bonds that involve one in the qm atoms and one outside
                if bond[0].index not in qm_atoms or bond[1].index not in qm_atoms:
                    qm_atom = {}
                    mm_atom = {}
                    if bond[0].index in qm_atoms:
                        qm_atom = bond[0]
                        mm_atom = bond[1]
                    else:
                        qm_atom = bond[1]
                        mm_atom = bond[0]
                    self.qmmm_boundary_bonds.append((qm_atom, mm_atom))

    def edit_qm_atoms(self, qm_atoms=None, solvent='water')

        if qm_atoms is None:
            qm_atoms = self.qm_atoms

        top = self.topology
        residues = [] 
        for i in qm_atoms:
            idx = top.atom(i).residue.index
            # make sure just go through each residues once
            if idx not in residues:
                residues.append(idx)
                res = top.residue(idx)
                if res.is_water:
                ## add any hydrogens that have oxygen inside of qm region 
                    if top.atom(i).element.symbol == 'O':
                        for a in res.atoms:
                            if (a.element.symbol =='H' and a.index not in qm_atoms):
                                qm_atoms.append(a.index)
                ## remove any hydrogens that have oxygen outside of qm region from qm region
                    elif top.atom(i).element.symbol == 'H':
                        for a in res.atoms:
                            if (a.element.symbol =='O' and a.index not in qm_atoms):
                                qm_atoms.remove(i)

        qm_atoms.sort()
        return qm_atoms


    def prepare_link_atom(self, RC=False):
        """
        Saves the qm and mm atom associated with the bond being cut across the QM/MM boundary
        and computes the g scaling factor for the link atom.

        Parameters
        ----------
        RC: a bool specifying whether to find the indices of the atoms bonded to mm atom of 
            bond being cut. Default is false

        Returns
        -------
        None

        Examples
        --------
        prepare_link_atom(RC=True)
        """

        self.link_atoms = {}

        for i, bond in enumerate(self.qmmm_boundary_bonds):

            qm = bond[0]
            mm = bond[1]

            self.link_atoms[i] = {}

            # saving id because id does not change between sys and modeller
            self.link_atoms[i]['qm_atom'] = qm

            self.link_atoms[i]['mm_atom'] = mm

            self.link_atoms[i]['link_atom'] = self.link_atom_element
            g = self.system.compute_scale_factor_g(qm.element.symbol, mm.element.symbol, self.link_atom_element)
            self.link_atoms[i]['scale_factor'] = g 
            # this is in nm
            self.link_atoms[i]['link_positions'] = self.get_link_atom_position(self.positions[qm.index], self.positions[mm.index], g)

            # MAYBE CAN USE MDTRAJ FIND NEIGHBORS FOR THIS!!!!!!!!
            if RC is True:
                bonds = []
                # find index of atoms bonded to mm atom
                for bond in self.topology.bonds:
                    if bond[0]== mm or bond[1] == mm:
                        if bond[0] != qm and bond[1] != qm:
                            if bond[0] != mm:
                                bonds.append(bond[0].index)
                            elif bond[1] != mm:
                                bonds.append(bond[1].index)

                self.link_atoms[i]['bonds_to_mm'] = bonds


    def create_primary_subsys_trajectory(self, qm_atoms=None):
        '''
        Creates an OpenMM modeller object that includes any link atoms.
        Note: Currently adds a very specfic link atom, need to expand 

        Parameters
        ----------
        mod: modeller object of primary system without link atom added
        atom: dictionary containing information for the link atom 

        Returns
        -------
        A OpenMM modeller object

        Examples
        --------
        create_link_atom_modeller(mod=modeller, atom=link)
        '''

        if qm_atoms is None:
            qm_atoms = self.qm_atoms

        qm_atoms = self.edit_qm_atoms(qm_atoms)

        traj = self.trajectory.atom_slice(qm_atoms)

        traj.find_boundary_bonds(qm_atoms):

        if self.qmmm_boundary_bonds():
            self.prepare_link_atoms()

            for i, link in self.link_atoms.items():

                link_element = md.element.Element.getBySymbol(link['link_atom'])

                for atom in traj.topology.atoms:
                    if atom.serial == link['qm'].serial:
                        traj.topology.add_atom(name='H', element=link, residue=atom.residue, serial='link')

                        for atom2 in traj.topology.atoms:
                            if atom2.serial == 'link':
                                traj.topology.add_bond(atom2, atom)


            positions = np.append(positions[0], [link['link_positions']], axis=0)

