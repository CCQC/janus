from copy import deepcopy
import numpy as np
import mdtraj as md
from .system import System
"""
QMMM class for QMMM computations
"""
class QMMM(object):

    def __init__(self, config, qm_wrapper, mm_wrapper):
        
        self.qm_wrapper = qm_wrapper
        self.mm_wrapper = mm_wrapper
        self.qm_geometry = None
        self.run_ID = 0

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

        self.systems = {}

    def run_qmmm(main_info):

        self.update_traj(main_info['positions'], main_info['pdb'])

        system = System(self.qm_atoms, self.run_ID)

        if self.embedding_method =='Mechanical':
            self.mechanical(system, main_info)
        elif self.embedding_method =='Electrostatic':
            self.electrostatic(system, main_info)
        else:
            print('only mechanical and electrostatic embedding schemes implemented at this time')
            
        self.systems.[self.run_ID] = {}
        self.systems[self.run_ID][system.partition_ID] = system
        self.systems[self.run_ID]['qmmm_forces'] = system.qmmm_forces
        self.systems[self.run_ID]['qmmm_energy'] = system.qmmm_energy
        self.run_ID += 1


    def update_traj(self, position, topology):
        
        # later can think about saving instead of making new instance
        # convert openmm topology to mdtraj topology
        if mm_wrapper.program == 'OpenMM':
            top = md.Topology.from_openmm(topology)
        self.traj = md.Trajectory(position, top)


    def mechanical(self, system, main_info):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a subtractive mechanical embedding scheme
        """
        if self.qmmm_scheme == 'subtractive':
            # Get MM energy on whole system
            system.entire_sys = deepcopy(main_info)

            # Get MM energy on QM region
            traj_ps = self.make_primary_subsys_trajectory()
            system.primary_subsys['trajectory'] = traj_ps
            system.primary_subsys_mm = mm_wrapper.compute_mm(traj_ps, include_coulomb=None)

            # Get QM energy
            self.qm_geometry = self.get_qm_positions(traj)
            charges = self.get_external_charges(system)
            self.qm_wrapper.set_external_charges(charges)
            system.qm_info = self.qm_wrapper.get_qm(self.qm_geometry)

            # Compute the total QM/MM energy based on
            # subtractive Mechanical embedding
            system.qmmm_energy = system.entire_sys['energy']\
                        - system.primary_subsys_mm['energy']\
                        + system.qm_info['energy']

            self.compute_gradients(system)
        else:
            print('only a subtractive scheme is implemented at this time')

    def electrostatic(self, system, main_info):
        """
        Gets energies of needed components and computes
        a qm/mm energy with a subtractive mechanical embedding scheme
        """
        if self.qmmm_scheme == 'subtractive':

            # Get MM energy on whole system
            system.entire_sys = deepcopy(main_info)

            # Get MM energy on QM region
            traj_ps = self.make_primary_subsys_trajectory()
            system.primary_subsys['trajectory'] = traj_ps
            system.primary_subsys_mm = mm_wrapper.compute_mm(traj_ps, include_coulomb=None)


            # Get MM coulomb energy on secondary subsystem
            traj_ss = self.make_second_subsys_trajectory()
            system.second_subsys['trajectory'] = traj_ss
            system.second_subsys_mm = mm_wrapper.compute_mm(traj_ss, include_coulomb='only')

            # Get QM energy
            self.qm_geometry = self.get_qm_positions(traj)
            charges = self.get_external_charges(system)
            self.qm_wrapper.set_external_charges(charges)
            system.qm_info = self.qm_wrapper.get_qm(self.qm_geometry)

            # Compute the total QM/MM energy based on
            # subtractive Mechanical embedding
            system.qmmm_energy = system.entire_sys['energy']\
                        - system.primary_subsys_mm['energy']\
                        + system.second_subsys_mm['energy']\
                        + system.qm_info['energy']

            self.compute_gradients(system)

        else:
            print('only a subtractive scheme is implemented at this time')

    def compute_gradients(self, system):
        # NEED TO MAKE SURE: am I working with GRADIENTS or FORCES? NEED TO MAKE SURE CONSISTENT!
        # NEED TO MAKE SURE UNITS CONSISTENT

        if scheme == 'subtractive':

            ps_mm_grad, qm_grad = system.primary_subsys['gradients'], self.qm['gradients']
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
                    for j, link in self.link_atoms.items():
                        q1 = link['qm_atom'].index
                        m1 = link['mm_atom'].index
                        link_index = link['link_atom_index']
                        g = link['scale_factor'] 
                        if atom == q1:
                            if self.boundary_treatment == 'link_atom':
                                # Project forces of link atoms onto the mm and qm atoms of the link atom bond
                                qmmm_force[q1] += -(1 - g) * ps_mm_grad[link_index] + (1 - g) * qm_grad[link_index]
                                qmmm_force[m1] += -g * ps_mm_grad[link_index] + g * qm_grad[link_index]
                                
                        # # Forces on M2 requires forces on point charges which I'm not sure about so need to double check
                        # if self.boundary_treatment == 'RC' or self.boundary_treatment == 'RCD':
                        #     qmmm_force[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                        #     qmmm_force[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]


            if system.second_subsys['gradients']:
                # iterate over list of mm atoms
                for i, atom in enumerate(self.mm_atoms):
                    qmmm_force[atom] = -1 * system.second_subsys['gradients'][i]

            system.qmmm_forces = qmmm_force
        

    def get_qm_geometry(self, qm_traj):
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
        if qm_traj is None:
            qm_traj = self.qm_trajectory

        out = ""
        line = '{:3} {: > 7.3f} {: > 7.3f} {: > 7.3f} \n '

        for i in range(qm_traj.n_atoms):
            x, y, z =   qm_traj.xyz[0][i][0],\
                        qm_traj.xyz[0][i][1],\
                        qm_traj.xyz[0][i][2]

            symbol = qm_traj.topology.atom(i).element.symbol
            
            out += line.format(symbol, x*10, y*10, z*10)

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
            self.link_atoms[i]['link_positions'] = self.positions[qm.index] + g*self.positions[mm.index]

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


    def make_primary_subsys_trajectory(self, qm_atoms=None):
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

        self.find_boundary_bonds(qm_atoms)
        traj = self.trajectory.atom_slice(qm_atoms)

        if self.qmmm_boundary_bonds():
            self.prepare_link_atoms()

            for i, link in self.link_atoms.items():

                link_element = md.element.Element.getBySymbol(link['link_atom'])

                for atom in traj.topology.atoms:
                    if atom.serial == link['qm_atom'].serial:
                        traj.topology.add_atom(name='H', element=link_element, residue=atom.residue, serial='link')

                        for atom2 in traj.topology.atoms:
                            if atom2.serial == 'link':
                                traj.topology.add_bond(atom2, atom)
                                link['link_atom_index'] = atom2.index

                traj.xyz = np.append(traj.xyz[0], [link['link_positions']], axis=0)
        

        return traj

    def make_second_subsys_trajectory(self, qm_atoms=None):
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
    
        self.mm_atoms = [i for i in range(self.traj.n_atoms) if i not in qm_atoms]

        traj = self.trajectory.atom_slice(mm_atoms)

        return traj


    def get_forces(self):

        return self.systems[self.run_ID - 1]['qmmm_forces']


'''
PUT THIS IN QMMM
'''
    def get_external_charges(self, system):
        """
        Gets the point charges of atoms from secondary subsystem for electrostatic embedding 

        Note: check to make sure positions are in angstroms

        Parameters
        ----------
        link: a bool specifying whether the link atom scheme is used for determining point charges. 
              default is False.
        RC: a bool specifying whether the RC scheme is used for determining point charges. 
              default is False.
        RCD: a bool specifying whether the RCD scheme is used for determining point charges. 
              default is False.

        Returns
        -------
        A list of charges and cooresponding positions as  xyz coordinates

        Examples
        --------
        get_external_charge(link=True)
        get_external_charge(RC=True)
        """

        if self.embedding_method == 'Mechanical':
            return None
        
        charges = []
        es = system.entire_sys

        elif self.boundary_treatment == 'link':
            for i, chrg in enumerate(es['charges']):
                # add every atom not in qm system 
                if i not in self.qm_atoms:
                    charges.append([chrg, es['positions'][i][0], es['positions'][i][1], es['positions'][i][2]])
        
        # This is for the RC and RCD schemes
        else: 
            for i, atom in self.link_atoms.items():
                mm_index = atom['mm_atom'].index
                bonds = atom['bonds_to_mm']

                # get q0
                q0 = es['charges'][mm_index] / len(bonds)

                # get positions
                positions = self.get_redistributed_positions(es['positions'], bonds, mm_index)

                if RC is True:
                    for j, chrg in enumerate(es['charges']):
                        # add every atom not in qm system or the M1 atom 
                        if j not in self.qm_atoms and j != mm_index:
                                charges.append([chrg, es['positions'][j][0], es['positions'][j][1], es['positions'][j][2]])
                    for pos in positions:
                        charges.append([q0, pos[0], pos[1], pos[2]])

                elif RCD is True:

                    q0_RCD = q0 * 2
                    for j, chrg in enumerate(es['charges']):
                        # add every atom not in qm system or the M1 atom 
                        if j not in self.qm_atoms and j != mm_index:
                            if j in bonds:
                            # modified M2 charges in RCD scheme
                                charges.append([chrg - q0, es['positions'][j][0], es['positions'][j][1], es['positions'][j][2]])
                            else:
                                charges.append([chrg, es['positions'][j][0], es['positions'][j][1], es['positions'][j][2]])
                    for pos in positions:
                        charges.append([q0_RCD, pos[0], pos[1], pos[2]])
                
        return charges

'''
PUT THE FOLLOWING IN QMMM!!
'''
    def get_redistributed_positions(self, positions, bonds, mm):
        """
        Gets the positions for the redistributed point charges in the RC and RCD schemes

        Parameters
        ----------
        positions: a list of the positions
        bonds: a list of indices of all atoms (in secondary subsystem) bonded to M1  
        mm: the index of M1

        Returns
        -------
        List of positions for the redistributed charges

        Examples
        --------
        get_redistributed_positions(positions=pos, bonds=bond, mm=mm_index)
        """
        
        pos = []
    
        for bond in bonds:
            new_pos = (positions[bond] + positions[mm]) / 2
            pos.append(new_pos)
        
        return pos
