from copy import deepcopy
import numpy as np
import mdtraj as md
from .system import System
"""
QMMM class for QMMM computations
"""
class QMMM(object):

    def __init__(self, param, hl_wrapper, ll_wrapper, md_simulation_program=None):
        """
        Initializes QMMM class with parameters given in param

        Parameters
        ----------
        param : dict
            contains parameters for QM/MM and adaptive QM/MM computations
            Individual parameters include:
            - qmmm_scheme: str of scheme for computing QM/MM energies and gradients, 
                            only substractive available(default)
            - embedding_method: str of embedding method to use for QM/MM. 
                                Mechanical(default) and Electrostatic available
            - boundary_treatment: str of method for treating dangling bonds in the QM region,
                                    link_atom(default), RC, and RCD available
            - link_atom_element: str of element to use for link atom,
                                    default is H. Beware of using others (not all functionality tested)
            - qm_atoms: list of indices that define the qm_region. This is not static for adaptive QM/MM computations
            - run_aqmmm: bool to specify whether adaptive QM/MM wrappers are called,
                            default is True. Regular QM/MM run if False

        hl_wrapper : a qm_wrapper or mm_wrapper object, depending on the user input
        ll_wrapper : a mm_wrapper object only for now
        md_simulation_program : str
            The program that performs the MD time step integration

        """
        
        self.class_type = 'QMMM'
        self.hl_wrapper = hl_wrapper
        self.ll_wrapper = ll_wrapper
        self.md_simulation_program = md_simulation_program
        self.qm_geometry = None
        self.run_ID = 0

        self.traj = self.convert_input(param['system_info'], param['system_info_format'])
        self.topology = self.traj.topology
        self.positions = self.traj.xyz[0]

        self.qm_atoms = param['qm_atoms']
        self.qmmm_scheme = param['qmmm_scheme']
        self.embedding_method = param['embedding_method']
        self.boundary_treatment = param['boundary_treatment']
        self.link_atom_element = param['link_atom_element']

        self.systems = {}

    def run_qmmm(self, main_info):
        """
        Updates the positions and topology given in main_info,
        and determines the QM/MM energy and gradients

        Parameters
        ----------
        main_info : dict 
            contains the energy and forces for the whole system

        """

        self.update_traj(main_info['positions'], main_info['topology'])

        system = System(qm_indices=self.qm_atoms, qm_residues=None, run_ID=self.run_ID)

        if self.embedding_method =='Mechanical':
            self.mechanical(system, main_info)
        elif self.embedding_method =='Electrostatic':
            self.electrostatic(system, main_info)
        else:
            print('only mechanical and electrostatic embedding schemes implemented at this time')
            
        self.systems[self.run_ID] = {}
        self.systems[self.run_ID][system.partition_ID] = system
        self.systems[self.run_ID]['qmmm_forces'] = system.qmmm_forces 
        self.systems[self.run_ID]['qmmm_energy'] = system.qmmm_energy 
        self.systems[self.run_ID]['kinetic_energy'] = main_info['kinetic']

        if self.run_ID % 10 == 0:
            print(self.run_ID, self.systems[self.run_ID]['qmmm_energy'])
        # updates current step count
        self.run_ID += 1

    def update_traj(self, position, topology):
        """
        Updates the positions and topology of self.traj,
        a MDtraj trajectory object

        Parameters
        ----------
        positions : list 
            positions in nm
        topology : OpenMM topology object
            If the mm program is  OpenMM

       """ 
        # later can think about saving instead of making new instance
        # convert openmm topology to mdtraj topology
        if self.md_simulation_program == 'OpenMM':
            top = md.Topology.from_openmm(topology)
        for atom in top.atoms:
            atom.serial = atom.index + 1

        self.traj = md.Trajectory(position, top)
        self.topology = self.traj.topology
        self.positions = self.traj.xyz[0]


    def mechanical(self, system, main_info):
        """
        Gets energies of needed components and computes
        a QM/MM energy with a subtractive mechanical embedding scheme
        using the formula

        E(QM/MM) = E(MM)_entire_sys - E(MM)_primary_subsys + E(QM)_primary_subsys

        Parameters
        ----------
        system : System object
        main_info : dict 
            contains the energy and forces for the whole system

        """

        if self.qmmm_scheme == 'subtractive':
            # Get MM energy on whole system
            system.entire_sys = deepcopy(main_info)

            #print(system.entire_sys['energy'])
            # Get MM energy on QM region
            traj_ps, link_indices = self.make_primary_subsys_trajectory(qm_atoms=system.qm_atoms)
            system.primary_subsys['trajectory'] = traj_ps
            system.primary_subsys['ll'] = self.ll_wrapper.get_energy_and_gradient(traj_ps, include_coulomb='no_link', link_atoms=link_indices)

            # Get QM energy
            system.primary_subsys['hl'] = self.hl_wrapper.get_energy_and_gradient(traj_ps)

            # Compute the total QM/MM energy based on
            # subtractive Mechanical embedding
            system.qmmm_energy = system.entire_sys['energy']\
                        - system.primary_subsys['ll']['energy']\
                        + system.primary_subsys['hl']['energy']

            self.compute_gradients(system)
        else:
            print('only a subtractive scheme is implemented at this time')

    def electrostatic(self, system, main_info):
        """
        Gets energies of needed components and computes
        a QM/MM energy with a subtractive electrostatic embedding scheme

        E(QM/MM) = E(MM no coulomb)_entire_sys - E(MM no coulomb)_primary_subsys 
                 + E(QM)_primary_subsys + E(MM just coulomb)_secondary_subsys

        Parameters
        ----------
        system : System object
        main_info : dict 
            contains the energy and forces for the whole system

        """ 

        if self.qmmm_scheme == 'subtractive':

            # Get MM energy on whole system
            system.entire_sys = deepcopy(main_info)
            #print(system.entire_sys['energy'])

            # Get MM energy on QM region
            traj_ps, link_indices = self.make_primary_subsys_trajectory(qm_atoms=system.qm_atoms)
            system.primary_subsys['trajectory'] = traj_ps
            system.primary_subsys['ll'] = self.ll_wrapper.get_energy_and_gradient(traj_ps, include_coulomb=None)

            # Get MM coulomb energy on secondary subsystem
            traj_ss = self.make_second_subsys_trajectory()
            system.second_subsys['trajectory'] = traj_ss
            system.second_subsys['ll'] = self.ll_wrapper.get_energy_and_gradient(traj_ss, include_coulomb='only')

            # Get QM energy
            charges = self.get_external_charges(system)
            system.primary_subsys['hl'] = self.hl_wrapper.get_energy_and_gradient(traj_ps, charges=charges)

            # Compute the total QM/MM energy based on
            # subtractive Mechanical embedding
            system.qmmm_energy = system.entire_sys['energy']\
                        - system.primary_subsys['ll']['energy']\
                        + system.second_subsys['ll']['energy']\
                        + system.primary_subsys['hl']['energy']

            self.compute_gradients(system)

        else:
            print('only a subtractive scheme is implemented at this time')

    def compute_gradients(self, system):
        """
        Computes the QM/MM gradients 
        TODO:RCD gradients need work

        Note
        ----
        RCD gradients currently not implemented

        Parameters
        ----------
        system : System object 
            Contains the gradients at  high level and low level  

        """
        # NEED TO MAKE SURE: am I working with GRADIENTS or FORCES? NEED TO MAKE SURE CONSISTENT!
        # NEED TO MAKE SURE UNITS CONSISTENT

        if self.qmmm_scheme == 'subtractive':

            ps_mm_grad, qm_grad = system.primary_subsys['ll']['gradients'], system.primary_subsys['hl']['gradients']
            print('ps_mm', ps_mm_grad)
            print('qm', qm_grad)
            qmmm_force = {}
           # print('sys qm_atom', system.qm_atoms)
                
            # iterate over list of qm atoms
            for i, atom in enumerate(system.qm_atoms):

                # compute the qmmm gradient for the qm atoms: 
                # mm_entire - mm_primary - qm
                qmmm_force[atom] = np.zeros(3)
                # these are in units of au_bohr, convert to openmm units in openmm wrapper
                qmmm_force[atom] += -1 * (- ps_mm_grad[i] + qm_grad[i])
                
                # treating gradients for link atoms
                if self.qmmm_boundary_bonds:
                    for j, link in self.link_atoms.items():
                        if isinstance(j, int):
                            q1 = link['qm_atom'].index
                            m1 = link['mm_atom'].index
                            link_index = link['link_atom_index']
                            g = link['scale_factor'] 
                            if atom == q1:
                                if self.boundary_treatment == 'link_atom':
                                    # Project forces of link atoms onto the mm and qm atoms of the link atom bond
                                    # need to make sure sign is correct
                                    qmmm_force[q1] += -(1 - g) * ps_mm_grad[link_index] + (1 - g) * qm_grad[link_index]
                                    qmmm_force[m1] = -g * ps_mm_grad[link_index] + g * qm_grad[link_index]
                                    
                            # # Forces on M2 requires forces on point charges which I'm not sure about so need to double check
                            # if self.boundary_treatment == 'RC' or self.boundary_treatment == 'RCD':
                            #     qmmm_force[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                            #     qmmm_force[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]


            if 'll' in system.second_subsys:
                # iterate over list of mm atoms
                for i, atom in enumerate(self.mm_atoms):
                    qmmm_force[atom] = -1 * system.second_subsys['ll']['gradients'][i]

            system.qmmm_forces = qmmm_force
        

    def find_boundary_bonds(self, qm_atoms=None):
        """
        Identified any covalent bonds that the QM/MM boundary cuts across

        Parameters
        ----------
        qm_atoms : list 
            atom indicies corresponding to the atoms in
            the primary subsystem. Default is None and uses self.qm_atoms

        Examples
        --------
        find_boundary_bonds()
        find_boundary_bonds(qm_atoms=[0,1,2,3])
        """

        if qm_atoms is None:
            qm_atoms = self.qm_atoms

        qm_atoms = self.edit_qm_atoms(qm_atoms)

        self.qmmm_boundary_bonds = []
        # determining if there are bonds that need to be cut
        for bond in self.topology.bonds:
            # find any bonds that involve the qm atoms
            # iterates over tuple of mdtraj atom objects
            if (bond[0].index in qm_atoms or bond[1].index in qm_atoms):
                # isolate bonds that involve one in the qm atoms and one outside
                if (bond[0].index not in qm_atoms or bond[1].index not in qm_atoms):
                    qm_atom = {}
                    mm_atom = {}
                    if bond[0].index in qm_atoms:
                        qm_atom = bond[0]
                        mm_atom = bond[1]
                    else:
                        qm_atom = bond[1]
                        mm_atom = bond[0]
                    self.qmmm_boundary_bonds.append((qm_atom, mm_atom))

    def edit_qm_atoms(self, qm_atoms=None, solvent='water'):
        """
        Cleans up the qm_atoms. If the definition of the primary subsystem 
        cuts across a solvent molecule, will move the whole molecule in the qm_atoms 
        if COM of solvent molecule in qm_atoms already, or delete parts of the solvent
        molecule from qm_atoms if COM of solvent molecule not in qm_atoms already

        Parameters
        ----------
        qm_atoms : list 
            atom indicies corresponding to the atoms in
            the primary subsystem. Default is None and uses self.qm_atoms

        solvent : str 
            identifies what solvent needs to be edited.
            Only water supported for now(default)

        Returns
        -------
        list
            edited qm_atoms

        Examples
        --------
        atoms = edit_qm_atoms()
        atoms = edit_qm_atoms(qm_atoms=[0,1,2])
        """

        if qm_atoms is None:
            qm_atoms = self.qm_atoms
        
        top = self.topology
        residue_tracker = [] 
        self.qm_residues = [] 
        qm_atoms_copy = deepcopy(qm_atoms)
        for i in qm_atoms_copy:
            idx = top.atom(i).residue.index
            # make sure just go through each residues once
            if idx not in residue_tracker:
                residue_tracker.append(idx)
                self.qm_residues.append(idx)
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
                                self.qm_residues.remove(idx)
                                for a1 in res.atoms:
                                    if (a1.element.symbol =='H' and a1.index in qm_atoms):
                                        qm_atoms.remove(a1.index)

        qm_atoms.sort()
        return qm_atoms

    def prepare_link_atom(self):
        """
        Saves the qm and mm atom associated with the bond being cut across the QM/MM boundary
        and computes the g scaling factor for the link atom.
        """

        self.link_atoms = {}
        self.link_atoms['all_mm'] = []
        self.link_atoms['all_outer_bonds'] = []

        for i, bond in enumerate(self.qmmm_boundary_bonds):

            qm = bond[0]
            mm = bond[1]

            self.link_atoms[i] = {}

            self.link_atoms[i]['qm_atom'] = qm
            self.link_atoms[i]['mm_atom'] = mm
            self.link_atoms['all_mm'].append(mm.index)

            self.link_atoms[i]['link_atom'] = self.link_atom_element
            g = System.compute_scale_factor_g(qm.element.symbol, mm.element.symbol, self.link_atom_element)
            self.link_atoms[i]['scale_factor'] = g 
            # this is in nm
            self.link_atoms[i]['link_positions'] = (1-g) * self.positions[qm.index] + g*self.positions[mm.index]

            if self.boundary_treatment == 'RC' or self.boundary_treatment == 'RCD':
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
                self.link_atoms['all_outer_bonds'].append(bonds)


    def make_primary_subsys_trajectory(self, qm_atoms=None):
        '''
        Creates a MDtraj trajectory object with just the 
        primary subsystem, and adds in any link atoms

        Note
        ----
        Currently adds just H as link atom, need to expand. Also, H is added
        as a very specific H1 atom, or else the connectivity in the create_new_residue_templates
        function in openmm gets messed up and gives a "set of atoms match but bonds are different
        error"

        Parameters
        ----------
        qm_atoms : list 
            atom indicies corresponding to the atoms in
            the primary subsystem. Default is None and uses self.qm_atoms

        Returns
        -------
        MDtraj trajectory object
        list
            The link atom indices in traj

        Examples
        --------
        make_primary_subsys_trajectory([0,1,2])
        make_primary_subsys_trajectory()
        '''

        if qm_atoms is None:
            qm_atoms = self.qm_atoms
        #print(qm_atoms)

        self.find_boundary_bonds(qm_atoms)
        traj = self.traj.atom_slice(qm_atoms)

        link_indices = []
        if self.qmmm_boundary_bonds:
            self.prepare_link_atom()

            for i, link in self.link_atoms.items():
                if isinstance(i, int):
                
                    link_element = md.element.Element.getBySymbol(link['link_atom'])

                    for atom in traj.topology.atoms:
                        if atom.serial == link['qm_atom'].serial:
                            traj.topology.add_atom(name='H1', element=link_element, residue=atom.residue, serial='link')
                            for atom2 in traj.topology.atoms:
                                if atom2.serial == 'link':
                                    traj.topology.add_bond(atom2, atom)
                                    link['link_atom_index'] = atom2.index

                    link_indices.append(link['link_atom_index'])
                    traj.xyz = np.append(traj.xyz[0], [link['link_positions']], axis=0)
        
        return traj, link_indices

    def make_second_subsys_trajectory(self, qm_atoms=None):
        '''
        Creates a MDtraj trajectory object with just the 
        secondary subsystem

        Parameters
        ----------
        qm_atoms : list 
            atom indicies corresponding to the atoms in
            the primary subsystem. Default is None and uses self.qm_atoms

        Returns
        -------
        MDtraj trajectory object

        Examples
        --------
        make_second_subsys_trajectory([0,1,2])
        make_second_subsys_trajectory()
        '''

        if qm_atoms is None:
            qm_atoms = self.qm_atoms
    
        self.mm_atoms = [i for i in range(self.traj.n_atoms) if i not in qm_atoms]

        traj = self.traj.atom_slice(self.mm_atoms)

        return traj


    def get_forces(self, run_ID=None):
        """
        Function to return qmmm forces

        Parameters
        ----------
        run_ID : int
            identifies which step to get forces from
        
        Returns
        -------
        dict
            qmmm forces in au/bohr

        Examples
        --------
        forces = get_forces()
        """
        if run_ID is None:
            run_ID = self.run_ID - 1

        return self.systems[run_ID]['qmmm_forces']


    def get_external_charges(self, system):
        #TODO: at some point maybe need to migrate this to mm_wrapper
        """
        Gets the point charges of atoms from secondary subsystem for electrostatic embedding 

        Parameters
        ----------
        system : System object 

        Returns
        -------
        list
            charges and corresponding positions in angstroms as xyz coordinates

        """
        charges = []
        # in angstroms
        es_pos = 10*system.entire_sys['positions']
        charge = self.ll_wrapper.get_main_charges()

        if self.embedding_method == 'Mechanical':
            return None

        elif self.boundary_treatment == 'link_atom':
            for i, chrg in enumerate(charge):
                # add every atom not in qm system 
                if i not in system.qm_atoms:
                    # save positions in angstroms
                    charges.append([chrg, es_pos[i][0], es_pos[i][1], es_pos[i][2]])
        
        # This is for the RC and RCD schemes
        elif self.boundary_treatment == 'RC':

            for i, chrg in enumerate(charge):
                # add every atom not in qm system or the M1 atom 
                if i not in system.qm_atoms and i not in self.link_atoms['all_mm']:
                        charges.append([chrg, es_pos[i][0], es_pos[i][1], es_pos[i][2]])
            
            bonds = self.link_atoms['all_outer_bonds']

            for i, index in enumerate(self.link_atoms['all_mm']):

                # check to see that the M1 atom is attached to any M2 atoms
                if bonds[i]:
                    # get q0
                    q0 = charge[index] / len(bonds[i])

                    # get positions in angstroms
                    positions = self.get_redistributed_positions(es_pos, bonds[i], index)

                    for pos in positions:
                        charges.append([q0, pos[0], pos[1], pos[2]])

        elif self.boundary_treatment == 'RCD':

            bonds = self.link_atoms['all_outer_bonds']
            m1_m2 = []

            for i, index in enumerate(self.link_atoms['all_mm']):

                m1_m2.append(index)
                # get q0

                if bonds[i]:

                    q0 = charge[index] / len(bonds[i])
                    q0_RCD = q0 * 2

                    # get positions in angstroms
                    positions = self.get_redistributed_positions(es_pos, bonds[i], index)

                    for pos in positions:
                        charges.append([q0_RCD, pos[0], pos[1], pos[2]])

                    for bond in bonds[i]:
                        m1_m2.append(bond)
                        charges.append([charge[bond] - q0, es_pos[bond][0], es_pos[bond][1], es_pos[bond][2]])

            for i, chrg in enumerate(charge):
                # add every atom not in qm system or the M1 atom 
                if i not in system.qm_atoms and i not in m1_m2:
                        charges.append([chrg, es_pos[i][0], es_pos[i][1], es_pos[i][2]])

                
        return charges

    def get_redistributed_positions(self, positions, bonds, mm):
        """
        Gets the positions for the redistributed point charges in the RC and RCD schemes

        Parameters
        ----------
        positions : list 
        bonds : list 
            indices of all atoms (in secondary subsystem) bonded to M1  
        mm : int 
            the index of M1

        Returns
        -------
        list
            positions for the redistributed charges

        """
        
        pos = []
    
        for bond in bonds:
            new_pos = ((positions[bond] + positions[mm]) / 2)
            pos.append(new_pos)
        
        return pos

            
    def convert_input(self, fil, form):

        if form == 'pdb':
            traj = md.load(fil)

        return traj
