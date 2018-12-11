import numpy as np
from copy import deepcopy
import mdtraj as md
from janus.partition import Partition

class HystereticPartition(Partition):

    def __init__(self, trajectory, topology, Rmin_qm, Rmax_qm, Rmin_bf, Rmax_bf):

        self.Rmin_qm = Rmin_qm  
        self.Rmax_qm = Rmax_qm
        self.Rmin_bf = Rmin_bf 
        self.Rmax_bf = Rmax_bf

        super().__init__(trajectory, topology, 'hysteretic')

    def define_buffer_zone(self, qm_center, prev_qm, prev_bf):
        """
        Determines buffer group atoms.
        Gets the buffer groups in the buffer zone based on a distance 
        partitioning scheme, saves each buffer group as a 
        :class:`~janus.system.Buffer` object,
        and saves all buffer groups in the dictionary self.buffer_groups.
        For water as a solvent, considers the whole water molecule as a buffer group.

        Note
        ----
        Currently only worked with explicit solvent based systems. 
        Cannot treat buffer atoms that are part of large molecular structures (e.g., proteins).

        Parameters
        ----------
        qm_center : list 
            the indicies that make up the qm center

        """

        self.qm_residues = []
        self.buffer_groups = {}
        residue_tracker = [] 
        top = self.topology

        self.find_buffer_atoms(qm_center)

        for i in self.buffer_atoms:
            idx = top.atom(i).residue.index

            if idx not in residue_tracker:
                residue_tracker.append(idx)
                buf = self.get_residue_info(idx)

                # all within Rmin_qm considered QM
                if buf.r_i < self.Rmin_qm:
                    self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, add=True)

                # Between Rmin_qm and Rmax_qm, only considered QM if previously QM
                if (buf.r_i >= self.Rmin_qm and buf.r_i < self.Rmax_qm):

                    if (buf.ID in prev_qm and buf.ID not in prev_bf):
                        self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, add=True)
                    elif (buf.ID not in prev_qm and buf.ID in prev_bf):
                        self.buffer_groups[idx] = buf
                        self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)
                    else:
                        raise Exception('Inconsistent group definition. Group was both QM and BZ atom')

                # Between Rmax_qm and Rmin_bf considered buffer
                if (buf.r_i >= self.Rmax_qm and buf.r_i < self.Rmin_bf):
                    self.buffer_groups[idx] = buf
                    self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)

                # Between Rmin_bf and Rmax_bf only considered as BZ atom if previously BZ atom
                if (buf.r_i >= self.Rmin_bf and buf.r_i < self.Rmax_bf):
                    if (buf.ID not in prev_qm and buf.ID not in prev_bf):
                        self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)
                    elif (buf.ID not in prev_qm and buf.ID in prev_bf):
                        self.buffer_groups[idx] = buf
                        self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)
                    else:
                        raise Exception('Inconsistent group definition.')

                # Beyond Rmax_bf not buffer atom
                if buf.r_i >= self.Rmax_bf:
                    self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)
                    
        qm_atoms = deepcopy(self.qm_atoms)
        # tracking qm_residues and cleaning up qm
        for i in qm_atoms:
            idx = top.atom(i).residue.index
            if idx not in self.qm_residues:
                res = self.get_residue_info(idx)

                # Cleaning up additional qm residues not catched previously
                if res.r_i >= self.Rmax:
                    self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, remove=True)
                elif res.r_i < self.Rmin:
                    self.edit_atoms(atoms=self.qm_atoms, res_idx=idx, add=True)
                    self.qm_residues.append(idx)

    def find_buffer_atoms(self, qm_center):
        """
        Find the buffer groups whose COM falls in between Rmin and Rmax

        Parameters
        ----------
        qm_center : list 
            the indicies that make up the qm center

        """

        temp_traj, qm_center_idx = self.compute_qm_center_info(qm_center)

        rmin_atoms = md.compute_neighbors(temp_traj, self.Rmin_qm/10, qm_center_idx)
        rmax_atoms = md.compute_neighbors(temp_traj, self.Rmax_bf/10, qm_center_idx)
        self.buffer_atoms = np.setdiff1d(rmax_atoms, rmin_atoms)
        print('buffer atoms identified by find_buffer_atom function:')
        print(self.buffer_atoms)
        self.qm_atoms = rmin_atoms[0].tolist()

        if self.COM_as_qm_center is False:
            self.qm_atoms.append(qm_center[0])

        print('qm_atoms identified by the find_buffer_atom function: ' )
        print(self.qm_atoms)

    def get_Rmin_qm(self):
        return self.Rmin_qm

    def get_Rmin_bf(self):
        return self.Rmin_bf

    def get_Rmax_qm(self):
        return self.Rmax_qm

    def get_Rmax_bf(self):
        return self.Rmax_bf

    def set_Rmin_qm(self, Rmin):
        self.Rmin_qm = Rmin

    def set_Rmin_bf(self, Rmin):
        self.Rmin_bf = Rmin

    def set_Rmax_qm(self, Rmax):
        self.Rmax_qm = Rmax

    def set_Rmax_bf(self):
        self.Rmax_bf = Rmax

