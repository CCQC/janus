from copy import deepcopy
"""
QMMM class for QMMM computations
"""
class QMMM(object):

    def __init__(self, qm_wrapper):
        
        self.qm_wrapper = qm_wrapper
        self.qm_atoms = qm_wrapper._system.qm_atoms
        self.boundary_treatment = qm_wrapper._system.boundary_treatment
        
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
        if qm_positions == None:
            qm_positions = mm_wrapper.get_qm_positions() 
        self.qm = qm_wrapper.get_qm(qm_positions)

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
        self.entire_sys = mm_wrapper.main_info

        # Get MM energy on QM region
        self.primary_subsys = mm_wrapper.get_primary_subsys(link=True)

        # Get position and identity of link atom for QM computation if relevant
        self.boundary_info = mm_wrapper.get_boundary_info()

        # Get QM energy
        # get QM positions from pdb
        if qm_positions == None:
            qm_positions = mm_wrapper.get_qm_positions() 
        self.qm = qm_wrapper.get_qm(qm_positions)

        # Compute the total QM/MM energy based on
        # subtractive Mechanical embedding
        self.qmmm_energy = self.entire_sys['energy']\
                      - self.primary_subsys['energy']\
                      + self.qm['energy']


    def compute_gradients(self, scheme='subtractive'):
        # NEED TO MAKE SURE: am I working with GRADIENTS or FORCES? NEED TO MAKE SURE CONSISTENT!
        # NEED TO MAKE SURE UNITS CONSISTENT

        if scheme == 'subtractive':

            all_mm_grad, ps_mm_grad, qm_grad = self.entire_sys['gradients'], self.primary_subsys['gradients'], self.qm['gradients']
            qmmm_grad = np.zeros((len(all_mm_grad),3))
                
            # iterate over list of qm atoms
            for i, atom in enumerate(self.qm_atoms):

                # compute the qmmm gradient for the qm atoms: 
                # mm_entire - mm_primary - qm
                qmmm_grad[atom] += - ps_mm_grad[i] + qm_grad[i]
                
                if self.boundary_info:
                    q1 = int(self.boundary_info[0]['qm_id']) - 1
                    m1 = int(self.boundary_info[0]['mm_id']) - 1
                    g = self.boundary_info[0]['g_factor']
                    if atom == q1:
                        if self.boundary_treatment == 'link_atom':
                            # Project forces of link atoms onto the mm and qm atoms of the link atom bond
                            qmmm_grad[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                            qmmm_grad[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]
                            
                        # Forces on M2 requires forces on point charges which I'm not sure about so need to double check
                        if self.boundary_treatment == 'RC' or self.boundary_treatment == 'RCD':
                            qmmm_grad[atom] += -(1 - g) * ps_mm_grad[-1] + (1 - g) * qm_grad[-1]
                            qmmm_grad[m1] += -g * ps_mm_grad[-1] + g * qm_grad[-1]

            self.qmmm_forces = -1 * qmmm_grad
        
    def get_info(self, scheme, mm_wrapper, partition=None):

        if scheme =='subtractive':

            self.subtractive(mm_wrapper)
            self.compute_gradients(scheme='subtractive')
            
        if scheme == 'additive':
            print("Additive scheme needs some work and is not available yet") 


