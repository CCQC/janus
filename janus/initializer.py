import os
import json
from .psi4_wrapper import Psi4_wrapper 
from .openmm_wrapper import OpenMM_wrapper 
from .qmmm import QMMM
from .oniom_xs import ONIOM_XS
from .hot_spot import HotSpot
from .pap import PAP
from .sap import SAP

class Initializer(object):

    def __init__(self, param, as_file=True):

        if as_file is True:
            self.param = self.load_param(param)
        else:
            self.param = param

        self.qmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/qmmm.json'
        self.aqmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/aqmmm.json'
        self.psi4_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/psi4.json'
        self.openmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/openmm.json'

        self.qmmm_param = self.load_param(self.qmmm_paramfile)
        self.aqmmm_param = self.load_param(self.aqmmm_paramfile)

        try:
            print('Using pdb file ', self.param['system']['mm_pdb_file'])
        except KeyError:
            print('No pdb file given')
            
        try:
            self.qm_program = self.param['qmmm']['qm_program']
        except:
            print("No QM program was specified. Psi4 will be used")
            self.qm_program = 'Psi4'

        try:
            self.mm_program = self.param['qmmm']['mm_program']
        except:
            print("No MM program was specified. OpenMM will be used")
            self.mm_program = 'OpenMM'

        self.update_param()

    def update_param(self):

        if self.qm_program == "Psi4":
            self.qm_param = self.load_param(self.psi4_paramfile)
        else:
        # add other options for qm program here
            print("Only Psi4 currently available")

        if self.mm_program == "OpenMM":
            self.mm_param = self.load_param(self.openmm_paramfile)
        else:
        # add other options for mm program here
            print("Only OpenMM currently available")

        try:
            self.qm_param.update(self.param['qm'])
        except:
            print("No QM parameters given. Using defaults")
        try:
            self.mm_param.update(self.param['mm'])
        except:
            print("No MM parameters given. Using defaults")

        self.mm_param.update(self.param['system'])

        try:
            self.qmmm_param.update(self.param['qmmm'])
        except: 
            print("No QMMM parameters given. Using defaults")

        self.qmmm_param.update(self.param['system'])

        try:
            self.aqmmm_param.update(self.param['aqmmm'])
            self.aqmmm_param.update(self.param['qmmm'])
        except: 
            self.aqmmm_param.update(self.qmmm_param)
            print("No AQMMM parameters given. Using defaults")
        

    def initialize_wrappers(self):
        """
        Initializes the programs to use for computations
        """

        # create qm_wrapper object
        if self.qm_program == "Psi4":
            qm_wrapper = Psi4_wrapper(qm_param)
        else:
        # add other options for qm program here
            print("Only Psi4 currently available")

        # create mm_wrapper object
        if self.mm_program == "OpenMM":
            mm_wrapper = OpenMM_wrapper(mm_param)
        else:
        # add other options for mm program here
            print("Only OpenMM currently available")


        if self.qmmm_param['run_aqmmm'] is False: 
            qmmm = QMMM(qmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmm_param['aqmmm_scheme'] == 'ONIOM-XS':
            qmmm = ONIOM_XS(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmm_param['aqmmm_scheme'] == 'Hot-Spot':
            qmmm = HotSpot(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmm_param['aqmmm_scheme'] == 'PAP':
            qmmm = PAP(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmm_param['aqmmm_scheme'] == 'SAP':
            qmmm = SAP(aqmmm_param, qm_wrapper, mm_wrapper)
        else:
            print("Only ONIOM_XS, Hot Spot, PAP, and SAP currently implemented")

        # initialize mm_wrapper with information about initial system
        mm_wrapper.initialize(self.param['qmmm']['embedding_method'])
        
        return mm_wrapper, qmmm

    def load_param(self, filename):

        with open(filename) as parameter_file:
            param = json.load(parameter_file)
        
        return param
