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

    def _init__(self, paramfile):

        self.param = self.load_param(paramfile)

        self.qmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/qmmm.json'
        self.aqmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/aqmmm.json'
        self.psi4_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/psi4.json'
        self.openmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/openmm.json'

        self.qmmm_param = self.load_param(self.qmmm_paramfile)
        self.aqmmm_param = self.load_param(self.aqmmm_paramfile)

        self.qm_program = self.param['qmmm']['qm_program']
        self.mm_program = self.param['qmmm']['mm_program']

        if 'aqmmm' in self.param:
            self.aqmmm_scheme = self.param['aqmmm']['aqmmm_scheme']
        else:
            self.aqmmm_scheme = None
            

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

        self.qm_param.update(self.param['qm'])
        self.mm_param.update(self.param['mm'])
        self.mm_param.update(self.param['system'])
        self.qmmm_param.update(self.param['qmmm'])
        self.qmmm_param.update(self.param['system'])
        self.aqmmm_param.update(self.param['aqmmm'])
        self.aqmmm_param.update(self.param['qmmm'])
        

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


        if self.aqmmm_scheme is None:
            qmmm = QMMM(qmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_scheme == 'ONIOM-XS':
            qmmm = ONIOM_XS(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_scheme == 'Hot-Spot':
            qmmm = HotSpot(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_scheme == 'PAP':
            qmmm = PAP(aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_scheme == 'SAP':
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
