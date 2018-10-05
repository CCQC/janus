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
    """
    Class that initializes the initial parameters 
    and necessary wrapppers
    """

    def __init__(self, param, as_file=True):
        """
        Initializes the Initializer class with parameters and updates 
        defaults with user defined parameters

        Parameters
        ----------
        param : str 
            filename containing a json input file, if as_file is False, param is a dict
        as_file : bool
             param is a file if True, dict if False

        """

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
            self.steps = self.param['system']['md_steps']
        except:
            print("Number of steps is not specified")

        try:
            self.md_sim_prog = self.param['system']['md_simulation_program']
        except:
            print("No program for an MD simulation is specified. Not performing a simulation")
            self.md_sim_prog = None
            
        try:
            self.hl_program = self.param['qmmm']['hl_program']
        except:
            print("No QM program was specified. Psi4 will be used")
            self.hl_program = 'Psi4'

        try:
            self.mm_program = self.param['qmmm']['ll_program']
        except:
            print("No MM program was specified. OpenMM will be used")
            self.ll_program = 'OpenMM'


        self.update_param()

    def update_param(self):
        """
        updates default parameters with user defined parameters

        """

        if self.hl_program == "Psi4":
            self.hl_param = self.load_param(self.psi4_paramfile)
        elif self.hl_program == "OpenMM"
            self.hl_param = self.load_param(self.openmm_paramfile)
        else:
            print("Only Psi4 and OpenMM currently available")

        if self.ll_program == "OpenMM":
            self.ll_param = self.load_param(self.openmm_paramfile)
        else:
            print("Only OpenMM currently available")

        self.ll_param.update(self.param['system']

        if self.md_sim_prog == "OpenMM":
            self.md_sim_param = self.load_param(self.openmm_paramfile)
        else:
            print("Only OpenMM currently available")

        self.md_sim_param.update(self.param['system']

        try:
            self.hl_param.update(self.param['hl'])
        except:
            print("No QM parameters given. Using defaults")
        try:
            self.ll_param.update(self.param['ll'])
        except:
            print("No MM parameters given. Using defaults")

        try:
            self.qmmm_param.update(self.param['qmmm'])
        except: 
            print("No QMMM parameters given. Using defaults")

        self.qmmm_param.update(self.param['system'])

        try:
            self.aqmmm_param.update(self.param['aqmmm'])
        except: 
            print("No AQMMM parameters given. Using defaults")

        self.aqmmm_param.update(self.qmmm_param)
        

    def initialize_wrappers(self, simulation=True):
        """
        Instantiates qm, mm, qmmm, and/or aqmmm wrapper objects 
        used for computation based on input parameters


        Returns
        -------
        MM_wrapper object
        QMMM_wrapper object

        """

        # create qm_wrapper object
        if self.hl_program == "Psi4":
            hl_wrapper = Psi4_wrapper(self.qm_param)
        if self.hl_program == "OpenMM":
            hl_wrapper = OpenMM_wrapper(self.qm_param)
        else:
        # add other options for qm program here
            print("Only Psi4 and OpenMM currently available")

        # create mm_wrapper object
        if self.mm_program == "OpenMM":
            ll_wrapper = OpenMM_wrapper(self.mm_param)
        else:
        # add other options for mm program here
            print("Only OpenMM currently available")

        if self.md_sim_prog == "OpenMM":
            md_sim_wrapper = OpenMM_wrapper(self.md_sim_param)
        else:
            print("Only OpenMM currently available")


        if self.qmmm_param['run_aqmmm'] is False: 
            qmmm = QMMM(self.qmmm_param, hl_wrapper, ll_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'ONIOM-XS':
            qmmm = ONIOM_XS(self.aqmmm_param, hl_wrapper, ll_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'Hot-Spot':
            qmmm = HotSpot(self.aqmmm_param, hl_wrapper, ll_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'PAP':
            qmmm = PAP(self.aqmmm_param, hl_wrapper, ll_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'SAP':
            qmmm = SAP(self.aqmmm_param, hl_wrapper, ll_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'DAS':
            qmmm = SAP(self.aqmmm_param, hl_wrapper, ll_wrapper)
        else:
            print("Only ONIOM_XS, Hot Spot, PAP, SAP, and DAS currently implemented")

        if simulation is True:
            # initialize mm_wrapper with information about initial system
            md_sim_wrapper.initialize(self.qmmm_param['embedding_method'])
            return md_sim_wrapper, qmmm

        else:
            return mm_wrapper, qmmm
        

    def load_param(self, filename):
        """
        Converts a json file into a dictionary
    
        Parameters
        ----------
        filename : str
            name of json file

        Returns
        -------
        dict
            parameters contained in filename

        """

        with open(filename) as parameter_file:
            param = json.load(parameter_file)
        
        return param

