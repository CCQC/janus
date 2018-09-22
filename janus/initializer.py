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
        param: a string of a filename containing a json input file
               If as_file is False, param is a dict
        as_file: a bool specifying whether param is a file or dict,
                 default is True

        Returns
        -------
        An initializer object

        Examples
        --------
        init = Initializer('input.dat')
        init = Initializer(param, as_file=False)
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
        """
        updates default parameters with user defined parameters

        Parameters
        ----------
        None

        Returns
        -------
        None

        Examples
        --------
        update_param()
        """

        if self.qm_program == "Psi4":
            self.qm_param = self.load_param(self.psi4_paramfile)
        else:
            print("Only Psi4 currently available")

        if self.mm_program == "OpenMM":
            self.mm_param = self.load_param(self.openmm_paramfile)
        else:
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
        print(self.qmmm_param)

        try:
            self.aqmmm_param.update(self.param['aqmmm'])
        except: 
            print("No AQMMM parameters given. Using defaults")

        self.aqmmm_param.update(self.qmmm_param)
        

    def initialize_wrappers(self):
        """
        Instantiates qm, mm, qmmm, and/or aqmmm wrapper objects 
        used for computation based on input parameters

        Parameters
        ----------
        None

        Returns
        -------
        mm_wrapper object, qmmm object

        Examples
        --------
        mm_wrapper, qmmm = initialize_wrapers()
        """

        # create qm_wrapper object
        if self.qm_program == "Psi4":
            qm_wrapper = Psi4_wrapper(self.qm_param)
        else:
        # add other options for qm program here
            print("Only Psi4 currently available")

        # create mm_wrapper object
        if self.mm_program == "OpenMM":
            mm_wrapper = OpenMM_wrapper(self.mm_param)
        else:
        # add other options for mm program here
            print("Only OpenMM currently available")


        if self.qmmm_param['run_aqmmm'] is False: 
            qmmm = QMMM(self.qmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'ONIOM-XS':
            qmmm = ONIOM_XS(self.aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'Hot-Spot':
            qmmm = HotSpot(self.aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'PAP':
            qmmm = PAP(self.aqmmm_param, qm_wrapper, mm_wrapper)
        elif self.aqmmm_param['aqmmm_scheme'] == 'SAP':
            qmmm = SAP(self.aqmmm_param, qm_wrapper, mm_wrapper)
        else:
            print("Only ONIOM_XS, Hot Spot, PAP, and SAP currently implemented")

        # initialize mm_wrapper with information about initial system
        mm_wrapper.initialize(self.qmmm_param['embedding_method'])
        
        return mm_wrapper, qmmm

    def load_param(self, filename):
        """
        Converts a json file into a dictionary
    
        Parameters
        ----------
        filename: string of json file

        Returns
        -------
        a dictionary

        Examples
        --------
        dict = load_param('input.json')
        """

        with open(filename) as parameter_file:
            param = json.load(parameter_file)
        
        return param

