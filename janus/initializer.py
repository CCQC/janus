import json
import os
from janus.qm_wrapper import Psi4Wrapper 
from janus.mm_wrapper import OpenMMWrapper 
from janus.qmmm import QMMM, OniomXS, HotSpot, PAP, SAP

class Initializer(object):
    """
    Class that initializes the initial parameters 
    and necessary wrapppers
    """

    def __init__(self, parameters, as_file=True):
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

        self.system = None
        self.qmmm = None
        self.aqmmm = None
        self.md = None
        self.ll = None
        self.hl = None

        self.system_info_format = 'pdb'
        self.run_qmmm = True
        self.run_aqmmm = False
        self.run_md = False
        self.aqmmm_scheme = None
        self.md_restart = False
        self.ll_program = 'OpenMM'
        self.hl_program = 'Psi4'
        self.md_simulation_program = "OpenMM"

        if as_file is True:
            self.param = self.load_param(param)
        else:
            self.param = param

        self.set_attributes(self.param, self)
        self.aqmmm['qmmm_param'] = self.qmmm

        if self.system is not None and 'system_info' in self.system:
            self.set_attributes(self.system, self)
        else:
            raise KeyError("system_info needs to be specified")

        self.get_wrappers()
    
    def get_wrappers(self):
        """
        updates default parameters with user defined parameters

        """

        if self.hl_program == "Psi4":
            self.hl_wrapper = Psi4Wrapper
        elif self.hl_program == "OpenMM":
            self.hl_wrapper = OpenMMWrapper
        else:
            raise ValueError("Only Psi4 and OpenMM currently available to be used in high level computations")

        if self.ll_program == "OpenMM":
            self.ll_wrapper = OpenMMWrapper
        else:
            raise ValueError("Only OpenMM currently available to be used in low level computations")

        if self.run_md is True:

            if self.md_simulation_program == "OpenMM":
                self.md_sim_wrapper = OpenMMWrapper
            else:
                raise ValueError("Only OpenMM currently available")

        #self.qmmm_param.update(self.system_param)
        #self.aqmmm_param.update(self.system_param)
        #self.md_sim_param.update(self.system_param)


    def initialize_wrappers(self):
        """
        Instantiates qm, mm, qmmm, and/or aqmmm wrapper objects 
        used for computation based on input parameters


        Returns
        -------
        MMWrapper object
        QMMM wrapper object

        """

        "TODO: figure out input for software wrappers"
        # create hl_wrapper object
        hl_wrapper = self.hl_wrapper(self.hl_param)
        # create ll_wrapper object
        ll_wrapper = self.ll_wrapper(self.ll_param)
        
        if self.aqmmm_scheme is None:
            qmmm = QMMM(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, **self.qmmm)
        elif self.aqmmm_scheme == 'ONIOM-XS':
            qmmm = OniomXS(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, aqmmm_param=self.aqmmm)
        elif self.aqmmm_scheme == 'Hot-Spot':
            qmmm = HotSpot(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, aqmmm_param=self.aqmmm)
        elif self.aqmmm_scheme == 'PAP':
            qmmm = PAP(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, aqmmm_param=self.aqmmm)
        elif self.aqmmm_scheme == 'SAP':
            qmmm = SAP(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, aqmmm_param=self.aqmmm)
        elif self.aqmmm_scheme == 'DAS':
            qmmm = DAS(hl_wrapper, ll_wrapper, self.system_info, self.system_info_format, aqmmm_param=self.aqmmm)
        else:
            raise ValueError("{} not recognized as a currently implemented method".format(self.aqmmm_param['aqmmm_scheme']))

        if self.run_md is True:

            md_sim_wrapper = self.md_sim_wrapper(self.md_sim_param)

            if self.md_restart is True:
                # initialize mm_wrapper with information about initial system
                md_sim_wrapper.initialize(self.qmmm_param['embedding_method'])
            elif self.md_restart is False:
                md_sim_wrapper.restart(self.qmmm_param['embedding_method'])

            return md_sim_wrapper, qmmm
 
        else:
            return ll_wrapper, qmmm
        

    def load_param(self, fname):
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

        if fname.endswith(".json"):
            rfunc = json.load
        else:
            raise TypeError("Did not understand file type {}.".format(fname))

        with open(fname, 'r') as handle:
            ret = rfunc(handle)

        return ret


    def set_attributes(self, dictionary, obj):

        for k, v in dictionary.items():
            setattr(obj, k, v)


