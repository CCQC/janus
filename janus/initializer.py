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

        self.qmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/qmmm.json'
        self.aqmmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/aqmmm.json'
        self.md_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/md.json'
        self.system_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/system.json'
        self.psi4_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/psi4.json'
        self.openmm_paramfile = os.path.abspath(os.path.dirname(__file__)) + '/default_input/openmm.json'

        self.system_param = self.load_param(self.system_paramfile)
        self.qmmm_param = self.load_param(self.qmmm_paramfile)
        self.aqmmm_param = self.load_param(self.aqmmm_paramfile)
        self.md_sim_param = self.load_param(self.md_paramfile)

        if as_file is True:
            self.param = self.load_param(param)
        else:
            self.param = param

        self.update_param()
    
    def update_param(self):
        """
        updates default parameters with user defined parameters

        """

        if ('system' in self.param and 'system_info' in self.param['system']):
            self.system_param.update(self.param['system'])
        else:
            raise KeyError("system_info needs to be specified")
            
        if 'qmmm' in self.param:
            self.qmmm_param.update(self.param['qmmm'])

        if 'aqmmm' in self.param:
            self.aqmmm_param.update(self.param['aqmmm'])
            self.aqmmm_param.update(self.qmmm_param)
            #print("No AQMMM parameters given, running traditional QM/MM")

        if self.qmmm_param['hl_program'] == "Psi4":
            self.hl_param = self.load_param(self.psi4_paramfile)
            self.hl_wrapper = Psi4Wrapper
        elif self.qmmm_param['hl_program'] == "OpenMM":
            self.hl_param = self.load_param(self.openmm_paramfile)
            self.hl_wrapper = OpenMMWrapper
        else:
            print("Only Psi4 and OpenMM currently available to be used in high level computations")

        if self.qmmm_param['ll_program'] == "OpenMM":
            self.ll_param = self.load_param(self.openmm_paramfile)
            self.ll_wrapper = OpenMMWrapper
        else:
            print("Only OpenMM currently available to be used in low level computations")

        if 'hl' in self.param:
            self.hl_param.update(self.param['hl'])

        if 'll' in self.param:
            self.ll_param.update(self.param['ll'])

        if 'md' in self.param:
            self.md_sim_param.update(self.param['md'])

        self.run_md = self.md_sim_param['run_md']
        self.md_restart = self.md_sim_param['restart']

        if 'return_forces_interval' in self.md_sim_param:
            self.return_forces_interval = self.md_sim_param['return_forces_interval']
        else:
            self.return_forces_interval = self.md_sim_param['return_checkpoint_interval']
    
        self.return_forces_filename = self.md_sim_param['return_forces_filename']
            
        self.md_sim_wrapper = None

        if self.run_md is True:
            self.md_sim_prog = self.md_sim_param['md_simulation_program']
            self.start_qmmm =  self.md_sim_param['start_qmmm']
            self.end_qmmm =    self.md_sim_param['end_qmmm']
            self.qmmm_steps = self.end_qmmm - self.start_qmmm

            if type(self.md_sim_param['md_steps']) is int:
                self.end_steps = self.md_sim_param['md_steps'] - self.end_qmmm
            elif type(self.md_sim_param['md_steps']) is list:
                self.end_steps = self.md_sim_param['md_steps'][-1] - self.end_qmmm

            if self.md_sim_prog == "OpenMM":
                print('reading openmm as md wrapper')
                self.md_sim_wrapper = OpenMMWrapper
            else:
                print("Only OpenMM currently available")

        self.qmmm_param.update(self.system_param)
        self.aqmmm_param.update(self.system_param)
        self.md_sim_param.update(self.system_param)
        self.md_sim_param.update(self.ll_param)
        self.ll_param.update(self.md_sim_param)

        #print('System information read from {}'.format(self.param['system']['system_info']))

    def initialize_wrappers(self, simulation=False, restart=False):
        """
        Instantiates qm, mm, qmmm, and/or aqmmm wrapper objects 
        used for computation based on input parameters


        Returns
        -------
        MMWrapper object
        QMMM wrapper object

        """

        # create hl_wrapper object
        hl_wrapper = self.hl_wrapper(self.hl_param)
        # create ll_wrapper object
        ll_wrapper = self.ll_wrapper(self.ll_param)
        
        if self.run_md is True:
            # create md wrapper
            md_sim_wrapper = self.md_sim_wrapper(self.md_sim_param)

        if self.qmmm_param['run_aqmmm'] is False:
            qmmm = QMMM(self.qmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        elif self.aqmmm_param['aqmmm_scheme'] == 'ONIOM-XS':
            qmmm = OniomXS(self.aqmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        elif self.aqmmm_param['aqmmm_scheme'] == 'Hot-Spot':
            qmmm = HotSpot(self.aqmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        elif self.aqmmm_param['aqmmm_scheme'] == 'PAP':
            qmmm = PAP(self.aqmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        elif self.aqmmm_param['aqmmm_scheme'] == 'SAP':
            qmmm = SAP(self.aqmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        elif self.aqmmm_param['aqmmm_scheme'] == 'DAS':
            qmmm = SAP(self.aqmmm_param, hl_wrapper, ll_wrapper, self.md_sim_prog)
        else:
            print("Only ONIOM-XS, Hot Spot, PAP, SAP, and DAS currently implemented")

        if (simulation is True and restart is False):
            # initialize mm_wrapper with information about initial system
            md_sim_wrapper.initialize(self.qmmm_param['embedding_method'])
            return md_sim_wrapper, qmmm
 
        elif (simulation is True and restart is True):
            # initialize mm_wrapper with information about initial system
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


