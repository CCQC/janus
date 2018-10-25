"""
This is the qmmm driver module
"""
import json
from .initializer import Initializer
    
def run_janus(filename='input.json'):

    initializer = Initializer(filename)

    if initialize.run_md is True:
        run_simulation(initializer)

    else:
        run_single_point(initializer)

def run_simulation(initializer):
    """
    Drives QM/MM with MD time step integration
    
    Parameters
    ----------

    equilibrate : int
        number of steps to run MD steps before QM/MM, default is 5000

    """

    # initialize wrappers
    md_simulation_wrapper, qmmm = initializer.initialize_wrappers(simulation=True)

    md_simulation_wrapper.take_step(initializer.start_qmmm)

    for step in range(initializer.qmmm_steps):

        #get MM information for entire system
        main_info = md_simulation_wrapper.get_main_info()

        qmmm.run_qmmm(main_info)
        
        # get aqmmm forces 
        forces = qmmm.get_forces()
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        md_simulation_wrapper.take_updated_step(force=forces)

    md_simulation_wrapper.take_step(initializer.end_steps)

    # here put return pdb??

def run_single_point(initializer):
    """
    Drives single QM/MM computation

    Parameters
    ----------
    filename : str 
        contains the filename of the input file, default is 'input.json'

    """
    ll_wrapper, qmmm = initializer.initialize_wrappers()

    #get MM information for entire system
    main_info = ll_wrapper.get_main_info()

    qmmm.run_qmmm(main_info)


