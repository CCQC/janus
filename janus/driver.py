"""
This is the qmmm driver module
"""
import json
from .initializer import Initializer

def run_simulation(filename='input.json', equilibrate=5000):
    """
    Drives QM/MM with MD time step integration
    
    Parameters
    ----------
    filename : str 
        contains the filename of the input file, default is 'input.json'

    equilibrate : int
        number of steps to run MD steps before QM/MM, default is 5000

    """

    initializer = Initializer(filename)

    # initialize wrappers
    md_simulation_wrapper, qmmm = initializer.initialize_wrappers(simulation=True)

    md_simulation_wrapper.equilibrate(equilibrate)

    for step in range(initializer.steps):

        #get MM information for entire system
        main_info = md_simulation_wrapper.get_main_info()

        qmmm.run_qmmm(main_info)
        
        # get aqmmm forces 
        forces = qmmm.get_forces()
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        md_simulation_wrapper.take_step(force=forces)

def run_single_point(filename='input.json'):
    """
    Drives single QM/MM computation

    Parameters
    ----------
    filename : str 
        contains the filename of the input file, default is 'input.json'

    """
    initializer = Initializer(filename)
    ll_wrapper, qmmm = initializer.initialize_wrappers()

    #get MM information for entire system
    main_info = ll_wrapper.get_main_info()

    qmmm.run_qmmm(main_info)

