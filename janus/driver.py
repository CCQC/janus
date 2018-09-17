"""
This is the qmmm driver module
"""
import json
from .initializer import Initializer

def run_janus(filename):

    initializer = Initializer(filename)

    # initialize wrappers
    mm_wrapper, qmmm = initializer.initialize_wrappers()


    for step in range(config['steps']):

        #get MM information for entire system
        main_info = mm_wrapper.get_main_info()

        qmmm.run_qmmm(main_info)
        
        # get aqmmm forces 
        forces = qmmm.get_forces()
        print('forces', forces)
    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        mm_wrapper.take_step(force=forces)

