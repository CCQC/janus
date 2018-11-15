"""
This is the qmmm driver module
"""
import pickle
from janus import Initializer
    
def run_janus(filename='input.json'):

    initializer = Initializer(filename)

    if initializer.run_md is True:
        run_simulation(initializer, restart=initializer.md_restart)

    else:
        run_single_point(initializer)

def run_simulation(initializer, restart):
    """
    Drives QM/MM with MD time step integration
    
    Parameters
    ----------

    equilibrate : int
        number of steps to run MD steps before QM/MM, default is 5000

    """

    print('Initializing')
    # initialize wrappers
    md_simulation_wrapper, qmmm = initializer.initialize_wrappers(simulation=True, restart=restart)

    print('Equilibrating with {} steps'.format(initializer.start_qmmm))
    md_simulation_wrapper.take_step(initializer.start_qmmm)

    for step in range(initializer.qmmm_steps):

        print('Taking step {}'.format(step + 1))
        #get MM information for entire system
        main_info = md_simulation_wrapper.get_main_info()

        print('Running QMMM for step {}'.format(step + 1))
        qmmm.run_qmmm(main_info)
        
        # get aqmmm forces 
        forces = qmmm.get_forces()

        if (step + 1) % initializer.return_forces_interval == 0:
            with open(initializer.return_forces_filename, 'wb') as f:
                pickle.dump(forces, f)

    
        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        md_simulation_wrapper.take_updated_step(force=forces)

    print('QMMM finished')

    md_simulation_wrapper.take_step(initializer.end_steps)

    main_info = md_simulation_wrapper.get_main_info()
    md_simulation_wrapper.write_pdb(main_info)


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


