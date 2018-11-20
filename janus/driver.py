"""
This is the qmmm driver module
"""
import pickle
from janus import Initializer
    
def run_janus(filename='input.json'):

    initializer = Initializer(filename)

    print('Initializing')
    # initialize wrappers
    ll_wrapper, qmmm_wrapper = initializer.initialize_wrappers()

    if initializer.run_md is True:
        run_simulation(ll_wrapper, qmmm_wrapper)
    else:
        run_single_point(ll_wrapper, qmmm_wrapper)

def run_simulation(md_sim_wrapper, qmmm_wrapper):
    """
    Drives QM/MM with MD time step integration
    
    Parameters
    ----------

    equilibrate : int
        number of steps to run MD steps before QM/MM, default is 5000

    """

    print('Equilibrating with {} steps'.format(md_sim_wrapper.start_qmmm))
    md_sim_wrapper.take_step(mm_sim_wrapper.start_qmmm)

    for step in range(md_sim_wrapper.qmmm_steps):

        print('Taking step {}'.format(step + 1))
        run_single_point(md_sim_wrapper, qmmm_wrapper)
        
        # get aqmmm forces 
        forces = qmmm_wrapper.get_forces()

        if (step + 1) % md_sim_wrapper.return_forces_interval == 0:
            with open(md_sim_wrapper.return_forces_filename, 'wb') as f:
                pickle.dump(forces, f)

        # feed forces into md simulation and take a step
        # make sure positions are updated so that when i get information on entire system 
        # getting it on the correct one
        md_sim_wrapper.take_updated_step(force=forces)

    print('QMMM finished')

    md_sim_wrapper.take_step(md_sim_wrapper.end_steps)

    main_info = md_sim_wrapper.get_main_info()
    md_sim_wrapper.write_pdb(main_info)


def run_single_point(ll_wrapper, qmmm_wrapper):
    """
    Drives single QM/MM computation

    Parameters
    ----------
    filename : str 
        contains the filename of the input file, default is 'input.json'

    """
    #get MM information for entire system
    main_info = ll_wrapper.get_main_info()

    qmmm_wrapper.run_qmmm(main_info, ll_wrapper.class_type)


