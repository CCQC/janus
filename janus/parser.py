import configparser 

def parse_input(filename, system):

    config = configparser.ConfigParser()
    config.read(filename)

    system.qm_param = {}
    system.qm_param['basis'] = config['QM_PARAM']['basis']
    system.qm_param['reference'] = config['QM_PARAM']['reference']
#    system.qm_param['max_iter'] = config['QM_PARAM']['max_iter']
    system.qm_param['scf_type'] = config['QM_PARAM']['scf_type']
    system.qm_param['e_convergence'] = config['QM_PARAM']['e_convergence']
    system.qm_param['d_convergence'] =  config['QM_PARAM']['d_convergence']
    system.qm_method = config['QM_PARAM']['method']
    system.qm_molecule = config['MOLECULE']['molecule']
