import configparser 

def parse_input(filename, sys):
    """
    This function is a parser that reads an input and
    stores the information into a system object
    
    Parameters
    ----------
    filename: an input file 
    sys: System class object

    Returns
    -------
    None

    Examples
    --------
    parse_input('input.dat', system)
    """
    config = configparser.ConfigParser()
    config.read(filename)

    sys.qm_param = {}
    sys.qm_param['basis'] = config['QM_PARAM']['basis']
    sys.qm_param['reference'] = config['QM_PARAM']['reference']
#   sys.qm_param['max_iter'] = config['QM_PARAM']['max_iter']
    sys.qm_param['scf_type'] = config['QM_PARAM']['scf_type']
    sys.qm_param['e_convergence'] = config['QM_PARAM']['e_convergence']
    sys.qm_param['d_convergence'] =  config['QM_PARAM']['d_convergence']
    sys.qm_method = config['QM_PARAM']['method']
    sys.qm_molecule = config['MOLECULE']['molecule']
