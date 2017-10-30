"""
Testing for the parser.py module
"""
import pytest
import janus

@pytest.mark.datafiles('tests/examples/test_parser/input.dat')
def test_input_parser(datafiles):
    #assert (datafiles / 'input.dat').check(file=1)
    system = janus.parser.parse_input(datafiles / 'input.dat')
    assert system.qm_param['basis'] == 'STO-3G' 
    assert system.qm_param['reference'] == 'rhf'     
    assert system.qm_param['scf_type'] == 'df'
    assert system.qm_param['guess'] == 'sad'
    assert system.qm_param['e_convergence'] == '1e-8'
    assert system.qm_param['d_convergence'] == '1e-8'
    assert system.qm_method == 'scf'
    assert system.qm_molecule == '''
0 1
O
H 1 R
H 1 R 2 A
R = 1.0
A = 104.5
symmetry c1'''
