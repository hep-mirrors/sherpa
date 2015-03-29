from os import path
from test_model import test_model

sm   = test_model('sm',   path.abspath('./'),      'param_card_sm.dat')
heft = test_model('heft', path.abspath('./'),    'param_card_heft.dat')
mued = test_model('mued', path.abspath('./'),    'param_card_mued.dat')
mssm = test_model('mssm', path.abspath('./'), 'param_card_mssm.dat')
qagc = test_model('qagc', path.abspath('./'),    'param_card_qagc.dat')

