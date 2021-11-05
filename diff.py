import sys
import os
from quatics_lexers import QuanticsParsers
from logall import ParseLogAll
import numpy as np
import copy

def compare(n1, n2, err=1E-9):
    status=0
    assert (n1.shape==n2.shape)
    for i in range(n1.size):
        if np.abs(n1.flatten()[i] - n2.flatten()[i]) > err:
            print(f'A large error occured {n1.flatten()[i]} {n2.flatten()[i]} | err>{err}')
            status = 1
    return status

if len(sys.argv) < 3 :
    print(f'Use {sys.argv[0]} dir1/path/to/inp dir2/path/to/inp')
    sys.exit(-1)

_, dir1, dir2 = sys.argv

assert(os.path.exists(dir1))
assert(os.path.exists(dir2))

# Parse the quantics input
with open(dir1, 'r') as f:
    data1 = f.read()
q_inp_data_1=QuanticsParsers.parse_input(data1)
print(f'QINP file 1 parsed OK')
with open(dir2, 'r') as f:
    data2 = f.read()
q_inp_data_2=QuanticsParsers.parse_input(data2)
print(f'QINP file 2 parsed OK')

# Quick check that the output path exists
datadir1 = os.path.join(os.path.dirname(dir1), q_inp_data_1['data'])
print(datadir1)
assert(os.path.exists(datadir1))

datadir2 = os.path.join(os.path.dirname(dir2), q_inp_data_2['data'])
print(datadir2)

# Load in all GWP logalls
print('Parsing data 1')
data_gwpx1 = ParseLogAll.I_ImportLogalls(datadir1, q_inp_data_1['ngwp'], quantities=['case', 'ci', 'csf','fo'], fname='gwp{}_V1_gaussian_data.logall')
print('Parsing data 2')
data_gwpx2 = ParseLogAll.I_ImportLogalls(datadir2, q_inp_data_2['ngwp'], quantities=['case', 'ci', 'csf','fo'], fname='gwp{}_V1_gaussian_data.logall')

# Do the comparisons
state = 0
print('Compare CAS Energy')
state += compare(data_gwpx1['case'], data_gwpx2['case'])
print('Compare adibatic populations')
state += compare(data_gwpx1['adiabats'], data_gwpx2['adiabats'])
print('Compare diabatic populations')
state += compare(data_gwpx1['diabats'], data_gwpx2['diabats'])
print('Compare gradients')
state += compare(data_gwpx1['forces'], data_gwpx2['forces'])
sys.exit(state)