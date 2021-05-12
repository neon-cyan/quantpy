import mathutils
from glogpy.freqency_job import frequency_job
from quatics_lexers import QuanticsParsers
import logall
import numpy as np
import sys
import os

# This is the user-run script
#  Currently it plots all the properties available

# TODO : argparse - Seperate plots
# TODO : dataloader to avoid re-reading logall files if possible
# TODO : bond angles / dihedral angles


QINP = sys.argv[1]
STEPLIMS = int(sys.argv[2]) if sys.argv[2] != 'A' else None
OUTPUT = 'x11' # Either plot to a file or x11 window
ANNA_DIR = 'analysis'
NM_PLOTS = [int(x) for x in sys.argv[3].split(',')]
BOND_PAIRS = []
for x in sys.argv[4].split(','):
    a, b = [int(z) for z in x.split('-')]
    BOND_PAIRS.append([a,b])

CI_STATES = None
if sys.argv[5] != 'A': CI_STATES = [int(x) for x in sys.argv[5].split(',')]

CSF_STATES = None
if sys.argv[6] != 'A': CSF_STATES = [int(x) for x in sys.argv[5].split(',')]

assert(os.path.exists(QINP))
with open(QINP, 'r') as f:
    data = f.read()
q_inp_data=QuanticsParsers.parse_input(data)
print('Read Quantics input OK')

datadir = os.path.join(os.path.dirname(QINP), q_inp_data['data'])
assert(os.path.exists(datadir))

# Load up quantics output file
qoutf = os.path.join(os.path.dirname(QINP), q_inp_data['name'], 'output')
assert(os.path.exists(qoutf))
with open(qoutf, 'r') as f:
    data = f.read()
q_out_data=QuanticsParsers.parse_output(data)
if STEPLIMS!=None: q_out_data = q_out_data[:STEPLIMS]
print('Read Quantics output OK')

#  Load up frwquency file
freqf = os.path.join(datadir, q_inp_data['freqf'])
assert(os.path.exists(freqf))
with open(freqf, 'r') as f:
    data = f.read()
freq = frequency_job(data)
freq_data = freq.parse()
print('Read Freq file OK')


# Load in all GWP logalls
data_gwpx = logall.I_ImportLogalls(datadir, q_inp_data['ngwp'], step_lim=STEPLIMS)
print('Read GWP files OK')
nsteps = data_gwpx['steps']
assert(nsteps == len(q_out_data))

# print(freq_data['geom'])
geom_init = mathutils.MathUtils.dict_to_list(freq_data['geom'])
geom_init = np.array([x[1] for x in geom_init])

# Compute nm2xtyz and xyz2nm matrices
nm2xyz, xyz2nm = mathutils.NormModeUtils.nm_matrix(data_gwpx['atommasses'], 
freq_data['vibfreqs'], freq_data['vibdisps'])

# TODO Add option to save the matrices to avoid re-computing
nmdata = mathutils.NormModeUtils.xyz_to_nm(xyz2nm,geom_init, data_gwpx['geomx'])

# print(nmdata)
# print(nmdata.shape)
nnmode = nmdata.shape[2]
nmdata  = nmdata.transpose(1,0,2)

# print('REEE')
times = np.array([x['time'] for x in q_out_data])
gwp_sf = np.array([x['GGP'] for x in q_out_data])
gwp_sf /= 10 # Sum up to 10 - not 1

# I HAVE NO IDEA HOW TO USE NP.DOT // PLEASE NO BULLY
res = np.zeros((nsteps, freq_data['vibdisps'].shape[0])) # NM x S matrix
for i in range(res.shape[0]):
    res[i] = gwp_sf[i].dot(nmdata[i])

res = res.T
# print(res)
print(f'res = {res.shape}')
print(f'nm2xyz = {nm2xyz.shape}')

avegeom = np.zeros((nsteps, geom_init.shape[0], 3))
for i in range(nsteps):
    avegeom[i] = np.copy(geom_init)
    for j in range(nnmode):
        avegeom[i] += res[j, i] * nm2xyz.T[j].reshape(geom_init.shape[0], 3)

# print(avegeom)
####################### PLOTTING HAPPENS HERE

import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, num=QINP)
#   Looks like
#   FIG1 FIG2 FIG3
#   FIG4 FIG5 FIG6

for x in NM_PLOTS:
    axes[0,0].plot(times, res[x-1], label=f'NM{x}')
axes[0,0].set_title('Normal modes')
axes[0,0].set_ylabel('Normal mode excitation')
axes[0,0].set_xlabel('Time (fs)')
axes[0,0].legend(loc='lower left');

for a in BOND_PAIRS:
    dp = []
    for x in range(nsteps):
        dp.append(mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] ))
    axes[0,1].plot(times, dp, label=f'BD {a[0]} - {a[1]}')
axes[0,1].set_title('Bond lengths')
axes[0,1].set_ylabel('Bond length (Å)')
axes[0,1].set_xlabel('Time (fs)')
axes[0,1].legend(loc='lower left');

# CI State evolution
adiabats = data_gwpx['adiabats'].transpose(2,1,0)
print(f'ADB = {adiabats.shape}')
print(f'GWPS = {gwp_sf.shape}')
axes[1,0].set_title('Adiabatic [CI] state evolution')
for i in range(adiabats.shape[0]):
    if CI_STATES == None: pass
    else:
        if i+1 in CI_STATES: continue

    dp = np.zeros(nsteps)
    for x in range(nsteps):
        dp[x] = np.square(abs(adiabats[i, x]).dot(gwp_sf[x]))
    axes[1,0].plot(times, dp, label=f'CI {i+1}')
axes[1,0].set_ylabel('State population')
axes[1,0].set_xlabel('Time (fs)')
axes[1,0].legend(loc='lower right');

# CSF State evoluition
diabats = data_gwpx['diabats'].transpose(2,1,0)
axes[1,1].set_title('Diabatic [CSF] state evolution')
for i in range(diabats.shape[0]):
    if CSF_STATES == None: pass
    else:
        if i+1 in CSF_STATES: continue

    dp = np.zeros(nsteps)
    for x in range(nsteps):
        dp[x] = np.square(abs(diabats[i, x]).dot(gwp_sf[x]))
    axes[1,1].plot(times, dp, label=f'CSF {i+1}')
axes[1,1].set_ylabel('State population')
axes[1,1].set_xlabel('Time (fs)')
axes[1,1].legend(loc='lower right');

# MULLIKEN & SPINDENS
# Mulliken summed
msum = data_gwpx['mullikensum'].transpose(2,1,0)
axes[0,2].set_title('Mulliken charge evolution (H Summed)')
for i in range(msum.shape[0]):
    dp = np.zeros(nsteps)
    for x in range(nsteps):
        dp[x] = abs(msum[i, x]).dot(gwp_sf[x])
    axes[0,2].plot(times, dp, label='Atom {}'.format(data_gwpx['mullikenmap'][i]))
axes[0,2].set_ylabel('Mulliken charge')
axes[0,2].set_xlabel('Time (fs)')
axes[0,2].legend(loc='lower right');

# Spin densities summed
sdsum = data_gwpx['spindensum'].transpose(2,1,0)
axes[1,2].set_title('Spin density evolution (H Summed)')
for i in range(sdsum.shape[0]):
    dp = np.zeros(nsteps)
    for x in range(nsteps):
        dp[x] = abs(sdsum[i, x]).dot(gwp_sf[x])
    axes[1,2].plot(times, dp, label='Atom {}'.format(data_gwpx['spindenmap'][i]))
axes[1,2].set_ylabel('Spin density')
axes[1,2].set_xlabel('Time (fs)')
axes[1,2].legend(loc='lower right');

if OUTPUT=='x11' : plt.show()
else: plt.savefig(OUTPUT)
