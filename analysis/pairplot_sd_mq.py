import os
import numpy as np
import sys
from itertools import combinations
import json
from defs import ATOMICLABELS


try: _, mode, json_path, output = sys.argv
except:
    print(f'Use : python {sys.argv[0]} mq/sd path/to/manifest.json output')
    sys.exit(-1)

assert(os.path.exists(json_path))
with open(json_path, 'r') as f:
    manifest = json.load(f)

if mode=='sd':
    assert('sd' in manifest['quantities'])
elif mode=='mq':
    assert('mq' in manifest['quantities'])
else:
    raise Exception('Illegal mode [must be one of either nm/mq]')

basepath = os.path.dirname(json_path)
datafile = 'sd_ave'if mode=='sd' else 'mq_ave'
print(datafile)
data = np.load(os.path.join(basepath, datafile))

times = np.load(os.path.join(basepath, 'times'))
print(data.shape)

pairs = list(combinations(range(data.shape[0]), 2))

dim = 0
while len(pairs) > (2+dim)*(3+dim):
    dim += 1


import matplotlib.pyplot as plt
fig_idx = 0
fig, axes = plt.subplots(2+dim, 3+dim , num=f'{json_path} {mode}')

for p in pairs:
    # Figure out which plot to use
    fig_x = np.floor(fig_idx / (3+dim)).astype(int)
    fig_y = fig_idx % (3+dim)
    print(f'[{fig_y}, {fig_x}]')

    if mode=='sd':
        axes[fig_x, fig_y].set_ylabel('Spin density')
    elif mode=='mq':
        axes[fig_x, fig_y].set_ylabel('Mulliken charge')
    axes[fig_x, fig_y].set_xlabel('Time (fs)')

    # Work out the atomic symbols
    an1 = manifest['spindenmap'][p[0]]
    try: sym1 = ATOMICLABELS[manifest['atomnos'][str(an1)]-1]
    except: sym1 = '?'

    an2 = manifest['spindenmap'][p[1]]
    try: sym2 = ATOMICLABELS[manifest['atomnos'][str(an2)]-1]
    except: sym2 = '?'

    axes[fig_x, fig_y].plot(times, data[p[0]], label='{} [{}]'.format(sym1, an1))
    axes[fig_x, fig_y].plot(times, data[p[1]], label='{} [{}]'.format(sym2, an2))

    axes[fig_x, fig_y].legend(loc='upper right');
    fig_idx += 1
    

if output=='x11' : plt.show()
else: plt.savefig(output)
