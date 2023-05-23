import sys, os
import json
import numpy as np
from defs import ATOMICLABELS

if len(sys.argv) < 3 :
    print(f'Use: {sys.argv[0]} /path/to/manifest /path/to/output_traj')
    sys.exit()
_, path_to_manifest, output = sys.argv

assert os.path.exists(path_to_manifest)
with open(path_to_manifest, 'r') as f:
    manifest = json.load(f)

assert('xyz' in manifest['quantities'])
assert('t' in manifest['quantities'])

basepath = os.path.dirname(path_to_manifest)

xyz = np.load(os.path.join(basepath, 'xyz_ave'))
times = np.load(os.path.join(basepath, 'times'))
natoms = len(manifest['atomnos'])

with open(output, 'w') as f:
    for i, t in enumerate(times):
        f.write(f'{natoms}\nt=\t{t}\n')
        for na, a in enumerate(xyz[i]):
            f.write('{}\t\t{:>6.5f}\t\t{:>6.5f}\t\t{:>6.5f}\n'.format(ATOMICLABELS[manifest['atomnos'][str(na+1)]-1], a[0], a[1], a[2]))

print("Write xyz OK")