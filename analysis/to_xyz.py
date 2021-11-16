import sys, os
import json
import numpy as np

ATOMICLABELS = ['H',  'He',  'Li',  'Be',  'B',  'C',  'N',  'O',  'F',  'Ne',  'Na',  'Mg',  'Al',  'Si',  'P',  'S', 
    'Cl',  'Ar',  'K',  'Ca',  'Sc',  'Ti',  'V',  'Cr',  'Mn',  'Fe',  'Co',  'Ni',  'Cu',  'Zn',  'Ga',  'Ge',  'As',  'Se',  
    'Br',  'Kr',  'Rb',  'Sr',  'Y',  'Zr',  'Nb',  'Mo',  'Tc',  'Ru',  'Rh',  'Pd',  'Ag',  'Cd',  'In',  'Sn',  'Sb',  'Te', 
    'I',  'Xe',  'Cs',  'Ba',  'La',  'Ce',  'Pr',  'Nd',  'Pm',  'Sm',  'Eu',  'Gd',  'Tb',  'Dy',  'Ho',  'Er',  'Tm',  'Yb', 
    'Lu',  'Hf',  'Ta',  'W',  'Re',  'Os',  'Ir',  'Pt',  'Au',  'Hg',  'Tl',  'Pb',  'Bi',  'Po',  'At',  'Rn',  'Fr',  'Ra', 
    'Ac',  'Th',  'Pa',  'U',  'Np',  'Pu',  'Am',  'Cm',  'Bk',  'Cf',  'Es',  'Fm',  'Md',  'No',  'Lr',  'Rf',  'Db',  'Sg', 
    'Bh',  'Hs',  'Mt',  'Ds',  'Rg',  'Cn',  'Nh',  'Fl',  'Mc',  'Lv',  'Ts',  'Og'] # Yes *all* of the elements

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