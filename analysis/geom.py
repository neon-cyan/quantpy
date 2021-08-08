import mathutils
import numpy as np
import sys
import os
import json

# This is the user-run script
#  Currently it plots all the properties available

# TODO : argparse - Seperate plots
# TODO : dataloader to avoid re-reading logall files if possible
# TODO : bond angles / dihedral angles

if len(sys.argv) < 5:
    print(f"""Not enough arguemnts!\n Use : {sys.argv[0]}
    [/path/manifest.json] 
    [mode = (P)BL | BA | DA]
    [what_to_plot 1-2 | 1-2-3 | 1-2-3-4]
    [Output dims w,h]
    [Output x11 | filename.eps]
    """)
    sys.exit(-1)

ATOMICLABELS = ['H',  'He',  'Li',  'Be',  'B',  'C',  'N',  'O',  'F',  'Ne',  'Na',  'Mg',  'Al',  'Si',  'P',  'S', 
    'Cl',  'Ar',  'K',  'Ca',  'Sc',  'Ti',  'V',  'Cr',  'Mn',  'Fe',  'Co',  'Ni',  'Cu',  'Zn',  'Ga',  'Ge',  'As',  'Se',  
    'Br',  'Kr',  'Rb',  'Sr',  'Y',  'Zr',  'Nb',  'Mo',  'Tc',  'Ru',  'Rh',  'Pd',  'Ag',  'Cd',  'In',  'Sn',  'Sb',  'Te', 
    'I',  'Xe',  'Cs',  'Ba',  'La',  'Ce',  'Pr',  'Nd',  'Pm',  'Sm',  'Eu',  'Gd',  'Tb',  'Dy',  'Ho',  'Er',  'Tm',  'Yb', 
    'Lu',  'Hf',  'Ta',  'W',  'Re',  'Os',  'Ir',  'Pt',  'Au',  'Hg',  'Tl',  'Pb',  'Bi',  'Po',  'At',  'Rn',  'Fr',  'Ra', 
    'Ac',  'Th',  'Pa',  'U',  'Np',  'Pu',  'Am',  'Cm',  'Bk',  'Cf',  'Es',  'Fm',  'Md',  'No',  'Lr',  'Rf',  'Db',  'Sg', 
    'Bh',  'Hs',  'Mt',  'Ds',  'Rg',  'Cn',  'Nh',  'Fl',  'Mc',  'Lv',  'Ts',  'Og'] # Yes *all* of the elements

manifest_path = sys.argv[1]
mode = sys.argv[2]

BOND_PAIRS = []
ness_len = 2 # Assume (percentage) bond length
if mode == 'BA' : ness_len = 3
elif mode == 'DA' : ness_len = 4

for x in sys.argv[3].split(','):
    a = [int(z) for z in x.split('-')]
    assert(len(a) == ness_len)
    BOND_PAIRS.append(a)

OUTDIMS = [float(i) for i in sys.argv[4].split(',')]
OUTPUT = sys.argv[5] # Either plot to a file or x11 window

assert(os.path.exists(manifest_path))
with open(manifest_path, 'r') as f:
    manifest = json.load(f)
basepath = os.path.dirname(manifest_path)

times = np.load(os.path.join(basepath, 'times'))
avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
nsteps = manifest['steps']

####################### PLOTTING HAPPENS HERE

import matplotlib.pyplot as plt

fig =  plt.figure(num=manifest_path, figsize=OUTDIMS)

if mode == 'BL':
    for a in BOND_PAIRS:
        dp = []
        for x in range(nsteps):
            dp.append(mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] ))

        try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
        except: alab1 = '?'
        try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
        except: alab2 = '?'

        plt.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
    plt.ylabel('Bond length (Å)')


elif mode == 'PBL':
    for a in BOND_PAIRS:
        dp = []
        init_bl = mathutils.MathUtils.bond_length(avegeom[0, a[0]-1],avegeom[0, a[1]-1] )
        for x in range(nsteps):
            bl = mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] )
            dp.append((bl - init_bl) / init_bl)

        try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
        except: alab1 = '?'
        try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
        except: alab2 = '?'

        plt.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
    plt.ylabel('Fractional change in bond length')

elif mode == 'BA':
    for a in BOND_PAIRS:
        dp = []
        for x in range(nsteps):
            dp.append(mathutils.MathUtils.bond_angle(avegeom[x, a[0]-1],avegeom[x, a[1]-1], avegeom[x, a[2]-1] ))

        try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
        except: alab1 = '?'
        try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
        except: alab2 = '?'
        try: alab3 = ATOMICLABELS[manifest['atomnos'][str(a[2])]-1]
        except: alab3 = '?'

        plt.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}]')
    plt.ylabel('Bond angle (rad)')

elif mode == 'DA':
    for a in BOND_PAIRS:
        dp = []
        for x in range(nsteps):
            dha = mathutils.MathUtils.dihedral([avegeom[x, a[0]-1],avegeom[x, a[1]-1], avegeom[x, a[2]-1], avegeom[x, a[3]-1] ])
            dp.append(dha)

        try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
        except: alab1 = '?'
        try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
        except: alab2 = '?'
        try: alab3 = ATOMICLABELS[manifest['atomnos'][str(a[2])]-1]
        except: alab3 = '?'
        try: alab4 = ATOMICLABELS[manifest['atomnos'][str(a[3])]-1]
        except: alab4 = '?'

        plt.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}] - {alab4}[{a[3]}]')
    plt.ylabel('Dihedral angle (rad)')

plt.title(f'{mode}')
plt.xlabel('Time (fs)')
plt.legend(loc='lower right');

if OUTPUT=='x11' : plt.show()
else: plt.savefig(OUTPUT, dpi=600)
