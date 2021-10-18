import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
import mathutils

# Welcome to jank city - things may or may not work
# In code gwps are 0-indexed - like numpy arrays - but users see 1-indexed notation
ATOMICLABELS = ['H',  'He',  'Li',  'Be',  'B',  'C',  'N',  'O',  'F',  'Ne',  'Na',  'Mg',  'Al',  'Si',  'P',  'S', 
    'Cl',  'Ar',  'K',  'Ca',  'Sc',  'Ti',  'V',  'Cr',  'Mn',  'Fe',  'Co',  'Ni',  'Cu',  'Zn',  'Ga',  'Ge',  'As',  'Se',  
    'Br',  'Kr',  'Rb',  'Sr',  'Y',  'Zr',  'Nb',  'Mo',  'Tc',  'Ru',  'Rh',  'Pd',  'Ag',  'Cd',  'In',  'Sn',  'Sb',  'Te', 
    'I',  'Xe',  'Cs',  'Ba',  'La',  'Ce',  'Pr',  'Nd',  'Pm',  'Sm',  'Eu',  'Gd',  'Tb',  'Dy',  'Ho',  'Er',  'Tm',  'Yb', 
    'Lu',  'Hf',  'Ta',  'W',  'Re',  'Os',  'Ir',  'Pt',  'Au',  'Hg',  'Tl',  'Pb',  'Bi',  'Po',  'At',  'Rn',  'Fr',  'Ra', 
    'Ac',  'Th',  'Pa',  'U',  'Np',  'Pu',  'Am',  'Cm',  'Bk',  'Cf',  'Es',  'Fm',  'Md',  'No',  'Lr',  'Rf',  'Db',  'Sg', 
    'Bh',  'Hs',  'Mt',  'Ds',  'Rg',  'Cn',  'Nh',  'Fl',  'Mc',  'Lv',  'Ts',  'Og'] # Yes *all* of the elements

def plotnms(basepath, manifest, nms):
    assert('nm' in manifest['quantities'])
    times = np.load(os.path.join(basepath, 'times'))
    data = np.load(os.path.join(basepath, 'nm'))
    fig, ax = plt.subplots()
    for gwp in range(data.shape[0]):
        raw_data = data[gwp].T
        
        for i in nms:
            ax.plot(times, raw_data[i], label=f'NM{i+1}')
        ax.set_title(f'Normal mode evolution (for GWP{gwp+1})')
        ax.set_ylabel('Normal mode evolution')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join(basepath, f'dbg_nms_gwp{gwp+1}.png'))
        plt.cla()
    plt.close(fig)

def plotbl(basepath, manifest, BPS):
    assert('xyz' in manifest['quantities'])
    data = np.load(os.path.join(basepath, 'xyz'))
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for gwp in range(data.shape[0]):
        raw_data = data[gwp]
        for a in BPS:
            dp = []
            for x in range(manifest['steps']):
                dp.append(mathutils.MathUtils.bond_length(raw_data[x, a[0]-1],raw_data[x, a[1]-1] ))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'

            ax.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
        ax.set_title(f'BL evolution (for GWP{gwp+1})')
        ax.set_ylabel('Bond length (Å)')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join(basepath, f'dbg_bl_gwp{gwp+1}.png'))
        plt.cla()
    plt.close(fig)

if len(sys.argv) < 3:
    print(f'Use {sys.argv[0]} path/to/manifest.json nms bls')

_, manifestpath, nms, bls = sys.argv
basepath = os.path.dirname(manifestpath)
try:
    with open(manifestpath, 'r') as f:
        manifest = json.load(f)
except:
    print('Unable to open manifest!')
    sys.exit(-1)

if nms == 'None': nms = None
else: nms = [int(i)-1 for i in nms.split(',')]

if bls == 'None' : BPS = None
else:
    BPS = []
    for x in bls.split(','):
        a, b = [int(z) for z in x.split('-')]
        BPS.append([a,b])

if nms is not None:
    print(f'Plotting up normal modes {nms}')
    plotnms(basepath, manifest, nms)
if BPS is not None:
    print('Plotting up bond lengths {}'.format([f'{i[0]} - {i[1]}' for i in BPS]))
    plotbl(basepath, manifest,BPS)
