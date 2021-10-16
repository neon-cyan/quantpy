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


def find_probelm_gwps(basepath, manifest):
    assert('casde' in manifest['quantities'])
    assert('nucde' in manifest['quantities'])
    assert('maxf' in manifest['quantities'])
    timestep=manifest['tout']

    DEFAULT_CAST = 0.01
    DEFAULT_NUCT = 10.0 
    DEFAULT_FORCET = 10.0

    print(f"Please input CASSCF Delta-E thereshold (default = {DEFAULT_CAST})")
    casde_thresh=input()
    casde_thresh = DEFAULT_CAST if casde_thresh=='' else float(casde_thresh)

    print(f"Please input QuEh gradient thereshold (default = {DEFAULT_FORCET})")
    maxf_thresh=input()
    maxf_thresh = DEFAULT_FORCET if maxf_thresh=='' else float(maxf_thresh)

    print(f"Please input vMCG Delta-E thereshold (default = {DEFAULT_NUCT})")
    nucde_thresh=input()
    nucde_thresh = DEFAULT_NUCT if nucde_thresh=='' else float(nucde_thresh)

    prob_gwps=[]
    with open(os.path.join(basepath, 'casde') , 'rb') as f:
        data = np.load(f)
    for x, v in enumerate(data):
        for y, value in enumerate(v):
            if abs(value)>casde_thresh:
                print(f"CASDE {value:.4f} > {casde_thresh} [GWP={x+1} | STEP={y} | TIME={y*timestep:.3f}]")
                if x not in prob_gwps:
                    prob_gwps.append(x)

    with open(os.path.join(basepath, 'maxf') , 'rb') as f:
        data = np.load(f)
    for x, v in enumerate(data):
        for y, value in enumerate(v):
            if abs(value)>maxf_thresh:
                print(f"MAXF {value:.4f} > {maxf_thresh} [GWP={x+1} | STEP={y} | TIME={y*timestep:.3f}]")
                if x not in prob_gwps:
                    prob_gwps.append(x)

    with open(os.path.join(basepath, 'nucde') , 'rb') as f:
        data = np.load(f)
    for x, v in enumerate(data):
        if abs(v)>nucde_thresh:
            print(f"NUCDE {v:.4f} > {nucde_thresh} [STEP={x} | TIME={x*timestep:.3f}]")

    print(f'GWPs need further attention : {[i+1 for i in prob_gwps]}')

def plotgwpforce(basepath, manifest):
    assert('maxf' in manifest['quantities'])
    print('Which GWPs to plot? Use a comma seperated list')
    gwps = input()
    gwps = [int(i)-1 for i in gwps.split(',')]
    raw_data = np.load(os.path.join(basepath, 'maxf'))
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for i in gwps:
        ax.plot(times, raw_data[i], label=f'GWP{i+1}')
    ax.set_title('Maximum force (per GWP)')
    ax.set_ylabel('Force (au)')
    ax.set_xlabel('Time (fs)')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(basepath, 'dbg_force.png'))
    plt.show()

def plotcascon(basepath, manifest):
    assert('casde' in manifest['quantities'])
    print('Which GWPs to plot? Use a comma seperated list')
    gwps = input()
    gwps = [int(i)-1 for i in gwps.split(',')]
    raw_data = np.load(os.path.join(basepath, 'casde'))
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for i in gwps:
        ax.plot(times, raw_data[i], label=f'GWP{i+1}')
    ax.set_title('CASSCF convergence (per GWP)')
    ax.set_ylabel('Convergence')
    ax.set_xlabel('Time (fs)')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(basepath, 'dbg_casde.png'))
    plt.show()

def plotnms(basepath, manifest):
    assert('nm' in manifest['quantities'])
    print('Which GWP to plot?')
    gwp = input()
    gwp = int(gwp)-1
    print('Which NMs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]

    raw_data = np.load(os.path.join(basepath, 'nm'))[gwp].T
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for i in nms:
        ax.plot(times, raw_data[i], label=f'NM{i+1}')
    ax.set_title(f'Normal mode evolution (for GWP{gwp+1})')
    ax.set_ylabel('Normal mode evolution')
    ax.set_xlabel('Time (fs)')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(basepath, f'dbg_nms_gwp{gwp+1}.png'))
    plt.show()

def plotcsf(basepath, manifest):
    assert('csf' in manifest['quantities'])
    print('Which GWP to plot?')
    gwp = input()
    gwp = int(gwp)-1
    print('Which CSFs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]

    raw_data = np.load(os.path.join(basepath, 'csf'))[gwp].T
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for i in nms:
        ax.plot(times, raw_data[i], label=f'NM{i+1}')
    ax.set_title(f'CSF population evolution (for GWP{gwp+1})')
    ax.set_ylabel('CSF population evolution')
    ax.set_xlabel('Time (fs)')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(basepath, f'dbg_csf_gwp{gwp+1}.png'))
    plt.show()

def plotci(basepath, manifest):
    assert('ci' in manifest['quantities'])
    print('Which GWP to plot?')
    gwp = input()
    gwp = int(gwp)-1
    print('Which CIs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]

    raw_data = np.load(os.path.join(basepath, 'ci'))[gwp].T
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
    for i in nms:
        ax.plot(times, raw_data[i], label=f'NM{i+1}')
    ax.set_title(f'CI population evolution (for GWP{gwp+1})')
    ax.set_ylabel('CI population evolution')
    ax.set_xlabel('Time (fs)')
    ax.legend(loc='upper right')
    fig.savefig(os.path.join(basepath, f'dbg_ci_gwp{gwp+1}.png'))
    plt.show()

def plotbl(basepath, manifest):
    assert('xyz' in manifest['quantities'])
    print('Which GWP to plot?')
    gwp = input()
    gwp = int(gwp)-1
    print('Which BLs to plot? Give a dash, comma seperated list e.g 1-2,2-3,4-6')
    bls = input()
    BPS = []
    for x in bls.split(','):
        a, b = [int(z) for z in x.split('-')]
        BPS.append([a,b])

    raw_data = np.load(os.path.join(basepath, 'xyz'))[gwp]
    times = np.load(os.path.join(basepath, 'times'))
    fig, ax = plt.subplots()
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
    fig.savefig(os.path.join(basepath, f'dbg_ci_gwp{gwp+1}.png'))
    plt.show()

if len(sys.argv) < 3:
    print(f'Use {sys.argv[0]} path/to/manifest.json task')
    print('''Availble tasks:
    qs = QuickSearch help identify problem GWPs
    pforce = Plot Gradient (Force)
    pcasde = Plot CASSCF Convergence [GWP-wise]
    pnm = Plot normal modes [GWP-wise]
    pcsf = Plot csf populations [GWP-wise]
    pci = Plot ci populations [GWP-wise]
    pbl = Plot bond lengths [GWP-wise]
    ''')
    sys.exit(-1)

_, manifestpath, task = sys.argv
basepath = os.path.dirname(manifestpath)
try:
    with open(manifestpath, 'r') as f:
        manifest = json.load(f)
except:
    print('Unable to open manifest!')
    sys.exit(-1)

task = task.lower()
if task=='qs' : find_probelm_gwps(basepath, manifest)
if task=='pforce' : plotgwpforce(basepath, manifest)
if task=='pcasde' : plotcascon(basepath, manifest)
if task=='pnm' : plotnms(basepath, manifest)
if task=='pcsf' : plotcsf(basepath, manifest)
if task=='pci' : plotci(basepath, manifest)
if task=='pbl' : plotbl(basepath, manifest)