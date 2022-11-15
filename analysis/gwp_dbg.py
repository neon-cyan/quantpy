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

def get_nth_col(idx):
    cols =['#e6194B', '#3cb44b', '#FFC800', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#800000', '#aaffc3', '#000075', '#a9a9a9']
    return cols[idx%len(cols)]

# General probelm-GWP-finding function

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

# Plots all GWPs properties

def plotgwpforce(basepath, manifest):
    assert('maxf' in manifest['quantities'])
    print('Which GWPs to plot? Use a comma seperated list')
    raw_data = np.load(os.path.join(basepath, 'maxf'))
    gwps = input()
    if gwps == '*': gwps=list(range(0,raw_data.shape[0]))
    else: gwps = [int(i)-1 for i in gwps.split(',')]

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
    raw_data = np.load(os.path.join(basepath, 'casde'))
    gwps = input()
    if gwps == '*': gwps=list(range(0,raw_data.shape[0]))
    else: gwps = [int(i)-1 for i in gwps.split(',')]
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

# Singlewise GWP plots

def plotnms(basepath, manifest):
    assert('nm' in manifest['quantities'])
    print('Which GWP to plot? * for all')
    raw_data = np.load(os.path.join(basepath, 'nm'))
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which NMs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]

    for g in gwp:
        rd = raw_data[g].T
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for i in nms:
            ax.plot(times, rd[i], label=f'NM{i+1}')
        ax.set_title(f'Normal mode evolution (for GWP{g+1})')
        ax.set_ylabel('Normal mode evolution')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_nms_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotcsf(basepath, manifest):
    assert('csf' in manifest['quantities'])
    print('Which GWP to plot? * for all')
    raw_data = np.load(os.path.join(basepath, 'csf'))
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which CSFs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]
    for g in gwp:
        rd = raw_data[g].T
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for i in nms:
            ax.plot(times, np.square(np.abs(rd[i])), label=f'CSF{i+1}')
        ax.set_title(f'CSF population evolution (for GWP{g+1})')
        ax.set_ylabel('CSF population evolution')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_csf_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotci(basepath, manifest):
    assert('ci' in manifest['quantities'])
    print('Which GWP to plot?')
    raw_data = np.load(os.path.join(basepath, 'ci'))
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which CIs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]
    for g in gwp:
        rd = raw_data[g].T
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for i in nms:
            ax.plot(times, np.square(np.abs(rd[i])), label=f'CI{i+1}', color=get_nth_col(i))
        ax.set_title(f'CI population evolution (for GWP{g+1})')
        ax.set_ylabel('CI population evolution')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_ci_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotcsfv(basepath, manifest):
    assert('csf' in manifest['quantities'])
    print('Which GWP to plot? * for all')
    raw_data = np.load(os.path.join(basepath, 'csf'))
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which CSF vectors and labels? label1:0,0,1_label2:1,0,0')
    nms = input()
    to_plot={}
    for i in nms.split('_'):
        label, nums = i.split(':')
        nums = [float(j) for j in nums.split(',')]
        # assert(len(nums)==diabats.shape[1])
        vect = np.array(nums)
        to_plot[label] = vect / np.linalg.norm(vect)
        # Normalize the vector

    for g in gwp:
        rd = raw_data[g]
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for k,v in to_plot.items():
            ax.plot(times, np.square(np.abs(rd.dot(v))), label=k)
        ax.set_title(f'CSF/SD vector evolution (for GWP{g+1})')
        ax.set_ylabel('CSF/SD vector population')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_csf_vect_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotpes(basepath, manifest):
    assert('ci' in manifest['quantities'])
    print('Which GWP to plot?')
    raw_data = np.load(os.path.join(basepath, 'cies'))
    # print(raw_data.shape)
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which CIs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]
    for g in gwp:
        rd = raw_data[g].T
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for i in nms:
            ax.plot(times, rd[i], label=f'CI{i+1}', color=get_nth_col(i))
        ax.set_title(f'TD-PES (for GWP{g+1})')
        ax.set_ylabel('Energy / Ha')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_pes_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotbl(basepath, manifest):
    assert('xyz' in manifest['quantities'])
    print('Which GWP to plot?')
    raw_data = np.load(os.path.join(basepath, 'xyz'))
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which BLs to plot? Give a dash, comma seperated list e.g 1-2,2-3,4-6')
    bls = input()
    BPS = []
    for x in bls.split(','):
        a, b = [int(z) for z in x.split('-')]
        BPS.append([a,b])
    for g in gwp:
        rd = raw_data[g]
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()
        for a in BPS:
            dp = []
            for x in range(manifest['steps']):
                dp.append(mathutils.MathUtils.bond_length(rd[x, a[0]-1],rd[x, a[1]-1] ))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'

            ax.plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
        ax.set_title(f'BL evolution (for GWP{g+1})')
        ax.set_ylabel('Bond length (Å)')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.tight_layout()
        fig.savefig(os.path.join(basepath, f'dbg_bl_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotfnm(basepath, manifest):
    assert('nm' in manifest['quantities'])
    assert('forces' in manifest['quantities'])

    print('Which GWP to plot?')
    raw_data = np.load(os.path.join(basepath, 'forces'))[gwp]
    gwp = input()
    if gwp == '*': gwp = list(range(0, raw_data.shape[0]))
    else: gwp = [int(gwp)-1]

    print('Which NMs to plot? Give a comma seperated list')
    nms = input()
    nms = [int(i)-1 for i in nms.split(',')]
    for g in gwp:
        rd = raw_data[g]
        xyz2nm = np.load(os.path.join(basepath, 'xyz2nm'))
        times = np.load(os.path.join(basepath, 'times'))
        fig, ax = plt.subplots()

        for a in nms:
            dp = []
            xyz_nma = xyz2nm[a]
            for b in range(rd.shape[0]):
                dp.append(xyz_nma.dot(rd[b].reshape(xyz_nma.shape)))
            ax.plot(times, dp, label=f'NM{a+1}')
        ax.set_title(f'Gradient in terms of normal modes (for GWP{g+1})')
        ax.set_ylabel('Gradient as normal mode')
        ax.set_xlabel('Time (fs)')
        ax.legend(loc='upper right')
        fig.savefig(os.path.join(basepath, f'dbg_pfnm_gwp{g+1}.png'))
        if len(gwp) == 1 : plt.show()
    print('Plot OK')

def plotl118e(basepath, manifest):
    assert('nucde' in manifest['quantities'])
    assert(manifest['method'] == 'l118')
    fig, axes = plt.subplots(1, 3, figsize=(21,5))

    times = np.load(os.path.join(basepath, 'times'))
    nucde = np.load(os.path.join(basepath, 'nucde'))
    nucke = np.load(os.path.join(basepath, 'nucke'))
    nucpe = np.load(os.path.join(basepath, 'nucpe'))

    axes[0].plot(times, nucde)
    axes[0].set_title('Total energy')
    axes[0].set_ylabel('Energy / Ha')
    axes[0].set_xlabel('Time (fs)')

    axes[1].plot(times, nucke)
    axes[1].set_title('Kinetic energy')
    axes[1].set_ylabel('Energy / Ha')
    axes[1].set_xlabel('Time (fs)')

    axes[2].plot(times, nucpe)
    axes[2].set_title('Potential energy')
    axes[2].set_ylabel('Energy / Ha')
    axes[2].set_xlabel('Time (fs)')
    fig.tight_layout()
    fig.savefig(os.path.join(basepath, f'L118_Energies.png'))
    print('Plot OK')

if len(sys.argv) < 3:
    print(f'Use {sys.argv[0]} path/to/manifest.json task')
    print('''Availble tasks:
    qs = QuickSearch (helps identify problem GWPs)
    pforce = Plot Max Gradient (Force)
    pcasde = Plot CASSCF Convergence [GWP-wise]
    pnm = Plot normal modes [GWP-wise]
    pcsf = Plot csf populations [GWP-wise]
    pci = Plot ci populations [GWP-wise]
    pcsfv = Plot ci vector [GWP-wise]
    pbl = Plot bond lengths [GWP-wise]
    pfnm = Plot gradient in normal modes
    ppes = Plot TD-PES [GWP-wise]
    pl118e = Plot the L118 nucelar energies (EKin / EPot / ETot)
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
elif task=='pforce' : plotgwpforce(basepath, manifest)
elif task=='pcasde' : plotcascon(basepath, manifest)
elif task=='pnm' : plotnms(basepath, manifest)
elif task=='pcsf' : plotcsf(basepath, manifest)
elif task=='pci' : plotci(basepath, manifest)
elif task=='pcsfv' : plotcsfv (basepath, manifest)
elif task=='pbl' : plotbl(basepath, manifest)
elif task=='pfnm' : plotfnm(basepath, manifest)
elif task=='ppes' : plotpes(basepath, manifest)
elif task=='pl118e' : plotl118e(basepath, manifest)

else: print('Invalid job!')