import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os
import mathutils
from defs import ATOMICLABELS
# Welcome to jank city - things may or may not work


def get_nth_col(idx):
    cols =['#e6194B', '#3cb44b', '#FFC800', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#800000', '#aaffc3', '#000075', '#a9a9a9']
    return cols[idx%len(cols)]

# General probelm-GWP-finding function

def find_problem_gwps(basepath, manifest):
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
# Force + CASCON plotting is depracated
# NM, BL, CI + CSF(V) is deprecated - use plotter with GWP <- notation

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
if task=='qs' : find_problem_gwps(basepath, manifest)
elif task=='ppes' : plotpes(basepath, manifest)
elif task=='pl118e' : plotl118e(basepath, manifest)
else: print('Invalid job!')