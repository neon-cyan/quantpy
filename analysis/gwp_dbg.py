import numpy as np
import matplotlib.pyplot as plt
import sys
import json
import os

# Welcome to jank city - things may or may not work
# In code gwps are 0-indexed to align with arrays - but users use 1-indexed notation

def find_probelm_gwps(basepath, manifest):
    assert('casde' in manifest['quantities'])
    assert('nucde' in manifest['quantities'])
    assert('maxf' in manifest['quantities'])
    timestep=manifest['tout']

    DEFAULT_CAST = 0.1
    DEFAULT_NUCT = 2.5
    DEFAULT_FORCET = 5.0

    print(f"Please input CASSCF Delta-E thereshold (default = {DEFAULT_CAST})")
    casde_thresh=input()
    casde_thresh = DEFAULT_CAST if casde_thresh=='' else float(casde_thresh)

    print(f"Please input QuEh hessian thereshold (default = {DEFAULT_FORCET})")
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
            if value>casde_thresh:
                print(f"CASDE {value:.4f} > {casde_thresh} [GWP={x+1} | STEP={y} | TIME={y*timestep:.3f}]")
                if x not in prob_gwps:
                    prob_gwps.append(x)

    with open(os.path.join(basepath, 'maxf') , 'rb') as f:
        data = np.load(f)
    for x, v in enumerate(data):
        for y, value in enumerate(v):
            if value>maxf_thresh:
                print(f"MAXF {value:.4f} > {maxf_thresh} [GWP={x+1} | STEP={y} | TIME={y*timestep:.3f}]")
                if x not in prob_gwps:
                    prob_gwps.append(x)

    with open(os.path.join(basepath, 'nucde') , 'rb') as f:
        data = np.load(f)
    for x, v in enumerate(data):
        if v>nucde_thresh:
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
    ax.legend(loc='upper right');
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
    ax.legend(loc='upper right');
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
    ax.legend(loc='upper right');
    fig.savefig(os.path.join(basepath, f'dbg_nms_gwp{gwp}.png'))
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
    ax.legend(loc='upper right');
    fig.savefig(os.path.join(basepath, f'dbg_csf_gwp{gwp}.png'))
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
    ax.legend(loc='upper right');
    fig.savefig(os.path.join(basepath, f'dbg_ci_gwp{gwp}.png'))
    plt.show()

if len(sys.argv) < 2:
    print(f'Use {sys.argv[0]} path/to/manifest.json')
    sys.exit(-1)

manifestpath=sys.argv[1]
basepath = os.path.dirname(manifestpath)
try:
    with open(manifestpath, 'r') as f:
        manifest = json.load(f)
except:
    print('Unable to open manifest!')
    sys.exit(-1)

print("""
#####################################
#       All in one debugging        #
#####################################

Options:
QS : Quicksearch for problem GWPs based on MaxForce / CAS DE / Quantics DelE
PFORCE : Plot MaxForce (Hessain)
PCASDE : Plot CASSCF Convergence
PNM : Plot normal mode evolutions
PCSF : Plot CSF populations
PCI : Plot CI populations

""")
command = input()
command = command.lower()
if command=='qs' : find_probelm_gwps(basepath, manifest)
if command=='pforce' : plotgwpforce(basepath, manifest)
if command=='pcasde' : plotcascon(basepath, manifest)
if command=='pnm' : plotnms(basepath, manifest)
if command=='pcsf' : plotcsf(basepath, manifest)
if command=='pci' : plotci(basepath, manifest)
