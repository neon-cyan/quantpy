import mathutils
import numpy as np
import sys
import os
import json

if len(sys.argv) < 5:
    print(f"""Not enough arguemnts!\n Use : {sys.argv[0]}
    [/path/manifest.json] 
    => [BL=1-2,2-3]
    => [PBL=1-2,2-3]
    => [NM=1,2,3]
    => [CIs=A|1,2,3]
    => [CSFs=A|1,2,3]
    => [SD=A|1,2]
    => [MQ=A|1,2]
    => [FFT=cd|mq|csf:CHOP:[START:END]|A:1,2,3|A]
    [width/panel, height]
    [Output x11 | filename.png]
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
commands = sys.argv[2:-2]
wpp, height = [float(i) for i in sys.argv[-2].split(',')]
OUTPUT = sys.argv[-1] # Either plot to a file or x11 window

assert(os.path.exists(manifest_path))
with open(manifest_path, 'r') as f:
    manifest = json.load(f)
basepath = os.path.dirname(manifest_path)

times = np.load(os.path.join(basepath, 'times'))
nms = np.load(os.path.join(basepath, 'nm_ave'))
diabats = np.load(os.path.join(basepath, 'csf_ave'))
adiabats = np.load(os.path.join(basepath, 'ci_ave'))
avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
mq = np.load(os.path.join(basepath, 'mq_ave'))
sd = np.load(os.path.join(basepath, 'sd_ave'))
nsteps = manifest['steps']


# DO PLOTTING

import matplotlib.pyplot as plt

# Define custom default colours - keep consistent between plots
# List lossely based on https://sashamaps.net/docs/resources/20-colors/
# Do this for the CSFs/CIs/FFT/NMs/MQ/SD
def get_nth_col(idx):
    cols =['#e6194B', '#3cb44b', '#FFC800', '#4363d8', '#f58231', '#42d4f4', '#f032e6', '#fabed4', '#469990', '#dcbeff', '#9A6324', '#800000', '#aaffc3', '#000075', '#a9a9a9']
    return cols[idx%len(cols)]

fig, axes = plt.subplots(1, len(commands), num=manifest_path, figsize=(wpp * len(commands), height))
if len(commands) == 1 : axes = [axes] # MPL messes with array if only one plot => Need to re-array

for n, c in enumerate(commands):
    cmd, ins = c.split('=')
    # GEOMETRICS
    if cmd == 'BL':
        BPS = []
        for x in ins.split(','):
            a, b = [int(z) for z in x.split('-')]
            BPS.append([a,b])
        for a in BPS:
            dp = []
            for x in range(nsteps):
                dp.append(mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] ))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
        axes[n].set_title('Bond lengths')
        axes[n].set_ylabel('Bond length (Å)')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'PBL':
        BPS = []
        for x in ins.split(','):
            a, b = [int(z) for z in x.split('-')]
            BPS.append([a,b])
        for a in BPS:
            dp = []
            init_bl = mathutils.MathUtils.bond_length(avegeom[0, a[0]-1],avegeom[0, a[1]-1] )
            for x in range(nsteps):
                bl = mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] )
                dp.append((bl - init_bl) / init_bl)

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
        axes[n].set_title('Bond lengths (fractional)')
        axes[n].set_ylabel('Fractional change')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');


    elif cmd == 'NM':
        for x in [int(i) for i in ins.split(',')]:
            axes[n].plot(times, nms[x-1], label=f'NM{x}', color=get_nth_col(x-1))
        axes[n].set_title('Normal mode evolution')
        axes[n].set_ylabel('Normal mode excitation')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'BA':
        BPS = []
        for x in ins.split(','):
            a, b, c = [int(z) for z in x.split('-')]
            BPS.append([a,b,c])

        for a in BPS:
            dp = []
            for x in range(nsteps):
                dp.append(mathutils.MathUtils.bond_angle(avegeom[x, a[0]-1],avegeom[x, a[1]-1], avegeom[x, a[2]-1] ))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'
            try: alab3 = ATOMICLABELS[manifest['atomnos'][str(a[2])]-1]
            except: alab3 = '?'

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}]')
        axes[n].set_ylabel('Bond angle (rad)')
        axes[n].set_title('Bond angle')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'DA':
        BPS = []
        for x in ins.split(','):
            a, b, c, d = [int(z) for z in x.split('-')]
            BPS.append([a,b,c,d])

        for a in BPS:
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

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}] - {alab4}[{a[3]}]')
        axes[n].set_ylabel('Dihedral angle (rad)')
        axes[n].set_title('Bond angle')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');
        
    # ELECTRONICS
    elif cmd == 'CIs':
        CI_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(adiabats.shape[0]):
            if CI_STATES == None: pass
            else:
                if i+1 not in CI_STATES: continue
            axes[n].plot(times, adiabats[i], label=f'CI {i+1}', color=get_nth_col(i))
        axes[n].set_title('Adiabatic [CI] state evolution')
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'CSFs':
        CSF_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(diabats.shape[0]):
            if CSF_STATES == None: pass
            else:
                if i+1 not in CSF_STATES: continue
            axes[n].plot(times, diabats[i], label=f'CSF {i+1}', color=get_nth_col(i))
        axes[n].set_title('Diabatic [CSF] state evolution')
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'SD':
        SDS = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(len(manifest['spindenmap'])):
            atom_number = manifest['spindenmap'][i]
            if SDS == None : pass
            elif atom_number not in SDS: continue
            
            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'
            axes[n].plot(times, sd[i], label='{} [{}]'.format(symbol, atom_number), color=get_nth_col(atom_number))
        axes[n].set_title('Spin density evolution (H Summed)')
        axes[n].set_ylabel('Spin density')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');

    elif cmd == 'MQ':
        MQS = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(len(manifest['mullikenmap'])):
            atom_number = manifest['mullikenmap'][i]
            if MQS == None : pass
            elif atom_number not in MQS: continue
            
            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'
            axes[n].plot(times, mq[i], label='{} [{}]'.format(symbol, atom_number),  color=get_nth_col(atom_number))
        axes[n].set_title('Mulliken charge evolution (H Summed)')
        axes[n].set_ylabel('Mulliken charge')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right');
    
    # FFT
    elif cmd == 'FFT':
        # [FFT=cd|mq|csf:CHOP:[START-END]|A]
        mode, CHOP, RANGE, selector = ins.split(':')
        CHOP=int(CHOP)
        selector = None if selector == 'A' else [int(i) for i in selector.split(',')]

        if mode == 'csf':  data = diabats
        elif mode == 'mq': data = mq
        elif mode == 'sd': data = sd
        else: raise Exception('Illegal FFT mode')

        print(data.shape, len(times))
        assert(data.shape[1] == len(times)) # Make sure extract worked

        if RANGE != 'A':
            s_idx, e_idx = [int(i) for i in RANGE.split('-')]
            times_fft = times[s_idx:e_idx]
            data_fft  = data.T[s_idx:e_idx].T
        else:
            times_fft = times
            data_fft  = data

        # Do FFT and plot up
        N = data_fft.shape[1]
        fig =  plt.figure(num=manifest_path, figsize=(20.0, 15.0))

        for i in range(data_fft.shape[0]):
            if mode=='csf': # CSFs are picked by index
                if selector == None : pass
                elif i+1 not in selector: continue

            ft = np.fft.fft(data_fft[i])
            ft = ft.real**2 + ft.imag**2
            freq = np.fft.fftfreq(N, d=times_fft[1]-times_fft[0])

            if mode == 'sd' or mode == 'mq': # SD/MQ are picked based on atom number
                if mode == 'sd' : atom_number = manifest['spindenmap'][i]
                else : atom_number = manifest['mullikenmap'][i]
                
                if selector == None : pass
                elif atom_number not in selector: continue

                try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
                except: symbol = '?'
                label = f'{symbol}[{atom_number}]'
                colour = get_nth_col(atom_number)
            else: 
                label = f'CSF {i+1}'
                colour = get_nth_col(i)

            axes[n].plot(freq[CHOP:int(N/2)], ft[CHOP:int(N/2)], label=label, color=colour)
        axes[n].set_title(f'FFT {mode}')
        axes[n].set_ylabel('Intensity')
        axes[n].set_xlabel('Frequency PHz')
        axes[n].legend(loc='upper right');
    

    else:
        raise Exception(f'Illegal mode {cmd}')
fig.tight_layout() 
if OUTPUT=='x11' : plt.show()
else: plt.savefig(OUTPUT, dpi=500)
