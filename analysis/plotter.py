import mathutils
import numpy as np
import sys
import os
import json

if len(sys.argv) < 5:
    print(f"""Not enough arguemnts!\n Use : {sys.argv[0]}
    [/path/manifest.json] 

    *** GEOMETRIC ***
    => [BL=1-2,2-3]
    => [PBL=1-2,2-3]
    => [NM=1,2,3]
    => [BA=1-2-3]
    => [DA=1-2-3-4]

    *** ELECTRONIC ***
    => [CIs=A|1,2,3]
    => [CSFs=A|1,2,3]
    => [CSFv=label:1,1,0,0_label:0,0,1,1]
    => [AVCSFs=A|1,2,3:av_window]
    => [AVCSFv=av_window_label:1,0,0_label:0,0,1]
    => [SD=A|1,2]
    => [MQ=A|1,2]
    => [HM=sd|mq:atom_num:nbins:min_y:max_y]
    
    *** DISCRETE FOURIER TRANSFORMS ***
    => [FFT=cd|mq|csf(+|-|p):CHOP:1,2,3|A:[ZOOM_RANGE/None]]
    [width/panel, height]
    [Output x11 | filename.png]
    """)
    sys.exit(0)

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
    avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
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
        axes[n].legend(loc='upper right')

    elif cmd == 'PBL':
        avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
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
        axes[n].legend(loc='upper right')


    elif cmd == 'NM':
        nms = np.load(os.path.join(basepath, 'nm_ave'))
        for x in [int(i) for i in ins.split(',')]:
            axes[n].plot(times, nms[x-1], label=f'NM{x}', color=get_nth_col(x-1))
        axes[n].set_title('Normal mode evolution')
        axes[n].set_ylabel('Normal mode excitation')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'BA':
        avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
        BPS = []
        for x in ins.split(','):
            a, b, c = [int(z) for z in x.split('-')]
            BPS.append([a,b,c])

        for a in BPS:
            dp = []
            for x in range(nsteps):
                dp.append(mathutils.MathUtils.bond_angle(avegeom[x, a[0]-1],avegeom[x, a[1]-1], avegeom[x, a[2]-1], mode='deg'))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'
            try: alab3 = ATOMICLABELS[manifest['atomnos'][str(a[2])]-1]
            except: alab3 = '?'

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}]')
        axes[n].set_ylabel('Bond angle (deg)')
        axes[n].set_title('Bond angle')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'PBA':
        BPS = []
        for x in ins.split(','):
            a, b, c = [int(z) for z in x.split('-')]
            BPS.append([a,b,c])

        for a in BPS:
            dp = []
            for x in range(nsteps):
                dp.append(mathutils.MathUtils.bond_angle(avegeom[x, a[0]-1],avegeom[x, a[1]-1], avegeom[x, a[2]-1], mode='deg'))

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'
            try: alab3 = ATOMICLABELS[manifest['atomnos'][str(a[2])]-1]
            except: alab3 = '?'
            dp = np.array(dp)
            dp0 = dp[0]
            dpx = (dp-dp0)/dp0
            axes[n].plot(times, dpx, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}]')
        axes[n].set_title('Bond angles (fractional)')
        axes[n].set_ylabel('Fractional change')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'DA':
        avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
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
        axes[n].legend(loc='upper right')
        
    # ELECTRONICS
    elif cmd == 'CIs':
        adiabats = np.load(os.path.join(basepath, 'ci_ave'))
        CI_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(adiabats.shape[0]):
            if CI_STATES == None: pass
            else:
                if i+1 not in CI_STATES: continue
            axes[n].plot(times, adiabats[i], label=f'CI {i+1}', color=get_nth_col(i))
        axes[n].set_title('Adiabatic [CI] state evolution')
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'CSFs':
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        CSF_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(diabats.shape[0]):
            if CSF_STATES == None: pass
            else:
                if i+1 not in CSF_STATES: continue
            axes[n].plot(times, diabats[i], label=f'CSF {i+1}', color=get_nth_col(i))
        axes[n].set_title('Diabatic [CSF] state evolution')
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'CSFv':
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        # Expect a list of label:1,1,0,0_label:0,0,1,1
        to_plot={}
        for i in ins.split('_'):
            label, nums = i.split(':')
            nums = [float(j) for j in nums.split(',')]
            assert(len(nums)==diabats.shape[0])
            to_plot[label] = np.array(nums)
        for k, v in to_plot.items():
            data = np.dot(v, diabats)
            axes[n].plot(times, data, label=k)
        axes[n].set_title('Diabatic [CSF] state vector evolution')
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'AVCSFs':
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        csfs, av_window = ins.split(':')
        av_window = int(av_window)
        CSF_STATES = None if csfs=='A' else [int(i) for i in csfs.split(',')]
        for i in range(diabats.shape[0]):
            if CSF_STATES == None: pass
            else:
                if i+1 not in CSF_STATES: continue
            mav = mathutils.MathUtils.moving_avg(diabats[i], av_window)
            axes[n].plot(times[:len(mav)], mav, label=f'AVCSF {i+1}', color=get_nth_col(i))
        axes[n].set_title(f'{av_window} point moving average diabatic [CSF] state evolution')
        axes[n].set_ylabel('Averaged state population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'AVCSFv':
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        # Expect format av_window_label:1,0,0_label:0,0,1
        av_window, *toplt = ins.split('_')
        av_window = int(av_window)

        to_plot={}
        for i in toplt:
            print(i)
            label, nums = i.split(':')
            nums = [float(j) for j in nums.split(',')]
            assert(len(nums)==diabats.shape[0])
            to_plot[label] = np.array(nums)
        mavs = np.array([mathutils.MathUtils.moving_avg(i, av_window) for i in diabats])
        print(mavs.shape)
        for k, v in to_plot.items():
            data = np.dot(v, mavs)
            axes[n].plot(times[:len(data)], data, label=k)
        axes[n].set_title(f'{av_window} point moving average diabatic [CSF] state custom vector evolution')
        axes[n].set_ylabel('Averaged state population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'SD':
        sd = np.load(os.path.join(basepath, 'sd_ave'))
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
        axes[n].legend(loc='upper right')

    elif cmd == 'AVSD':
        csfs, av_window = ins.split(':')
        av_window = int(av_window)
        SDS = None if csfs=='A' else [int(i) for i in csfs.split(',')]
        for i in range(len(manifest['spindenmap'])):
            atom_number = manifest['spindenmap'][i]
            if SDS == None : pass
            elif atom_number not in SDS: continue
            
            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'
            mav = mathutils.MathUtils.moving_avg(sd[i], av_window)
            axes[n].plot(times[:len(mav)], mav, label='{} [{}]'.format(symbol, atom_number), color=get_nth_col(atom_number))
        axes[n].set_title(f'{av_window}-point moving average spin density (H Summed)')
        axes[n].set_ylabel('Spin density')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'MQ':
        mq = np.load(os.path.join(basepath, 'mq_ave'))
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
        axes[n].legend(loc='upper right')
        
    # Heatmaps - currently SD/MQ (may want to add BL)
    elif cmd == 'HM':
        def sbin_vals(nb, vals, val_scale, max_val, min_val):
            x = np.zeros(nb)
            step = (max_val-min_val)/nb
            for n, v in enumerate(vals):
                for i in range(nb):
                    if (min_val + i * step) < v < (min_val + (i+1) * step):
                        x[i] += val_scale[n]
            return x[::-1]

        def bin_vals(nb, vals, max_val, min_val):
            val_len = len(vals)
            return sbin_vals (nb, vals, np.repeat(1, val_len), max_val, min_val)

        mode, selector, nbins, minval, maxval = ins.split(':')
        nbins = int(nbins)
        minval = float(minval)
        maxval = float(maxval)
        bindata = np.zeros((len(times), nbins))

        if mode == 'sd':
            atom_number = int(selector)
            mapper = manifest['spindenmap'].index(atom_number)

            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'

            axes[n].set_title(f'Spin denisity (H Summed) heatmap for atom {symbol}[{atom_number}]')
            axes[n].set_xlabel('Spin density (H Summed)')

            ave_data = np.load(os.path.join(basepath, 'sd_ave'))[mapper]
            unbinned_data=np.load(os.path.join(basepath, 'sd'))[mapper]

        elif mode == 'mq':
            atom_number = int(selector)
            mapper = manifest['mullikenmap'].index(atom_number)

            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'

            axes[n].set_title(f'Mulliken charge (H Summed) heatmap for atom {symbol}[{atom_number}]')
            axes[n].set_xlabel('Mulliken charge (H Summed)')

            ave_data = np.load(os.path.join(basepath, 'mq_ave'))[mapper]
            unbinned_data=np.load(os.path.join(basepath, 'mq'))[mapper]

        else:
            raise Exception(f"Illegal mode {mode} for heatmap")

        for i, _ in enumerate(times):
            bindata[i] = bin_vals(nbins, unbinned_data[i], maxval, minval)

        axes[n].set_xlabel('Time (fs)')
        timewidth = (times[1]-times[0])/2 # Tiny fudging to get heatmap to align nicely
        axes[n].imshow(bindata.T, cmap='inferno', extent=(-timewidth, times[-1]+timewidth, minval, maxval), aspect='auto')
        axes[n].plot(times, ave_data, color='white', linestyle='dotted',linewidth=2)

    # FFT
    elif cmd == 'FFT':
        # [FFT=cd|mq|csf(+ Do Phase Correction | - Do cos correction | 2 Power spectra):Cutoff:1,2,3|A:limxmin,limxmax]
        mode, CHOP, selector, lims = ins.split(':')
        CHOP=int(CHOP)
        if lims == 'None' : lims = None
        else: lims = [float(x) * 1E15 for x in lims.split(',')]
        selector = None if selector == 'A' else [int(i) for i in selector.split(',')]
        
        doPhase=False
        # Phase modulation is gone based on the sign of FT value
        if '+' in mode : 
            doPhase=True
            mode=mode.replace('+', '')

        doCosCorrect=False
        if '-' in mode : 
            doCosCorrect=True
            mode=mode.replace('-', '')

        doPowerSpectra=False
        if '2' in mode : 
            doPowerSpectra=True
            mode=mode.replace('2', '')

        # Do a phase plot
        phase_mode =False
        if 'p' in mode : 
            phase_mode=True
            mode=mode.replace('p', '')

        if mode == 'csf':  data = np.load(os.path.join(basepath, 'csf_ave'))
        elif mode == 'mq': data = np.load(os.path.join(basepath, 'mq_ave'))
        elif mode == 'sd': data = np.load(os.path.join(basepath, 'sd_ave'))
        else: raise Exception('Illegal FFT mode')

        # print(data.shape, len(times))
        assert(data.shape[1] == len(times)) # Make sure extract worked

        data_fft  = data

        if doCosCorrect:
            N = data_fft.shape[1] * 2
        else:
            N = data_fft.shape[1]

        freq = np.fft.fftfreq(N, d=times[1]*1E-15-times[0]*1E-15)[CHOP:int(N/2)]
        fftdata = {}
        colourdata = {}
        for i in range(data_fft.shape[0]):
            if mode=='csf': # CSFs are picked by index
                if selector == None : pass
                elif i+1 not in selector: continue
                label = f'CSF {i+1}'
                colour = get_nth_col(i)
            
            elif mode == 'sd' or mode == 'mq': # SD/MQ are picked based on atom number
                if mode == 'sd' : atom_number = manifest['spindenmap'][i]
                else : atom_number = manifest['mullikenmap'][i]
                
                if selector == None : pass
                elif atom_number not in selector: continue

                try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
                except: symbol = '?'
                label = f'{symbol}[{atom_number}]'
                colour = get_nth_col(atom_number)
            else:
                raise Exception("Illegal FFT Mode")

            data = data_fft[i]
            if doCosCorrect:
                scalearray = [np.cos(np.pi * j / (2 * times[-1])) for j in times]
                data = data * scalearray
                data = np.concatenate((data[::-1], data))
            ft = np.fft.fft(data)
            fftdata[label] = ft
            colourdata[label] = colour

        if phase_mode:
            from itertools import combinations
            for i, j in combinations(fftdata.keys(), 2):
                phases = np.abs(np.angle(fftdata[i] / fftdata[j]))
                axes[n].plot(freq, phases[CHOP:int(N/2)], label=f'{i} - {j}')
            axes[n].set_title(f'FFT Phases {mode}')
            axes[n].set_ylabel('Relative phases (rad)')
            axes[n].set_xlabel('Frequency Hz')
            axes[n].axhline(y=np.pi, color="black", linestyle="--")

        else:
            for k, v in fftdata.items():
                magnitude = np.abs(v)
                if doPowerSpectra:
                    magnitude = magnitude * magnitude
                if doPhase:
                    phase = np.angle(v)
                    combined = magnitude * -np.sign(phase) # Minus as most freq will be -ve => visual clarity
                else:
                    combined = magnitude
                axes[n].plot(freq, combined[CHOP:int(N/2)], label=k, color=colourdata[k])
            title = f'[{mode}]'
            if doPhase : title = 'Phase Modulated ' + title

            if doPowerSpectra: title = 'Power Spectra '+ title
            else: title = title + 'FT Spectra '+ title

            axes[n].set_title(title)
            axes[n].set_ylabel('Intensity')
            axes[n].set_xlabel('Frequency Hz')
        axes[n].legend(loc='upper right')
        axes[n].axhline(y=0, color="black", linestyle="--")

        if lims != None:
            axes[n].axis(xmin = lims[0], xmax = lims[1])

    else:
        raise Exception(f'Illegal mode {cmd}')
fig.tight_layout() 
if OUTPUT=='x11' : plt.show()
else: plt.savefig(OUTPUT, dpi=500)
