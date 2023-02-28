import mathutils
import numpy as np
import sys
import os
import json
from defs import ATOMICLABELS, get_nth_col
import matplotlib.pyplot as plt

if len(sys.argv) < 5:
    print(f"""
    Not enough arguments!
    Use : {sys.argv[0]} [/path/manifest.json] [GRAPHS] [width/panel, height] [Output x11 | filename.png]

    Where graphs can be a space seperated list of any of

    ** BOND LENGTHS **
    => BL=[1]-[2],[3]-[4]                           - At least one dash seperated atom number pair <*>
    => PBL=[1]-[2],[3]-[4]                          - At least one dash seperated atom number pair <*>

    ** BOND ANGLES **
    => BA=[1]-[2]-[3],[2]-[3]-[4]                   - At least one dash seperated atom number trio
    => PBA=[1]-[2]-[3],[2]-[3]-[4]                  - At least one dash seperated atom number trio

    ** DIHEDRAL ANGLE **
    => DA=[A]-[B]-[C]-[D],[A]-[B]-[C]-[E]           - At least one dash seperated atom number quartet

    ** VIBRATIONAL NORMAL MODE **
    => NM=[A],[B],[C]                               - Comma seperated list of at least one NM <*>

    ** CI POPULATION **
    => CIs=[A|1,2,3]                                - Comma seperated list of at least one CI number or A for all <*>

    ** CSF/SD POPULATION **
    => CSFs=[A|1,2,3]                               - Comma seperated list of at least one CSF/SD number or A for all <*>

    ** MODIFIED CSF/SD POPULATION **
    => CSFv=[label:1,1,0,0_label:0,0,1,1]           - At least one label and real vector of length CSF (Auto-rescaled) <*>
    => SUMSCFV=[label:1,0,0+0,1,0_label:0,0,1]      - At least one label and real vectors of length CSF (Each is auto-rescaled) <*>
    => AVCSF=[A|1,2,3:av_window]                    - Comma seperated list of at least one CSF/SD number or A for all and average window size <*>
    => CSFvReIm=[label:1,1,0,0]                     - Real + Imaginary parts for a label and CSF vector <*>

    ** SPIN DENSITY AND MULLIKEN CHARGE (+ HEAT MAPS) **
    => SD=[A|1,2]                                   - At least one atom number or A for all
    => MQ=[A|1,2]                                   - At least one atom number or A for all
    => HM=[sd|mq:atom_num:nbins:min_y:max_y]        - Either spin density XOR mulliken charge, comma seperated lsit of atoms (or A for all)
    
    ** FOURIER TRANSFORMS **
    => FFT=[sd|mq|csf][+/-/2/p]:CHOP:SELECTOR:RANGE
    ===> + will do phase modulation (up/down)
    ===> - will do cosine zero-to-zero correction
    ===> 2 will square (power spectra)
    ===> p will do a pair phase plot
    ===> The CHOP will discard first CHOP FFT datapoints (close to zero freq)
    ===> The selector refers to either atoms or CSFs comma seperated list or A for all
    ===> Range is either None or given as min,max (auto * 1E15)

    Commands marked with <*> support some level of per-GWP functionality
    This is invoked using the arrow operator <- e.g. BL<-1=1-2,2-3
    NB: You many need to use a backtick prefix in your bash shell : \<-
    A few special GWP-only commands are included for completeness
    => maxf<-GWP=                                   - Max force (AU) - takes no ordinary parameters
    => casde<-GWP=                                  - CASSCF convergence for GWPs - takes no ordinary parameters
    => fnm<-GWP=[1,2,3]                             - Forces expressed in normal mode coordinates
    """)
    sys.exit(0)

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

fig, axes = plt.subplots(1, len(commands), num=manifest_path, figsize=(wpp * len(commands), height))
if len(commands) == 1 : axes = [axes] # MPL messes with array if only one plot => Need to re-array

for n, c in enumerate(commands):
    cmd, ins = c.split('=')

    cmdptx = None
    if '<-' in cmd:
        cmd, cmdptx = cmd.split('<-') # Deal with GWPs ('extra context')
        cmdptx = [int(i)-1 for i in cmdptx.split(',')]
    cmd=cmd.lower()

    # GENERAL GWP STUFF
    if cmd in ['maxf', 'mf']:
        assert(len(cmdptx) >= 1)
        forces = np.load(os.path.join(basepath, 'maxf'))
        # print(forces.shape)
        for gwp in cmdptx:
            data = forces[gwp]
            axes[n].plot(times, data, label=f'GWP{gwp+1}')
        axes[n].set_title('Max forces (per GWP)')
        axes[n].set_ylabel('Max force (au)')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd in ['casde', 'cascon']:
        assert('casde' in manifest['quantities'])
        assert(len(cmdptx) > 1)
        raw_data = np.load(os.path.join(basepath, 'casde'))
        for i in cmdptx:
            axes[n].plot(times, raw_data[i], label=f'GWP{i+1}')
        axes[n].set_title('CASSCF convergence (per GWP)')
        axes[n].set_ylabel('Convergence')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd=='fnm':
        assert('nm' in manifest['quantities'])
        assert('forces' in manifest['quantities'])
        xyz2nm = np.load(os.path.join(basepath, 'xyz2nm'))
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for FNM plots!")
            else: raw_data = np.load(os.path.join(basepath, 'forces'))[cmdptx[0]]
        nms = [int(i)-1 for i in ins.split(',')]
        # print(nms)
        for a in nms:
            dp = []
            xyz_nma = xyz2nm[a]
            for b in range(raw_data.shape[0]):
                dp.append(xyz_nma.dot(raw_data[b].reshape(xyz_nma.shape)))
            axes[n].plot(times, dp, label=f'NM{a+1}')
        axes[n].set_title(f'Gradient in normal modes (for GWP{cmdptx[0]+1})')
        axes[n].set_ylabel('Gradient as normal mode')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')
        
    # GEOMETRICS
    elif cmd in ['bl', 'pbl']:
        avegeom = np.load(os.path.join(basepath, 'xyz_ave'))
        # Code to deal with GWPs
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for [P]BL plots!")
            else: avegeom = np.load(os.path.join(basepath, 'xyz'))[cmdptx[0]]
            title = title + f' (for GWP {cmdptx[0]+1})'
        BPS = []
        for x in ins.split(','):
            a, b = [int(z) for z in x.split('-')]
            BPS.append([a,b])
        for a in BPS:
            dp = []
            init_bl = mathutils.MathUtils.bond_length(avegeom[0, a[0]-1],avegeom[0, a[1]-1] )
            for x in range(nsteps):
                bl = mathutils.MathUtils.bond_length(avegeom[x, a[0]-1],avegeom[x, a[1]-1] )
                dp.append((bl - init_bl) / init_bl if cmd=='pbl' else bl)

            try: alab1 = ATOMICLABELS[manifest['atomnos'][str(a[0])]-1]
            except: alab1 = '?'
            try: alab2 = ATOMICLABELS[manifest['atomnos'][str(a[1])]-1]
            except: alab2 = '?'

            axes[n].plot(times, dp, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}]')
        
        if cmd=='pbl':
            title = 'Bond lengths (fractional)' + title
            axes[n].set_ylabel('Fractional change')
        else:
            title = 'Bond lengths' + title
            axes[n].set_ylabel('Bond length (Å)')
        axes[n].set_title(title)
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'nm':
        nms = np.load(os.path.join(basepath, 'nm_ave'))
        # Code to deal with GWPs
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for NM plots!")
            else: nms = np.load(os.path.join(basepath, 'nm'))[cmdptx[0]].T
            title = title + f' (for GWP {cmdptx[0]+1})'

        for x in [int(i) for i in ins.split(',')]:
            axes[n].plot(times, nms[x-1], label=f'NM{x}', color=get_nth_col(x-1))
        axes[n].set_title('Normal mode evolution'+title)
        axes[n].set_ylabel('Normal mode excitation')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'ba':
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

    elif cmd == 'pba':
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
            dp = np.array(dp)
            dp0 = dp[0]
            dpx = (dp-dp0)/dp0
            axes[n].plot(times, dpx, label=f'{alab1}[{a[0]}] - {alab2}[{a[1]}] - {alab3}[{a[2]}]')
        axes[n].set_title('Bond angles (fractional)')
        axes[n].set_ylabel('Fractional change')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'da':
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
    elif cmd == 'ci' or cmd == 'cis':
        adiabats = np.load(os.path.join(basepath, 'ci_ave'))
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for CI plots!")
            adiabats = np.load(os.path.join(basepath, 'ci'))[cmdptx[0]].T
            adiabats=np.square(np.abs(adiabats))
            title = title + f' (for GWP {cmdptx[0]+1})'

        CI_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(adiabats.shape[0]):
            if CI_STATES == None: pass
            else:
                if i+1 not in CI_STATES: continue
            axes[n].plot(times, adiabats[i], label=f'CI {i+1}', color=get_nth_col(i))
        axes[n].set_title('Adiabatic [CI] state evolution'+title)
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd in ['csfs','csf','sdet','sdets']:
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for CSF/SD plots!")
            diabats = np.load(os.path.join(basepath, 'csf'))[cmdptx[0]].T
            diabats=np.square(np.abs(diabats))
            title = title + f' (for GWP {cmdptx[0]+1})'
        CSF_STATES = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(diabats.shape[0]):
            if CSF_STATES == None: pass
            else:
                if i+1 not in CSF_STATES: continue
            axes[n].plot(times, diabats[i], label=f'CSF {i+1}', color=get_nth_col(i))
        axes[n].set_title('Diabatic [CSF] state evolution'+title)
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd in ['sumcsfv','csfv']:
        # expect input of form label:0,0,1+1,0,0_label:1,1,1
        diabats = np.load(os.path.join(basepath, 'zcsf'))
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for CSF/SD Vector plots!")
            diabats = np.load(os.path.join(basepath, 'csf'))[cmdptx[0]].T
            title = title + f' (for GWP {cmdptx[0]+1})'
        # print(diabats.shape)
        # Expect a list of label:1,1,0,0_label:0,0,1,1
        to_plot={}
        for i in ins.split('_'):
            label, nums = i.split(':')
            sumcomps = nums.split('+')
            to_plot[label] = []
            for s in sumcomps:
                vect = np.array([float(j) for j in s.split(',')])
                vect = vect / np.linalg.norm(vect)
                # print(label, vect)
                to_plot[label].append(vect)
        for k, v in to_plot.items():
            sarr = np.array([np.square(np.abs(diabats.T.dot(j))) for j in v])
            # print(sarr.T.shape)
            # print(sorted(sarr.T[0]))
            data = np.sum(sarr, 0)
            axes[n].plot(times, data, label=k)
        axes[n].set_title('CSF/SD state vector ' + ('sum' if cmd=='sumcsfv'else '')+'evolution'+title)
        axes[n].set_ylabel('State population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'csfvreim':
        diabats = np.load(os.path.join(basepath, 'zcsf'))
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for CSF/SD ReIm plots!")
            diabats = np.load(os.path.join(basepath, 'csf'))[cmdptx[0]].T
            title = title + f' (for GWP {cmdptx[0]+1})'
        # print(diabats.shape)
        # Expect inp of form label:1,1,0,0
        label, nums = ins.split(':')
        nums = [float(j) for j in nums.split(',')]
        # assert(len(nums)==diabats.shape[1])
        vect = np.array(nums)
        to_plot[label] = vect / np.linalg.norm(vect)
            # Normalize the vector
        reals = np.real(diabats.T.dot(v))
        imgs = np.imag(diabats.T.dot(v))
        axes[n].plot(times, reals, label='Re')
        axes[n].plot(times, imgs, label='Im')

        axes[n].set_title(f'{label} state vector component evolution'+title)
        axes[n].set_ylabel('Coefficent')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'avcsf':
        diabats = np.load(os.path.join(basepath, 'csf_ave'))
        title = ''
        if cmdptx != None:
            if len(cmdptx) > 1: raise Exception("Only one GWP allowed for CSF/SD ReIm plots!")
            diabats = np.load(os.path.join(basepath, 'csf'))[cmdptx[0]].T
            diabats=np.square(np.abs(diabats))
            title = title + f' (for GWP {cmdptx[0]+1})'
        csfs, av_window = ins.split(':')
        av_window = int(av_window)
        CSF_STATES = None if csfs=='A' else [int(i) for i in csfs.split(',')]
        for i in range(diabats.shape[0]):
            if CSF_STATES == None: pass
            else:
                if i+1 not in CSF_STATES: continue
            mav = mathutils.MathUtils.moving_avg(diabats[i], av_window)
            axes[n].plot(times[:len(mav)], mav, label=f'AVCSF {i+1}', color=get_nth_col(i))
        axes[n].set_title(f'{av_window} point moving average CSF state evolution'+title)
        axes[n].set_ylabel('Averaged state population')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'sd':
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

    elif cmd == 'sdla':
        sd = np.load(os.path.join(basepath, 'sdla_ave'))
        SDS = None if ins=='A' else [int(i) for i in ins.split(',')]
        for i in range(len(manifest['spindenmapLA'])):
            atom_number = manifest['spindenmapLA'][i]
            if SDS == None : pass
            elif atom_number not in SDS: continue
            
            try: symbol = ATOMICLABELS[manifest['atomnos'][str(atom_number)]-1]
            except: symbol = '?'
            axes[n].plot(times, sd[i], label='{} [{}]'.format(symbol, atom_number), color=get_nth_col(atom_number))
        axes[n].set_title('Spin density')
        axes[n].set_ylabel('Spin density')
        axes[n].set_xlabel('Time (fs)')
        axes[n].legend(loc='upper right')

    elif cmd == 'avsd':
        sd = np.load(os.path.join(basepath, 'sd_ave'))
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

    elif cmd == 'mq':
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
    elif cmd == 'hm':
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
    elif cmd == 'fft':
        # [FFT=sd|sdla|mq|csf(+ Do Phase Correction | - Do cos correction | 2 Power spectra):Cutoff:1,2,3|A:limxmin,limxmax]
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
        elif mode == 'sdla': data = np.load(os.path.join(basepath, 'sdla_ave'))
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
            
            elif mode in ['sd', 'sdla', 'mq']: # SD/MQ are picked based on atom number
                if mode == 'sd' : atom_number = manifest['spindenmap'][i]
                elif mode == 'sdla' : atom_number = manifest['spindenmapLA'][i]
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
