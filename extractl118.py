from glogpy.l118_job import l118_job
import argparse
import sys, os
import mathutils
import numpy as np
import json
import glob
from glogpy.freqency_job import frequency_job

MANIFEST_NAME='manifest.json'

parser = argparse.ArgumentParser(description='logall extractor')
parser.add_argument('log', type=str,
                    help='Path to L118 log file, or if --glob is set, a glob expression')

parser.add_argument('--adir', type=str, default='analysis',
                    help='Where to put the extracted data files [defaults to LOGDIR/analysis]')

parser.add_argument('--glob', dest='glob', action='store_true',
                    help='Instead of a single file use a glob')

parser.add_argument('--ts', dest='ts',
                    help='Use a constant timestep interpolation')

parser.add_argument('--spindens', dest='spindens', action='store_true',
                    help='Extract spin density')

parser.add_argument('--freq', dest='freq', type=str,
                    help='Path to freq .log file')
args = parser.parse_args()
print(args)

manifest = {}
manifest['method'] = 'l118'
manifest['quantities'] = []

if args.glob:
    TIMESTEP = float(args.ts)
    
    glob_exp = glob.glob(args.log)
    print(f'Starting multi-trajectory glob\nIdentified {len(glob_exp)} candidate files')

    OUTDIR = os.path.join(os.path.dirname(glob_exp[0]), args.adir)
    # os.makedirs(OUTDIR)
    print(f'Starting multi-trajectory glob\nIdentified {len(glob_exp)} candidate files')
    datas = []

    for logfile in glob_exp:

        with open(logfile, 'r') as f:
            data = f.read()
        try:
            data = data.split('Initial command:\n')[-1] # Assume the L118 job is the last (NBOS->CAS->*TRAJ*)
            gj = l118_job(data)
            _, xns = gj.parse(spin_dens=args.spindens)
        except:
            print(f'Failed to parse {logfile}! Aborting extraction...')
            sys.exit(-1)
        datas.append(xns)
    NTRAJ = len(datas)
    print(f'Parsed {NTRAJ} L118 log files OK - starting analysis')

    # PART I : ESTABLISH NUMBER OF TIMESTEPS
    # Use the shortest trajectory as the limit
    maxtime = min([max([i[1][0]['time'] for i in j]) for j in datas])
    NSTEPS = int(np.floor(maxtime / TIMESTEP))
    print(f'TMAX = {maxtime}fs // TS = {TIMESTEP}fs // NSTEPS = {NSTEPS} [= {NSTEPS * TIMESTEP}fs]')

    # PART II : EXTRACT VIBRATIONAL INFO -> Maybe do a check that rotation has been projected out?
    geoms = np.array([[mathutils.MathUtils.dict_to_list(i[1][0]['geom']) for i in j] for j in datas])
    print(geoms)
    
    #  Load up frequency file
    assert(os.path.exists(args.freq))
    with open(args.freq, 'r') as f:
        data = f.read()
    freq = frequency_job(data)
    freq_data = freq.parse()
    print(f'FREQ file parsed')
    # Need to clean up the geom -> compute normal modes
    geom_init = mathutils.MathUtils.dict_to_list(freq_data['geom'])
    geom_init = np.array([x[1] for x in geom_init])
    # Compute nm2xyz and xyz2nm matrices
    nm2xyz, xyz2nm = mathutils.NormModeUtils.nm_matrix(freq_data['atomnos'], freq_data['vibfreqs'], freq_data['vibdisps'])

    nmdata = mathutils.NormModeUtils.xyz_to_nm(xyz2nm, geom_init, geoms)
    nnmode = nmdata.shape[2] # Number of normal modes
    nmdata = nmdata.transpose(1,0,2)

    print(nmdata)
    
    # TODO : INTP NM CSF / CI / SD / NM / XYZ

    sys.exit(-1)

else:
    LOGFILE = args.log
    assert(os.path.exists(LOGFILE))

    with open(LOGFILE, 'r') as f:
        data = f.read()
    data = data.split('Initial command:\n')[-1] # Assume the L118 job is the last
    gj = l118_job(data)

    OUTDIR = args.adir if args.adir[0] == '/' else os.path.join(os.path.dirname(LOGFILE), args.adir)
    os.makedirs(OUTDIR)

    l202, xns = gj.parse(spin_dens=args.spindens)
    manifest['atomnos'] = l202['proton_nums']
    manifest['steps'] = len(xns)

    # save times
    times = np.array([i[1][0]['time'] for i in xns])
    with open(os.path.join(OUTDIR, 'times'), 'wb') as f:
        np.save(f, times)
    manifest['quantities'].append('t')

    # save CI pop
    adiabats = np.abs(np.array([mathutils.MathUtils.dict_to_list(i[0]['adiabats']) for i in xns])).T
    with open(os.path.join(OUTDIR, 'ci_ave'), 'wb') as f:
        np.save(f, adiabats)
    manifest['quantities'].append('ci')

    # save CSF pop
    diabats = np.abs(np.array([mathutils.MathUtils.dict_to_list(i[0]['diabats']) for i in xns])).T
    with open(os.path.join(OUTDIR, 'csf_ave'), 'wb') as f:
        np.save(f, diabats)
    manifest['quantities'].append('csf')

    # save xyz
    xyz = np.array([mathutils.MathUtils.dict_to_list(i[1][0]['geom']) for i in xns])
    with open(os.path.join(OUTDIR, 'xyz_ave'), 'wb') as f:
        np.save(f, xyz)
    manifest['quantities'].append('xyz')

    if args.spindens:
        manifest['spindenmap'] = list(xns[0][2]['spinden_sum'].keys())
        sd_new = np.zeros((manifest['steps'], len(manifest['spindenmap'])))
        for i, d in enumerate([i[2]['spinden_sum'] for i in xns]):

            newlist = [d[j] for j in manifest['spindenmap']]
            sd_new[i] = newlist

        with open(os.path.join(OUTDIR, 'sd_ave'), 'wb') as f:
            np.save(f, sd_new.T)
        manifest['quantities'].append('sd')

    print('All tasks done! Writing manifest...')
    manifest_path = os.path.join(OUTDIR, MANIFEST_NAME)
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
