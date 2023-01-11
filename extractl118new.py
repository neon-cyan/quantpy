from glogpy.l118_traj_job import l118_job
from glogpy.freqency_job import frequency_job
from mathutils import Stitcher
import argparse
import sys, os
import mathutils
import numpy as np
import json

MANIFEST_NAME='manifest.json'

parser = argparse.ArgumentParser(description='L118 Eh job extractor')
parser.add_argument('log', type=str,
                    help='Path to L118 log file, or if --glob is set, a glob expression')

parser.add_argument('-a','--adir', type=str, default='analysis',
                    help='Where to put the extracted data files [defaults to LOGDIR/analysis]')

parser.add_argument('--glob', dest='glob', action='store_true',
                    help='Instead of a single file use a glob for multiple trajectories [NYI]')

parser.add_argument('-s','--stitch', dest='stitch', action='store_true',
                    help='Stitch the CI pops & energies to maintain state charecter throughout')

parser.add_argument('-t','--truncate', dest='allow_truncate', action='store_true',
                    help='Allow for mismatched route links (Allows parsing of partial logs)')      

parser.add_argument('-f','--freq', dest='freq',
                    help='Specify a frequency file for normal mode analysis')

args = parser.parse_args()
# print(args)

manifest = {}
manifest['method'] = 'l118'
manifest['quantities'] = []

if args.glob:
    OUTDIR = args.adir
    print('Globbing multiple logs is NYI')
    sys.exit(-1)

else:

    if args.freq != None:
        print("Reading FREQ file...")
        with open(args.freq, 'r') as f:
            data = f.read()
        data = data.split('Initial command:\n')[-1] # Assume the freq job is the last
        freq = frequency_job(data)
        freq_data = freq.parse()
        # print(freq_data)
        print(f'FREQ file read OK')

    LOGFILE = args.log
    assert(os.path.exists(LOGFILE))
    gj = l118_job(LOGFILE, allow_partial=args.allow_truncate)

    OUTDIR = args.adir if args.adir[0] == '/' else os.path.join(os.path.dirname(LOGFILE), f'{LOGFILE}.{args.adir}')
    try: os.makedirs(OUTDIR)
    except FileExistsError: pass

    l202, xns , l405= gj.parse(allow_truncate=args.allow_truncate)
    # print(l405) => Looks like {'slater': True, 'n_basis': 10, 'ndel': 0}
    manifest['atomnos'] = l202['proton_nums']
    manifest['steps'] = len(xns)

    # save times
    times = np.array([i[1]['time'] for i in xns])
    with open(os.path.join(OUTDIR, 'times'), 'wb') as f:
        np.save(f, times)
    manifest['quantities'].append('t')

    # CI state compositions, energies, populations & stitches
    ci_composition = np.array([mathutils.MathUtils.dict_to_list(i[0]['cic']) for i in xns])    
    # print(f'CICOMPSHAPES = {ci_composition.shape}')
    pop = lambda x, n: np.array([x[i] if i in x else 0.0 for i in range(n)])
    # print(l405)
    ci_composition = np.array([[[pop(j, l405['n_basis']) for j in i] for i in ci_composition]]).transpose((0,2,1,3))
    # print(f'CICOMPSHAPE = {ci_composition.shape}')

    cie_energies = np.array([[mathutils.MathUtils.dict_to_list(i[0]['cie']) for i in xns]]).transpose((0,2,1))
    # print('CIESSHPE=', cie_energies.shape)

    adiabats = np.array([[mathutils.MathUtils.dict_to_list(i[0]['adiabats']) for i in xns]])    
    adiabats = adiabats[:,:,:cie_energies.shape[1]].transpose(0,2,1)
    # print(adiabats)
    # print('ADBTS = ', adiabats.shape)
    if args.stitch:
        stitches = Stitcher.compute(ci_composition)
        ci_composition = Stitcher.run(ci_composition, stitches)
        cie_energies = Stitcher.run(cie_energies,stitches)
        adiabats = Stitcher.run(adiabats,stitches)
        manifest['stitched'] = True

    with open(os.path.join(OUTDIR, 'cies'), 'wb') as f:
        np.save(f, [cie_energies[0].T])
    manifest['quantities'].append('cies')
    # print(ci_composition.shape)
    with open(os.path.join(OUTDIR, 'cicomp'), 'wb') as f:
        np.save(f, ci_composition)
    manifest['quantities'].append('cicomp')

    with open(os.path.join(OUTDIR, 'ci_ave'), 'wb') as f:
        np.save(f, np.abs(adiabats[0]))
    with open(os.path.join(OUTDIR, 'zci'), 'wb') as f:
        np.save(f, adiabats[0])
    manifest['quantities'].append('ci')

    # save CSF pop
    diabats = np.abs(np.array([mathutils.MathUtils.dict_to_list(i[0]['diabats']) for i in xns])).T
    with open(os.path.join(OUTDIR, 'csf_ave'), 'wb') as f:
        np.save(f, diabats)
    with open(os.path.join(OUTDIR, 'zcsf'), 'wb') as f:
        np.save(f, np.array([np.array([mathutils.MathUtils.dict_to_list(i[0]['diabats']) for i in xns])])[0].T)
    manifest['quantities'].append('csf')

    # save xyz
    xyz = np.array([mathutils.MathUtils.dict_to_list(i[1]['geom']) for i in xns])
    with open(os.path.join(OUTDIR, 'xyz_ave'), 'wb') as f:
        np.save(f, xyz)
    manifest['quantities'].append('xyz')

    nucde = [i[1]['etot']for i in xns]
    nucke = [i[1]['ekin']for i in xns]
    nucpe = [i[1]['epot']for i in xns]
    with open(os.path.join(OUTDIR, 'nucde'), 'wb') as f:
        np.save(f, nucde)
    with open(os.path.join(OUTDIR, 'nucke'), 'wb') as f:
        np.save(f, nucke)
    with open(os.path.join(OUTDIR, 'nucpe'), 'wb') as f:
        np.save(f, nucpe)    
    manifest['quantities'].append('nucde')

    # print(xns[0])

    if args.freq != None:
        nm2xyz, xyz2nm = mathutils.NormModeUtils.nm_matrix(freq_data['atommasses'], freq_data['vibfreqs'], freq_data['vibdisps'])
        geom_init = mathutils.MathUtils.dict_to_list(freq_data['geom'])
        geom_init = np.array([x[1] for x in geom_init])
        # print(xyz, xyz.shape)
        nmdata = mathutils.NormModeUtils.xyz_to_nm(xyz2nm, geom_init, np.array([xyz]))
        # print(nmdata, nmdata.shape)
        with open(os.path.join(OUTDIR, 'nm_ave'), 'wb') as f:
            np.save(f, nmdata[0].T)
        manifest['quantities'].append('nm')


    if 'mulliken_sum' in xns[0][2]:
        manifest['mqmap'] = list(xns[0][2]['mulliken_sum'].keys())
        mq_new = np.zeros((manifest['steps'], len(manifest['mqmap'])))
        for i, d in enumerate([i[2]['mulliken_sum'] for i in xns]):

            newlist = [d[j] for j in manifest['mqmap']]
            mq_new[i] = newlist

        with open(os.path.join(OUTDIR, 'mq_ave'), 'wb') as f:
            np.save(f, mq_new.T)
        manifest['quantities'].append('mq')

    if 'spinden_sum' in xns[0][2]:
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
