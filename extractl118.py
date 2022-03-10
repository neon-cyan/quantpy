from glogpy.l118_job import l118_job
import argparse
import sys, os
import mathutils
import numpy as np
import json

MANIFEST_NAME='manifest.json'

parser = argparse.ArgumentParser(description='logall extractor')
parser.add_argument('log', type=str,
                    help='Path to L118 log file, or if --glob is set, a glob expression')

parser.add_argument('--adir', type=str, default='analysis',
                    help='Where to put the extracted data files [defaults to LOGDIR/analysis]')

parser.add_argument('--glob', dest='glob', action='store_true',
                    help='Instead of a single file use a ')

parser.add_argument('--spindens', dest='spindens', action='store_true',
                    help='Extract spin density')

args = parser.parse_args()
print(args)

manifest = {}
manifest['method'] = 'l118'
manifest['quantities'] = []

if args.glob:
    OUTDIR = args.adir
    print('Globbing multiple logs is NYI')
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
