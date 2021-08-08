import argparse
import os
import json
import mathutils
from glogpy.freqency_job import frequency_job
from quatics_lexers import QuanticsParsers
from logall import ParseLogAll
import numpy as np

MANIFEST_NAME='manifest.json'

parser = argparse.ArgumentParser(description='logall extractor')
parser.add_argument('qinp', type=str,
                    help='Path to quantics input file')

parser.add_argument('--adir', type=str, default='analysis',
                    help='Where to put the extracted data files [defaults to QINPDIR/analysis]')

parser.add_argument('--nsteps', default='A',
                    help='Number of steps to parse [defaults to all steps]')


#TODO : maxforce
#TODO : energy
#TODO : Dipole
#TODO : CAS Delta-E


#TODO MANIFEST

parser.add_argument('--tasks', type=str, nargs='*', default=['nm', 'ci', 'csf', 'mq', 'sd', 'xyz'],
                    help='''Things to extract. Options are
                    nm   = normal modes
                    xyz  = average geometeries 
                    ci   = adiabatic states
                    csf  = diabatic states
                    mq   = mulliken charges (H summed into HAs)
                    sd   = spin density
                    dp   = dipole [WIP]
                    fo   = forces
                    maxf = maxforce [WIP]
                    [defaults to nm, ci, csf, mq, sd, xyz]''')

parser.add_argument('--redo', dest='redo', action='store_true',
                    help='Ignore the extraction manifest (if it exists) [defaults to false]')

parser.add_argument('--quiet', dest='quiet', action='store_true',
                    help='Silence the output [WIP]')

args = parser.parse_args()
print(args)

QINP = args.qinp
STEPLIMS = int(args.nsteps) if args.nsteps != 'A' else None

# Quick sanity check that the quantics input is real
assert(os.path.exists(QINP))

# Prepare the output dir + load in manifest, if available
OUTDIR = args.adir if args.adir[0] == '/' else os.path.join(os.path.dirname(QINP), args.adir)
manifest_path = os.path.join(OUTDIR, MANIFEST_NAME)
manifest = None
if os.path.exists(OUTDIR):
    print('Output path exists')
    if os.path.exists(manifest_path):
        if args.redo:
            print('Ignore existing manifest')
        else:
            print('Attempting to load existig manifest')
            with open(manifest_path, 'r') as f:
                manifest = json.load(f)
    else:
        print('No manifest has been found')
else: os.makedirs(OUTDIR)

# Parse the quantics input
with open(QINP, 'r') as f:
    data = f.read()
q_inp_data=QuanticsParsers.parse_input(data)
print(f'QINP file parsed')

# Quick check that the output path exists
datadir = os.path.join(os.path.dirname(QINP), q_inp_data['data'])
assert(os.path.exists(datadir))

# Load up quantics output file & clean up needed properties
qoutf = os.path.join(os.path.dirname(QINP), q_inp_data['name'], 'output')
assert(os.path.exists(qoutf))
with open(qoutf, 'r') as f:
    data = f.read()
q_out_data=QuanticsParsers.parse_output(data)
if STEPLIMS!=None: q_out_data = q_out_data[:STEPLIMS] # Truncate the quantics log if necessary
times = np.array([x['time'] for x in q_out_data])
gwp_sf = np.array([x['GGP'] for x in q_out_data])
gwp_sf /= 10 # Sum up to 10 - not 1
print(f'QOUT file parsed')


# Load in all GWP logalls
tasks = ['an', 'am'] + args.tasks
data_gwpx = ParseLogAll.I_ImportLogalls(datadir, q_inp_data['ngwp'], step_lim=STEPLIMS, quantities=tasks)
nsteps = data_gwpx['steps']
assert(nsteps == len(q_out_data))

# Check against manifest
if manifest != None: 
    try: assert(nsteps == manifest['steps'])
    except: raise Exception('Conflicting step numbers in script/manifest! Run with --redo to rebuild')
else: manifest = {'steps' : nsteps}

print(f'LOGALL files parsed')


if 'nm' in args.tasks or 'xyz' in args.tasks:
    #  Load up frequency file
    freqf = os.path.join(datadir, q_inp_data['freqf'])
    assert(os.path.exists(freqf))
    with open(freqf, 'r') as f:
        data = f.read()
    freq = frequency_job(data)
    freq_data = freq.parse()
    print(f'FREQ file parsed')

    # Need to clean up the geom -> compute normal modes
    geom_init = mathutils.MathUtils.dict_to_list(freq_data['geom'])
    geom_init = np.array([x[1] for x in geom_init])

    # Compute nm2xtyz and xyz2nm matrices
    nm2xyz, xyz2nm = mathutils.NormModeUtils.nm_matrix(data_gwpx['atommasses'], 
    freq_data['vibfreqs'], freq_data['vibdisps'])

    nmdata = mathutils.NormModeUtils.xyz_to_nm(xyz2nm, geom_init, data_gwpx['geomx'])
    nnmode = nmdata.shape[2] # Number of normal modes
    nmdata = nmdata.transpose(1,0,2)

    # I HAVE NO IDEA HOW TO USE NP.DOT // PLEASE NO BULLY
    res = np.zeros((nsteps, freq_data['vibdisps'].shape[0])) # NM x S matrix
    for i in range(res.shape[0]):
        res[i] = gwp_sf[i].dot(nmdata[i])
    res = res.T

    avegeom = np.zeros((nsteps, geom_init.shape[0], 3)) # Average atom geometries (at each step)
    for i in range(nsteps):
        avegeom[i] = np.copy(geom_init)
        for j in range(nnmode):
            avegeom[i] += res[j, i] * nm2xyz.T[j].reshape(geom_init.shape[0], 3)

manifest['atomnos'] = data_gwpx['atomnos']
manifest['atommasses'] = data_gwpx['atommasses']

manifest['quantities'] = []
# Export of weights
with open(os.path.join(OUTDIR, 'weights'), 'wb') as f:
    np.save(f, gwp_sf)
manifest['quantities'].append('w')

# Export times
with open(os.path.join(OUTDIR, 'times'), 'wb') as f:
    np.save(f, times)
manifest['quantities'].append('t')


def weightscale_sq(datax, weight, ns):
    datax_ave = np.zeros([datax.shape[0], ns]) # [NDIABATS x NSTEPS]
    for i in range(datax_ave.shape[0]):
        for x in range(nsteps):
            datax_ave[i,x] = np.square(abs(datax[i, x]).dot(weight[x]))
    return datax_ave

def weightscale(datax, weight, ns):
    datax_ave = np.zeros([datax.shape[0], ns]) # [NDIABATS x NSTEPS]
    for i in range(datax_ave.shape[0]):
        for x in range(nsteps):
            datax_ave[i,x] = abs(datax[i, x]).dot(weight[x])
    return datax_ave


for task in args.tasks:
    if task=='ci':
        adiabats = data_gwpx['adiabats'].transpose(2,1,0)
        adiabats_ave = weightscale_sq(adiabats, gwp_sf, nsteps)
        with open(os.path.join(OUTDIR, 'ci'), 'wb') as f:
            np.save(f, adiabats)
        with open(os.path.join(OUTDIR, 'ci_ave'), 'wb') as f:
            np.save(f, adiabats_ave)
        manifest['quantities'].append('ci')

    elif task=='csf':
        diabats = data_gwpx['diabats'].transpose(2,1,0)
        diabats_ave = weightscale_sq(diabats, gwp_sf, nsteps)
        with open(os.path.join(OUTDIR, 'csf'), 'wb') as f:
            np.save(f, diabats)
        with open(os.path.join(OUTDIR, 'csf_ave'), 'wb') as f:
            np.save(f, diabats_ave)
        manifest['quantities'].append('csf')

    elif task=='mq':
        msum = data_gwpx['mullikensum'].transpose(2,1,0)
        msum_ave = weightscale(msum, gwp_sf, nsteps)
        with open(os.path.join(OUTDIR, 'mq'), 'wb') as f:
            np.save(f, msum)
        with open(os.path.join(OUTDIR, 'mq_ave'), 'wb') as f:
            np.save(f, msum_ave)
        manifest['quantities'].append('mq')
        manifest['mullikenmap'] = data_gwpx['mullikenmap']

    elif task=='sd':
        sdsum = data_gwpx['spindensum'].transpose(2,1,0)
        sdsum_ave = weightscale(sdsum, gwp_sf, nsteps)
        with open(os.path.join(OUTDIR, 'sd'), 'wb') as f:
            np.save(f, sdsum)
        with open(os.path.join(OUTDIR, 'sd_ave'), 'wb') as f:
            np.save(f, sdsum_ave)
        manifest['quantities'].append('sd')
        manifest['spindenmap'] = data_gwpx['spindenmap']

    # Nortmal mode & ave geom have seperate logic (see above)
    elif task=='nm':
        with open(os.path.join(OUTDIR, 'nm_ave'), 'wb') as f:
            np.save(f, res)
        with open(os.path.join(OUTDIR, 'nm'), 'wb') as f:
            np.save(f, nmdata)
        manifest['quantities'].append('nm')

    elif task=='xyz':
        with open(os.path.join(OUTDIR, 'xyz'), 'wb') as f:
            np.save(f, geom_init)
        with open(os.path.join(OUTDIR, 'xyz_ave'), 'wb') as f:
            np.save(f, avegeom)
        manifest['quantities'].append('xyz')

    elif task=='dp':
        # Currently scalar treatment -> vectors may need some fixing re-orentiation
        dps = np.array([[np.linalg.norm(x) for x in y] for y in data_gwpx['dipolemom']])
        dps = dps.T
        dps_ave = np.array([dps[i].dot(gwp_sf[i]) for i in range(dps.shape[0])])
        with open(os.path.join(OUTDIR, 'dps'), 'wb') as f:
            np.save(f, dps)
        with open(os.path.join(OUTDIR, 'dps_ave'), 'wb') as f:
            np.save(f, dps_ave)
        manifest['quantities'].append('dp')

    else:
        print(f'Unknown task {task}')
        # This is a warning for now - in next verson[s] this will throw an Exception

print('All tasks done! Writing manifest...')
with open(manifest_path, 'w') as f:
    json.dump(manifest, f, indent=2)