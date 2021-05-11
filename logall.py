from mathutils import  MathUtils
import os 
import numpy as np

import re
from glogpy.dynamics_job import dynamics_job as dj

class ParseLogAll():
    def parse(txt, step_lim=None):
        sections = re.split("\*{4} Time.{1,}\d{1,}\.\d{1,}.{1,}\*{4}", txt )
        res = []
        nsteps = 0
        for s in sections[1:]:
            if step_lim!=None:
                if nsteps == step_lim: break

            data = dj(s.strip())
            data = data.parse()
            res.append(data)
            nsteps += 1
        return res, nsteps

def I_ImportLogalls(basedir, ngwps, gwp_dir='gwp{}_V1', fname='gwp{}_V1_dd_data_nm.logall', step_lim=None):
    stepchecker = None
    datax = []

    for i in range(1, ngwps+1): 
        # Plus one is to compensate for weirdness of python range()
        # We want to loop over 1,2,3 ... ngwp-1, ngwp
    
        loagallf = os.path.join(basedir, gwp_dir.format(i), fname.format(i))
        assert(os.path.exists(loagallf))
        with open(loagallf, 'r') as f:
            data = f.read()
        parsed_data, nsteps = ParseLogAll.parse(data, step_lim = step_lim)

        # Make sure all GWP loagalls contain same num steps
        if stepchecker == None:
            stepchecker = nsteps
        else:
            assert(stepchecker == nsteps)
        datax.append(parsed_data)

    assert(len(datax) == ngwps) # Quick sanity check

    # Pull out constant params - atom numbers / atom weights / forces @ t=0
    atomnos     = datax[0][0]['atomnos']
    atommasses  = datax[0][0]['atommasses']
    init_forces = datax[0][0]['forces']
    #  Work way through the data to gen GWPx matrices
    diabats =  np.zeros([ngwps, stepchecker, len(datax[0][0]['diabats'])],   dtype=complex)
    adiabats = np.zeros([ngwps, stepchecker, len(datax[0][0]['adiabats'])],  dtype=complex)
    geomx    = np.zeros([ngwps, stepchecker, len(datax[0][0]['geom_init']), 3])

    mullikensum = np.zeros([ngwps,stepchecker,  len(datax[0][0]['mulliken_sum'])])
    mullikenmap = list(datax[0][0]['mulliken_sum'].keys())

    spindensum = np.zeros([ngwps,stepchecker,  len(datax[0][0]['spinden_sum'])])
    spindenmap = list(datax[0][0]['spinden_sum'].keys())

    biggest_maxf = np.zeros(ngwps)

    for i, gwp in enumerate(datax):
        for j, ts in enumerate(gwp):
            diabats[i,j] = np.array(MathUtils.dict_to_list(ts['diabats']))
            adiabats[i,j] = np.array(MathUtils.dict_to_list(ts['adiabats']))

            gtemp = MathUtils.dict_to_list(ts['geom_init'])
            gtemp = [x[1] for x in gtemp]
            geomx[i,j] = np.array(gtemp)

            for atomidx, mullsum in ts['mulliken_sum'].items():
                mullikensum[i,j,mullikenmap.index(atomidx)] = mullsum

            for atomidx, sdsum in ts['spinden_sum'].items():
                spindensum[i,j,spindenmap.index(atomidx)] = sdsum

        biggest_maxf[i] = max([x['maxforce'] for x in gwp])

    res = {    
        'atomnos'       : atomnos,
        'atommasses'    : atommasses,
        'init_forces'   : init_forces,
        'diabats'       : diabats,
        'adiabats'      : adiabats,
        'geomx'         : geomx,
        'mullikensum'   : mullikensum,
        'mullikenmap'   : mullikenmap,
        'spindensum'    : spindensum,
        'spindenmap'    : spindenmap,
        'biggest_maxf'  : biggest_maxf,
        'steps'         : stepchecker,
    }
    return res