from mathutils import  MathUtils
import os 
import numpy as np

import re
from glogpy.dynamics_job import dynamics_job as dj

class ParseLogAll():
    def ImportLogalls(basedir, ngwps, step_lim, gwp_dir_fmt='gwp{}_V1', fname='gwp{}_V1_dd_data_nm.logall', 
    quiet=False):
        first = True
        results = {}
        results['steps'] = step_lim
        for i in range(ngwps): 
            loagallf = os.path.join(basedir, gwp_dir_fmt.format(i+1), fname.format(i+1))
            assert(os.path.exists(loagallf))
            with open(loagallf, 'r') as f:
                txt = f.read()
            if not quiet: print(f'Parsing GWP {i+1}/{ngwps}')
            sections = re.split("\*{4} Time.{1,}\d{1,}\.\d{1,}.{1,}\*{4}", txt )
            nsteps = 0
            for j, s in enumerate(sections[1:]):
                if nsteps == step_lim: break

                try:
                    log = dj(s.strip())
                    data = log.parse(do_CI_States=True)
                except: raise Exception(f'An error occured parsing step {nsteps} in GWP {i}!')

                if first:
                    # First step : figure out avaialble quantities & construct arrays
                    # print(data)
                    # print([(x.number, x.iops) for x in log.link_list])

                    l510s = list(filter(lambda x: x.number == 510, log.link_list))
                    assert(len(l510s) == 1) # Only one L510 allowed per log
                    do_spindens = False if 72 not in l510s[0].iops else int(l510s[0].iops[72]) != 0
                    nbasis =  data['n_basis']

                    # First pull constants
                    results['atomnos'] = data['atomnos']
                    results['atommasses']  = data['atommasses']

                    # CI data (state pops, energies and coefs)
                    results['adiabats'] = np.zeros([ngwps, step_lim, len(data['adiabats'])], dtype=complex)
                    results['cic'] = np.zeros([ngwps, step_lim, len(data['cie']), data['n_basis']])  
                    results['cies'] = np.zeros([ngwps, step_lim, len(data['cie'])])

                    #  Work way through scalar data to gen [GWP x Step]
                    results['diabats'] =  np.zeros([ngwps, step_lim, len(data['diabats'])],   dtype=complex)
                    results['maxf'] = np.zeros([ngwps, step_lim])
                    results['rmsf'] = np.zeros([ngwps, step_lim])

                    # Vector quantities [GWP x Step x Dim]
                    results['geomx'] = np.zeros([ngwps, step_lim, len(data['geom_init']), 3])
                    results['forces'] = np.zeros((ngwps, step_lim, len(data['forces']), 3))
                    results['dipolemom'] = np.zeros((ngwps, step_lim, 3))
                    results['casde'] = np.zeros((ngwps, step_lim))
                    results['case'] = np.zeros((ngwps, step_lim))

                    results['mullikensum'] = np.zeros([ngwps, step_lim,  len(data['mulliken_sum'])])
                    results['mullikenmap'] = [int(i) for i in list(data['mulliken_sum'].keys())]

                    if do_spindens:
                        results['spindensum'] = np.zeros([ngwps, step_lim,  len(data['spinden_sum'])])
                        results['spindenmap'] = [int(i) for i in list(data['spinden_sum'].keys())]
                        results['spindenLA'] = np.zeros([ngwps, step_lim,  len(data['spinden'])])
                        results['spindenmapLA'] = [int(i) for i in list(data['spinden'].keys())]

                    first = False

                # Now populate arrays
                results['adiabats'][i,j] = np.array(MathUtils.dict_to_list(data['adiabats']))
                results['cies'][i,j] = np.array(MathUtils.dict_to_list(data['cie']))
                for nstate, ci in enumerate(MathUtils.dict_to_list(data['cic'])):
                    # print(ci)
                    for s in range(nbasis):
                        try:
                            results['cic'][i,j,nstate,s] = ci[s+1]
                        except KeyError:
                            results['cic'][i,j,nstate,s] = 0.0
                    # print(results['cic'][i,j,nstate])

                # print(data['diabats'], MathUtils.dict_to_list(data['diabats']) )

                results['diabats'][i,j] = np.array(MathUtils.dict_to_list(data['diabats']))
            
                gtemp = MathUtils.dict_to_list(data['geom_init'])
                gtemp = [x[1] for x in gtemp]
                results['geomx'][i,j] = np.array(gtemp)
                ftemp =  MathUtils.dict_to_list(data['forces'])
                results['forces'][i,j] = np.array(ftemp)
                results['dipolemom'][i,j] = np.array(data['dipole'][0])
                results['casde'][i,j] = data['casde']
                results['case'][i,j] = data['case']

                for atomidx, mullsum in data['mulliken_sum'].items():
                    results['mullikensum'][i,j,results['mullikenmap'].index(atomidx)] = mullsum

                results['maxf'][i,j] = data['maxforce']
                results['rmsf'][i,j] = data['rmsforce']                                

                if do_spindens:
                    for atomidx, sdsum in data['spinden_sum'].items():
                        results['spindensum'][i,j,results['spindenmap'].index(atomidx)] = sdsum
                    for atomidx, sd in data['spinden'].items():
                        results['spindenLA'][i,j,results['spindenmapLA'].index(atomidx)] = sd

                nsteps += 1
        return results