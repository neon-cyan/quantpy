from mathutils import  MathUtils
import os 
import numpy as np

import re
from glogpy.dynamics_job import dynamics_job as dj

class ParseLogAll():
    def parse(txt, step_lim=None, do_CI_States=True):
        sections = re.split("\*{4} Time.{1,}\d{1,}\.\d{1,}.{1,}\*{4}", txt )
        res = []
        nsteps = 0
        for s in sections[1:]:
            if step_lim!=None:
                if nsteps == step_lim: break
            try:
                data = dj(s.strip())
                data = data.parse(do_CI_States=do_CI_States)
                res.append(data)
            except: raise Exception(f'An error occured parsing step {nsteps}')
            nsteps += 1
        return res, nsteps

    def I_ImportLogalls(basedir, ngwps, gwp_dir_fmt='gwp{}_V1', fname='gwp{}_V1_dd_data_nm.logall', 
    step_lim=None, print_steps=True, pickle_if_fail=True,
    quantities=['xyz', 'ci', 'csf', 'mq', 'sd', 'an', 'am']):
        # Quantities options are
        # xyz   = geometries
        # ci    = CI Coeff
        # csf   = CSF Coeff
        # mq    = Mulliken Q         #TODO Split into sum/unsum?
        # sd    = Spin density       #TODO Split into sum/unsum?
        # dp    = Dipole moment
        # an    = Atom numbers
        # am    = Atom masses
        # fo    = Forces
        # maxf  = max + rms force
        # case  = CASSCF Energy
        # casde = CASSCF DE

        steps = None
        datax = []

        for i in range(1, ngwps+1): 
            if print_steps: print(f'Parsing GWP {i}/{ngwps}')
            # Plus one is to compensate for weirdness of python range()
            # We want to loop over 1,2,3 ... ngwp-1, ngwp
        
            loagallf = os.path.join(basedir, gwp_dir_fmt.format(i), fname.format(i))
            assert(os.path.exists(loagallf))
            with open(loagallf, 'r') as f:
                data = f.read()
            parsed_data, nsteps = ParseLogAll.parse(data, step_lim = step_lim, do_CI_States=('ci' in quantities))
            del(data)
            # Make sure all GWP loagalls contain same num steps
            if steps == None:
                steps = nsteps
            else:
                assert(steps == nsteps)
            datax.append(parsed_data)

        assert(len(datax) == ngwps) # Quick sanity check (GWP in = GWP out)

        results = {}
        results['steps'] = steps
        # First pull constants
        if 'an' in quantities:
            results['atomnos'] = datax[0][0]['atomnos']
        if 'am' in quantities:
            results['atommasses']  = datax[0][0]['atommasses']
        # CI data (state pops, energies and coefs)
        if 'ci' in quantities:
            results['adiabats'] = np.zeros([ngwps, steps, len(datax[0][0]['adiabats'])], dtype=complex)
            results['cic'] = np.zeros([ngwps, steps, len(datax[0][0]['cie']), len(datax[0][0]['cic'][1])])    
            results['cies'] = np.zeros([ngwps, steps, len(datax[0][0]['cie'])])

        #  Work way through scalar data to gen [GWP x Step]

        if 'csf' in quantities:
            results['diabats'] =  np.zeros([ngwps, steps, len(datax[0][0]['diabats'])],   dtype=complex)
        if 'maxf' in quantities:
            results['maxf'] = np.zeros([ngwps, steps])
            results['rmsf'] = np.zeros([ngwps, steps])

        # Vector quantities [GWP x Step x Dim]
        if 'xyz' in quantities:
            results['geomx'] = np.zeros([ngwps, steps, len(datax[0][0]['geom_init']), 3])

        if 'fo' in quantities:
            results['forces'] = np.zeros((ngwps, steps, len(datax[0][0]['forces']), 3))
            
        if 'dp' in quantities:
            results['dipolemom'] = np.zeros((ngwps, steps, 3))

        if 'casde' in quantities:
            results['casde'] = np.zeros((ngwps, steps))

        if 'case' in quantities:
            results['case'] = np.zeros((ngwps, steps))

        if 'mq' in quantities:
            results['mullikensum'] = np.zeros([ngwps,steps,  len(datax[0][0]['mulliken_sum'])])
            results['mullikenmap'] = [int(i) for i in list(datax[0][0]['mulliken_sum'].keys())]
        
        if 'sd' in quantities:
            results['spindensum'] = np.zeros([ngwps,steps,  len(datax[0][0]['spinden_sum'])])
            results['spindenmap'] = [int(i) for i in list(datax[0][0]['spinden_sum'].keys())]

        for i, gwp in enumerate(datax):
            for j, ts in enumerate(gwp):
                try:
                    if 'ci' in quantities:
                        results['adiabats'][i,j] = np.array(MathUtils.dict_to_list(ts['adiabats']))
                    if 'ci' in quantities:
                        results['cies'][i,j] = np.array(MathUtils.dict_to_list(ts['cie']))
                        results['cic'][i,j] = np.array([MathUtils.dict_to_list(x) for x in MathUtils.dict_to_list(ts['cic'])])

                    if 'csf' in quantities:
                        results['diabats'][i,j] = np.array(MathUtils.dict_to_list(ts['diabats']))
                    
                    if 'xyz' in quantities:
                        gtemp = MathUtils.dict_to_list(ts['geom_init'])
                        gtemp = [x[1] for x in gtemp]
                        results['geomx'][i,j] = np.array(gtemp)
                    
                    if 'fo' in quantities:
                        ftemp =  MathUtils.dict_to_list(ts['forces'])
                        results['forces'][i,j] = np.array(ftemp)

                    if 'dp' in quantities:
                        results['dipolemom'][i,j] = np.array(ts['dipole'][0])

                    if 'casde' in quantities:
                        results['casde'][i,j] = ts['casde']

                    if 'case' in quantities:
                        results['case'][i,j] = ts['case']

                    if 'mq' in quantities:
                        for atomidx, mullsum in ts['mulliken_sum'].items():
                            results['mullikensum'][i,j,results['mullikenmap'].index(atomidx)] = mullsum

                    if 'sd' in quantities:
                        for atomidx, sdsum in ts['spinden_sum'].items():
                            results['spindensum'][i,j,results['spindenmap'].index(atomidx)] = sdsum

                    if 'maxf' in quantities:
                        results['maxf'][i,j] = ts['maxforce']
                        results['rmsf'][i,j] = ts['rmsforce']
                except:
                    print(f'An error occured in reformating LOGALL files at GWP{i+1} / TS {j+1}!')
                    if pickle_if_fail:
                        import pickle
                        picklepath = os.path.join(basedir, 'data')
                        print(f'SAVING TO {picklepath}')
                        with open(picklepath, 'wb') as f:
                            pickle.dump(datax, f)
                    raise Exception('Failed to parse - halting')
        del(datax)
        return results