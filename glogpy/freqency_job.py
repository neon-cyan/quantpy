from glogpy.job import gaussian_job
from glogpy.linkparser import linkparsers
import numpy as np

class frequency_job(gaussian_job):
    def __init__(self, txt):
        super().__init__(txt)

    def parse(self):
        # Overall scheme 
        # 1 - Pull atoms from L202
        # 2 - Pull vibrtional data from L716
        l202 = None
        for x in self.link_list:
            if x.number == 202 : 
                l202 = x
                break
                
        if l202 == None : raise Exception("[FREQ] No link 202!")
        ans = linkparsers.L202(l202.text)

        NoSymm = False
        try:
            z = l202.iops[15]
            if z == '1' : NoSymm = True
        except KeyError : pass
        
        geom = ans['geom'] if NoSymm else ans['geom_symm']
        # Sometimes the keys are strored as strings - this is a quick fix
        if 1 not in geom:
            geom = {int(k) : v for k, v in geom.items()}
        proton_info = ans['proton_nums']

        res = None
        for x in self.link_list:
            if x.number == 716 and x.iops[8] == '11': 
                res = linkparsers.L716_hpmodes(x.text, len(geom))
                break

        if res == None: raise Exception("[FREQ] Unable to find a suitable L716")
        res['geom'] = geom
        res['proton_nums'] = proton_info
        res['natom'] = len(geom)
        # print(res)
        return res

# {
#     'vibdisps': array([NATOM x 3]), 
#     'vibfreqs': array([ 415.5995 ... 3231.3848]), 
#     'vibirs': array([ 0.    ... 41.8481]), 
#     'vibsyms': ['E2U' ... 'B2U'],
#     'vibredm': array([2.9171 ... 1.2477]), 
        
#     'atommasses': {1: 1.00783 ... 12: 12.0}, 
#     'atomnos': {1: 1 ... 12: 6}, 
        
#     'temperature': 298.15, 
#     'pressure': 1.0, 
        
#     'forces': {1: array([ 0.        ,  0.00013371, -0.        ]) ... 12: array([-0.00199148,  0.00114978,  0.        ])}, 
#     'maxforce': 0.002299566, 
#     'rmsforce': 0.00094038, 
        
#     'geom': {'1': [0, array([0.      , 2.484212, 0.      ])] ... '12': [0, array([ 1.209657, -0.698396,  0.      ])]}, 
#     'proton_nums': {1: 1, 2: 1, 3: 1, 4: 1, 5: 1, 6: 1, 7: 6, 8: 6, 9: 6, 10: 6, 11: 6, 12: 6}, 
#     'natom': 12
# }