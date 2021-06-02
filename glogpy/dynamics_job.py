from os import EX_CANTCREAT
from glogpy.job import gaussian_job
from glogpy.linkparser import linkparsers
import numpy as np

#  Jobs of the form
#  1/xxxx/1,18;
#  2/xxxx/2;          <- GEOM INIT
#  3/xxxx/1,2,3;
#  4/xxxx/1,5;
#  5/xxxx/10;         <- State composition
#  8/xxxx/1;
#  11/xxx/1;
#  10/xxx/3(-3);
#  6/xxxx/1;          <- Spin densiry / mulliken
#  7/xxxx/1,2,3,16;   <- Initial forces
#  1/xxxx/18(3);
#  2/xxxx/2;          <- GEOM FINAL
#  99/xxx/99;

class dynamics_job(gaussian_job):
    def __init__(self, txt):
        super().__init__(txt)

    def parse(self):
        l202_init = None
        for x in self.link_list:
            if x.number == 202 : 
                l202_init = x
                break
                
        if l202_init == None : raise Exception("[DYNX] No link 202!")
        geom = linkparsers.L202(l202_init.text)['geom']

        res = None
        for x in self.link_list[::-1]: # Look for the last L510 -> State composition
            if x.number == 510 and 97 in x.iops.keys():
                res = linkparsers.L510_TD(x.text)
                break
        if res==None : raise Exception("[DYNX] No dynamics-compatible L510 found")

        res2 = None
        for x in self.link_list[::-1]: # Look for the last L601 -> Miliken pop + sd
            if x.number == 601:
                res2 = linkparsers.L601(x.text, spin_dens=True, dipole=True)
                break
        if res2==None : raise Exception("[DYNX] No dynamics-compatible L601 found")
        res.update(res2)

        res2 = None
        for x in self.link_list[::-1]: # Look for the last L716 -> MaxForce + Forces
            if x.number == 716:
                res2 = linkparsers.L716_hpmodes(x.text, len(geom), False)
                break
        if res2==None : raise Exception("[DYNX] No dynamics-compatible L716 found")
        res.update(res2)

        res2 = None
        for x in self.link_list[::-1]: # Look for the last L202
            if x.number == 202:
                if x == l202_init : raise Exception("[DYNX] No second L202 found!")
                res2 = linkparsers.L202(x.text)
                break

        res['geom_final'] = res2['geom']

        res['geom_init'] = geom
        # print(res)
        return res

# {
# 'diabats': {1: (-0.22095+0.01429j), 2: (0.11627-0.00537j), 3: (0.16042+0.00451j), 4: (0.49084-0.1658j), 5: (0.24327+0.2859j), 6: (-0.05558+0.04399j), 7: (0.49918-0.09541j), 8: (-0.41637+0.25625j)}, 
# 'adiabats': {1: (0.73058-0.15666j), 2: (-0.02168+0.18205j), 3: (0.35406-0.05912j), 4: (-0.2881+0.02187j), 5: (-0.04753+0.18635j), 6: (-0.15433+0.09799j), 7: (-0.20731-0.05285j), 8: (-0.04674+0.27822j)}, 
# 
# 'muliken': {1: 0.860896, 2: -0.471127, 3: -0.333442, 4: -0.328013, 5: -0.732424, 6: 0.582688, 7: 0.315301, 8: 0.315301, 9: 0.39541, 10: 0.39541}, 
# 'mulliken_sum': {1: 0.860896, 2: -0.471127, 3: 0.249246, 4: 0.302588, 5: 0.058397}, 
# 
# 'spinden': {1: 0.208977, 2: 0.250161, 3: 0.452028, 4: 0.042886, 5: 0.045869, 6: 0.001223, 7: -6e-05, 8: -6e-05, 9: -0.000513, 10: -0.000513}, 
# 'spinden_sum': {1: 0.208977, 2: 0.250161, 3: 0.453252, 4: 0.042767, 5: 0.044843}, 
#
# 'dipole': [['-3.5701', '4.1027', '-0.0000'], '5.4385'],
#
# 'atommasses': {1: 12.0, 2: 15.99491, 3: 15.99491, 4: 12.0, 5: 14.00307, 6: 1.00783, 7: 1.00783, 8: 1.00783, 9: 1.00783, 10: 1.00783}, 
# 
# 'atomnos': {1: 6, 2: 8, 3: 8, 4: 6, 5: 7, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1}, 
# 'temperature': 298.15, 
# 'pressure': 1.0, 
#
# 'forces': {1: array([ 1.65471143e-01, -1.61748502e-01,  6.00000000e-08]), 2: array([-2.6836681e-02,  6.0881460e-03, -9.3000000e-08]), 3: array([-9.57545490e-02,  1.74179419e-01,  2.13000000e-07]), 4: array([ 6.999268e-03, -2.780492e-02, -2.400000e-08]), 5: array([ 1.2235064e-02,  8.0989678e-02, -5.0000000e-09]), 6: array([-4.3020526e-02, -2.9154888e-02, -1.1200000e-07]), 7: array([ 0.00334231, -0.00568359, -0.00953219]), 8: array([ 0.0033423 , -0.00568355,  0.00953219]), 9: array([-0.01288919, -0.01559086, -0.01395954]), 10: array([-0.01288915, -0.01559093,  0.01395951])}, 
# 'maxforce': 0.174179419, 
# 'rmsforce': 0.059306177, 
#
# 'geom_final': {1: [0, array([-0.090002,  0.512569,  0.      ])], 2: [0, array([1.342   , 0.766244, 0.      ])], 3: [0, array([-1.050758,  1.451089,  0.      ])], 4: [0, array([-0.641986, -0.894045,  0.      ])], 5: [0, array([ 0.39821 , -1.987024,  0.      ])], 6: [0, array([-0.87376 ,  2.454532,  0.      ])], 7: [0, array([-1.280347, -0.952623,  0.893902])], 8: [0, array([-1.280347, -0.952623, -0.893902])], 9: [0, array([ 0.997514, -1.794707,  0.80825 ])], 10: [0, array([ 0.997514, -1.794707, -0.80825 ])]}, 
# 'geom_init': {1: [0, array([-0.077164,  0.554472,  0.      ])], 2: [0, array([1.354838, 0.808147, 0.      ])], 3: [0, array([-1.03792 ,  1.492992,  0.      ])], 4: [0, array([-0.629148, -0.852142,  0.      ])], 5: [0, array([ 0.411048, -1.945121,  0.      ])], 6: [0, array([-0.860922,  2.496435,  0.      ])], 7: [0, array([-1.267509, -0.91072 ,  0.893902])], 8: [0, array([-1.267509, -0.91072 , -0.893902])], 9: [0, array([ 1.010352, -1.752804,  0.80825 ])], 10: [0, array([ 1.010352, -1.752804, -0.80825 ])]}
# }