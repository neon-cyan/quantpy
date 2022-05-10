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

    def parse(self, do_CI_States=False, spin_dens=True):
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
                res = linkparsers.L510_TD(x.text, do_CI_States=do_CI_States)
                break
        if res==None : raise Exception("[DYNX] No dynamics-compatible L510 found")

        res2 = None
        for x in self.link_list[::-1]: # Look for the last L601 -> Miliken pop + sd
            if x.number == 601:
                res2 = linkparsers.L601(x.text, spin_dens=spin_dens, dipole=True)
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
