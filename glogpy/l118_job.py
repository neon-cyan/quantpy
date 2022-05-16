from glogpy.job import gaussian_job
from glogpy.linkparser import linkparsers
import numpy as np

#  1/xxx/1,18;
#  2/xxx/2;         <- GEOM INIT
#  3/xxx/1,2,3;
#  4/xxx/1,5;
#  5/xxx/10;        <- L510 MCSCF INIT
#  8/xxx/1;
#  11/xxx/1;
#  10/xxx/3(-3);
#  6/xxx/1;         <- SD/MQ L601 INIT
#  7/xxx/1,2,3,16;
#  1/xxx/18(3);
#  2/xxx/2;
#  7/xxx/16;
#  99/xxx/99;
#  2/xxx/2;          
#  3/xxx/1,2,3;
#  4/xxx/1,5;
#  5/xxx/10;         <- L510s
#  8/xxx/1;
#  11/xxx/1;
#  10/xxx/3(-3);
#  6/xxx/1;          <- SD/MQ L601s
#  7/xxx/1,2,3,16;
#  1/xxx/18(-9);     <- GEOM + Time
#  2/xxx/2;
#  6/xxx/1;
#  7/xxx/16;
#  99/xxx/99;

# This is a definition for a single traj Eh L118 job
# TODO test with GradOnly fixed timestep

class l118_job(gaussian_job):
    def __init__(self, txt):
        super().__init__(txt)

    def parse(self, print_info=True, spin_dens=False,do_CI_States=True):
        l202_init = None
        for x in self.link_list:
            if x.number == 202 : 
                l202_init = x
                break
         
        if l202_init == None : raise Exception("[L118P] No link 202!")
        l202_init = linkparsers.L202(l202_init.text) # Need this for atom labels

        l510s = filter(lambda x: x.number==510, self.link_list)
        l510_init, *l510s = [linkparsers.L510_TD(x.text, do_CI_States=do_CI_States) for x in list(l510s)]

        l118s = filter(lambda x: x.number==118, self.link_list)
        l118s = [linkparsers.L118(x.text) for x in list(l118s)[1:-1]]
        l118_init = l118s[0][0]
        # print(len(l118s[0]))
        try:    l118s[0] = [l118s[0][1]]
        except IndexError: 
            l118s[0] = [l118s[0][0]]
            print("You have missed the first step - perhaps you specified init velocity?")

        l601s = filter(lambda x: x.number==601, self.link_list)
        l601_init, *l601s = [linkparsers.L601(x.text, spin_dens) for x in list(l601s)[:-1]]
        assert(len(l510s) == len(l118s) == len(l601s)) # The number of in-loop links must be consistent
        xn = [(l510_init, [l118_init], l601_init)] + list(zip(l510s, l118s, l601s)) # Prepend the first step

        return (l202_init, xn)