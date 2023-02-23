from glogpy.job import gaussian_job
from glogpy.linkparser import linkparsers
import numpy as np
import sys, copy

# This script tries to build all legal meta-CSFs via diagonalising the S**2 matrix
# To do this you need to supply either a valid L405-containing log file OR feed in the configs manually

# The first part of the code deals with reading in the configurations

class configlistjob(gaussian_job):
    def __init__(self, txt):
        super().__init__(txt)

    def parse(self):
        res = None
        for x in self.link_list[::-1]: # Look for the last L405
            if x.number == 405:
                res = linkparsers.L405(x.text, configs=True)
                break
        if res==None : raise Exception("No compatible L405 found")
        return res

a = []
if sys.argv[1] == 'MANUAL':
    while True:
        i = input()
        if i == '': break
        try:
            a.append(i.split()[-1])
        except:
            sys.exit(-1)
else:
    with open(sys.argv[1], 'r') as f:
        data = f.read()
    data = data.split('Initial command:\n')[-1]
    job = configlistjob(data)
    # job.prettyprint()
    data = job.parse()
    assert(data['slater'])
    a = [i.split()[-1] for i in data['config_list']]

def fullflip(s):
    ans = copy.deepcopy(s)
    ans = list(ans)
    for idx, i in enumerate(s):
        if i == 'a':
            ans[idx] = 'b'
        elif i == 'b':
            ans[idx] = 'a'
    return ans

n_configs, norbs = len(a), len(a[0])
print(f'About to process {n_configs} configs of {norbs} orbitals')

def posnlisttostr(length, posarray):
    ans = np.zeros(length)
    for pos, val in posarray:
        ans[pos] = val
    return ','.join([f'{i:.3f}' for i in ans])


seen = []
sings = []
trips = []
for i in a:
    complement = ''.join(fullflip(i))
    if i in seen or complement in seen:
        continue

    if i == complement:
        # i is self-complement (i.e. closed shell det e.g. 1100)
        sings.append(posnlisttostr(n_configs, [(a.index(i), 1.0)]))
    else:
        sings.append(posnlisttostr(n_configs, [(a.index(i), 1.0), (a.index(complement), 1.0)]))
        trips.append(posnlisttostr(n_configs, [(a.index(i), 1.0), (a.index(complement), -1.0)]))

    seen.append(i)

print('Singlet:{}_Triplet:{}'.format('+'.join(sings), '+'.join(trips)))
print(f'S:{len(sings)} T:{len(trips)}')