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


# print(a)
# Now the confiogs are loaded we can begin constructing the matrices !
# First we define 2nd q opertaors a_{p*} and a^+_{p*}
def ann(s, pos, spin):
    if s == 0: return 0
    ans = copy.deepcopy(s)
    if s[pos] == 'a' and spin == 'a':   ans = ans[:pos] + '0' + ans[pos+1:]
    elif s[pos] == 'b' and spin == 'b': ans = ans[:pos] + '0' + ans[pos+1:]
    elif s[pos] == '1' and spin == 'a': ans = ans[:pos] + 'b' + ans[pos+1:]
    elif s[pos] == '1' and spin == 'b': ans = ans[:pos] + 'a' + ans[pos+1:]
    else: return 0
    return ans

def crn(s, pos, spin):
    if s == 0: return 0
    ans = copy.deepcopy(s)
    if s[pos] == 'a' and spin == 'b':   ans = ans[:pos] + '1' + ans[pos+1:]
    elif s[pos] == 'b' and spin == 'a': ans = ans[:pos] + '1' + ans[pos+1:]
    elif s[pos] == '0' and spin == 'a': ans = ans[:pos] + 'a' + ans[pos+1:]
    elif s[pos] == '0' and spin == 'b': ans = ans[:pos] + 'b' + ans[pos+1:]
    else: return 0
    return ans
# Next we need S_z - trivial for slaterdets
def sz(s):
    tot = 0
    for i in s:
        if i == 'a': tot += 1
        elif i == 'b': tot -=1
    return tot/2

# Now we construct the matrices we will populate
n_configs, norbs = len(a), len(a[0])
print(f'About to process {n_configs} configs of {norbs} orbitals')
z = np.zeros([n_configs, n_configs])
spm = np.zeros([n_configs, n_configs])

for ni, i in enumerate(a): # Loop over config
    assert(len(i) == norbs)
    z[ni,ni] = sz(i)       # Populate S_z matrix
    for j in range(norbs):
        sms = crn(ann(i, j, 'a'),j,'b')      # S-
        for k in range(norbs):
            sps = crn(ann(sms,k,'b' ),k,'a') # S+S-
#            print(i, ni, j, sps)
            if sps != 0 : 
                if sps == i:
                    spm[a.index(sps), ni] = 1.0
                else:
                    spm[a.index(sps), ni] = -1.0 # Not sure why off-diag elems are negative

# Now just need to do the diagonalisation
ssq = spm + np.matmul(z, z-np.identity(n_configs))
evals, evects = np.linalg.eig(ssq)
evects = evects.T
print('Diag OK')

# Check + discard imag parts
assert(max(np.abs(np.imag(evals))) < 1e-12)
evals, evects = np.real(evals), np.real(evects)
evals = np.round(evals, 2)

# Print out the values
def cleanprint(x, cutoff):
    ans = ''
    for idx, i in enumerate(x):
        if np.abs(i) > cutoff:
            ans += f'{i:.3f} ({idx+1})  '
    return ans

for ev, evv in sorted(zip(evals, evects), key=lambda x : np.abs(x[0])):
    print(f'S**2 = {ev} >>> {cleanprint(evv, 1e-3)}')