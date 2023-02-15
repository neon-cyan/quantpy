import glob, sys, json
from glogpy.job import gaussian_job
from glogpy.linkparser import linkparsers
import numpy as np
from mathutils import Stitcher

class casjob(gaussian_job):
    # Define a simple CAS job to get the CI energies + pops
    def __init__(self, txt):
        super().__init__(txt)

    def parse(self):
        l510  = next(filter(lambda x : x.number==510 , self.link_list))
        l405  = next(filter(lambda x : x.number==405 , self.link_list))
        l9999 = next(filter(lambda x : x.number==9999, self.link_list))
        if l510==None or l405==None or l9999==None:
            raise Exception("Invalid job file!")
        else: 
            return linkparsers.L405(l405.text), linkparsers.L510_TD(l510.text, True), l9999.text
# print(glob.glob(sys.argv[1]))

elems = []

if len(sys.argv <2): 
    print(f"Use {sys.argv[0]} /glob/exp/to/logs/*.log")
    sys.exit(-1)

for f in glob.glob(sys.argv[1]):
    # Read in, construct and parse each log file
    with open(f, 'r') as f:
        data = f.read()
    data = casjob(data)
    l405, l510, l9999 = data.parse()
    nbasis = l405['n_basis']
    zeropop = lambda length, d : np.array([d[i+1] if i+1 in d else 0.0 for i in range(length)])
    cic_rc = {i: zeropop(nbasis, j) for i,j in l510['cic'].items()}
    cic = np.array([[cic_rc[i] for i in sorted(cic_rc.keys())]])
    cies = np.array([l510['cie'][i] for i in sorted(l510['cie'].keys())])
    coords = json.loads(''.join(l9999).split("""\\\\""")[2].replace("'", '"'))
    # elems contains a dict of coords, the energies and CICoeffs
    elems.append((coords, cies, cic))

# Clean up and sort by distance to minima to help with stitching
elems = sorted(elems, key=lambda x : np.linalg.norm(np.array([i for _,i in x[0].items()])))
s = Stitcher.compute(np.array([x[2] for x in elems]).transpose(1,2,0,3), quiet=True)
es = Stitcher.run(np.array([[x[1] for x in elems]]).transpose(0,2,1), s)
fixed_order = list(elems[0][0].keys())
nci_states = es.shape[1]

# Output just gets printed for now - pipe into a file if needed
print("{}\t{}".format('\t'.join(fixed_order), '\t'.join([f"CI{x}" for x in range(1, nci_states+1)])))
for coord, energy in zip([x[0] for x in elems],es[0].T):
    for i in fixed_order:
        print('{:>11.8f}'.format(coord[i]),end="\t")
    for c in energy:
        print('{:>11.8f}'.format(c),end="\t")
    print()
    pass